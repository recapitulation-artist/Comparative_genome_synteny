#!/usr/bin/env python3
"""Prepare curated or automatically selected genomic neighborhoods from source data."""

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile
from urllib.parse import unquote

from extract_regions import (ROOT, Feature, annotation_indices, descendants,
                             fasta_bytes, fasta_records, genome_file, open_text,
                             prepare_region, require, resolve_bedtools,
                             run_bedtools, table, write_table)


@dataclass
class Annotation:
    features: list
    children: dict
    genes: list
    traces: dict


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_annotation(path):
    """Retain source columns and attributes as well as the parent gene graph."""
    features, traces = [], {}
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == 9, f"Invalid GFF line {number}: {path}")
            attrs = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = unquote(value)
            feature = Feature(fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6], attrs)
            require(1 <= feature.start <= feature.end, f"Invalid GFF interval: {path}:{number}")
            features.append(feature)
            traces[id(feature)] = dict(contig=feature.contig, source=fields[1],
                                       feature_type=feature.kind, start_1based=feature.start,
                                       end_1based=feature.end, score=fields[5], strand=feature.strand,
                                       phase=fields[7], attributes=attrs)
    _, children = annotation_indices(features)
    genes = [feature for feature in features if feature.kind in {"gene", "pseudogene"}]
    identifiers = [gene.identifier for gene in genes]
    require(all(identifiers) and len(identifiers) == len(set(identifiers)),
            f"Missing or duplicate parent gene IDs: {path}")
    return Annotation(features, children, genes, traces)


def truth(value):
    return str(value).lower() in {"true", "yes", "1"}


def coding_gene(gene, children):
    """Count parent genes once, never CDS fragments, RNAs or pseudogenes."""
    nested = descendants(gene, children)
    if gene.kind != "gene" or "pseudo" in gene.attrs.get("gene_biotype", "").lower():
        return False
    if any(truth(feature.attrs.get("pseudo", "")) or feature.attrs.get("pseudogene")
           for feature in nested):
        return False
    biotype = gene.attrs.get("gene_biotype", "").lower()
    if biotype and biotype != "protein_coding":
        return False
    return any(feature.kind == "CDS" for feature in nested)


def match_target(gene, children, target):
    nested = descendants(gene, children)
    names = {feature.attrs.get("gene", "").casefold() for feature in nested}
    products = [feature.attrs.get("product", "") for feature in nested]
    query = target.casefold()
    if query == "nramp":
        reasons = (["gene=mntH"] if "mnth" in names else [])
        if any("nramp" in product.casefold() for product in products):
            reasons.append("product contains Nramp")
    else:
        reasons = ([f"gene={target}"] if query in names else [])
        pattern = re.compile(r"(?<!\w)" + re.escape(target) + r"(?!\w)", re.IGNORECASE)
        if any(pattern.search(product) for product in products):
            reasons.append(f"product contains token {target}")
    return "; ".join(reasons)


def select_neighbors(genes, target, children, count):
    """Select up to count coding parents on either genomic side of a target."""
    require(count >= 0, "Neighbor count must be nonnegative")
    require(coding_gene(target, children), f"Target is not protein-coding: {target.identifier}")
    ordered = sorted((gene for gene in genes if gene.contig == target.contig
                      and coding_gene(gene, children)),
                     key=lambda gene: (gene.start, gene.end, gene.identifier))
    require(any(gene is target for gene in ordered), "Target missing from source gene list")
    index = next(i for i, gene in enumerate(ordered) if gene is target)
    return ordered[max(0, index - count):index + count + 1], index, len(ordered) - index - 1


def compact_json(value):
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def gene_fields(gene, annotation, seq_id, assembly, start0, end0, target_gene, selected=None):
    nested = descendants(gene, annotation.children)
    proteins = sorted({feature.attrs["protein_id"] for feature in nested if feature.attrs.get("protein_id")})
    products = sorted({feature.attrs["product"] for feature in nested if feature.attrs.get("product")})
    local_start0, local_end0 = max(gene.start - 1, start0) - start0, min(gene.end, end0) - start0
    require(0 <= local_start0 < local_end0 <= end0 - start0, f"Invalid local interval: {seq_id}/{gene.identifier}")
    trace = annotation.traces[id(gene)]
    return dict(seq_id=seq_id, assembly_accession=assembly, contig=gene.contig,
                gene_id=gene.identifier, locus_tag=gene.attrs.get("locus_tag", ""),
                feature_type=gene.kind, gene=gene.attrs.get("gene", ""),
                biotype=gene.attrs.get("gene_biotype", ""), start_1based=gene.start,
                end_1based=gene.end, strand=gene.strand, local_start0=local_start0,
                local_end0=local_end0, local_start_1based=local_start0 + 1,
                local_end_1based=local_end0, clipped_left=gene.start - 1 < start0,
                clipped_right=gene.end > end0, is_target=gene is target_gene,
                protein_ids=";".join(proteins), products=";".join(products),
                is_protein_coding=coding_gene(gene, annotation.children),
                selected_for_plot=selected is not None and gene.identifier in selected,
                gff_source=trace["source"], gff_score=trace["score"], gff_phase=trace["phase"],
                gff_attributes=compact_json(gene.attrs),
                gff_features=compact_json([annotation.traces[id(feature)] for feature in nested]))


def extract_region(binary, output, work, seq_id, genome, annotation, contigs,
                   target, target_gene, start0, end0, header, selected=None, legacy=None):
    assembly = genome["assembly_accession"]
    require(target.contig in contigs, f"Missing source contig: {assembly}/{target.contig}")
    dna = contigs[target.contig]
    require(0 <= start0 < end0 <= len(dna), f"Region outside source contig: {seq_id}")
    bed = output / "regions" / f"{seq_id}.bed"
    bed.write_text(f"{target.contig}\t{start0}\t{end0}\t{header}\t0\t+\n", encoding="utf-8")
    fasta, extracted = work / "contig.fna", work / "region.fasta"
    fasta.write_bytes(fasta_bytes([(target.contig, dna)]))
    fasta.with_suffix(".fna.fai").unlink(missing_ok=True)
    run_bedtools(binary, "getfasta", "-fi", fasta, "-bed", bed, "-fo", extracted, "-nameOnly")
    records = list(fasta_records(extracted))
    require(records == [(header, dna[start0:end0])], f"BEDTools sequence differs from source interval: {seq_id}")
    genes_bed = output / "regions" / f"{seq_id}.genes.bed"
    genes_bed.write_text("".join(f"{gene.contig}\t{gene.start - 1}\t{gene.end}\t{gene.identifier}\t0\t{gene.strand}\n"
                                 for gene in annotation.genes), encoding="utf-8")
    overlaps = run_bedtools(binary, "intersect", "-a", genes_bed, "-b", bed, "-wa", "-wb")
    gene_index = {gene.identifier: gene for gene in annotation.genes}
    nearby = []
    for line in overlaps.splitlines():
        fields = line.split("\t")
        require(len(fields) == 12, f"Unexpected BEDTools intersection: {seq_id}")
        nearby.append(gene_fields(gene_index[fields[3]], annotation, seq_id, assembly,
                                  start0, end0, target_gene, selected))
    nearby.sort(key=lambda row: (row["start_1based"], row["end_1based"], row["gene_id"]))
    require(sum(row["is_target"] for row in nearby) == 1, f"Expected one selected target parent: {seq_id}")
    proteins = sorted({feature.attrs["protein_id"] for feature in descendants(target_gene, annotation.children)
                       if feature.attrs.get("protein_id")})
    protein = target.attrs.get("protein_id", "") or (proteins[0] if proteins else "")
    products = sorted({feature.attrs["product"] for feature in descendants(target_gene, annotation.children)
                       if feature.attrs.get("product")})
    region = dict(seq_id=seq_id, fasta_header=header, assembly_accession=assembly,
                  contig=target.contig, contig_length=len(dna), region_start0=start0,
                  region_end0=end0, region_start_1based=start0 + 1, region_end_1based=end0,
                  region_strand="+", length=end0 - start0, target_protein_id=protein,
                  target_start_1based=target.start, target_end_1based=target.end,
                  target_strand=target.strand, target_locus_tag=target_gene.attrs.get("locus_tag", ""),
                  target_product=target.attrs.get("product", "") or ";".join(products),
                  target_local_start_1based=target.start - start0, target_local_end_1based=target.end - start0,
                  legacy_target_start=legacy["start"] if legacy else "",
                  legacy_target_end=legacy["end"] if legacy else "",
                  legacy_annotations_verified=legacy["count"] if legacy else "",
                  assembly_status=genome.get("assembly_status", ""), genome_notes=genome.get("genome_notes", ""))
    return records[0], region, nearby


def write_sequences(path, sequences):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["seq_id", "length"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(sequences)


def source_name(path):
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def load_inputs(args):
    metadata = table(args.genome_metadata)
    required = {"assembly_accession", "organism_name"}
    if args.entry == "curated":
        required.add("seq_id")
    require(metadata and required <= set(metadata[0]),
            "Genome metadata needs " + ", ".join(sorted(required)))
    require(all(all(row[key] for key in required) for row in metadata),
            "Genome metadata contains empty required fields")
    if args.entry == "curated":
        require(len({row["seq_id"] for row in metadata}) == len(metadata), "Duplicate metadata sequence IDs")
    else:
        assembly_metadata = {}
        for row in metadata:
            # Historical region-specific fields are unused by automatic selection.
            assembly_fields = {key: value for key, value in row.items()
                               if key not in {"seq_id", "region_length", "sequence_id_match"}}
            assembly = row["assembly_accession"]
            require(assembly not in assembly_metadata or assembly_metadata[assembly] == assembly_fields,
                    f"Conflicting metadata for assembly: {assembly}")
            assembly_metadata[assembly] = assembly_fields
    manual = table(args.annotations) if args.entry == "curated" else None
    if manual is not None:
        needed = {"seq_id", "start", "end", "strand", "type", "feat_id", "geom_id", "note", "organisms"}
        require(manual and needed <= set(manual[0]), "Curated annotation columns are incomplete")
        require({row["seq_id"] for row in manual} <= {row["seq_id"] for row in metadata},
                "Curated annotations contain sequence IDs missing from genome metadata")
        metadata = [row for row in metadata if any(gene["seq_id"] == row["seq_id"] for gene in manual)]
    sources = [args.genome_metadata, Path(__file__).resolve(), ROOT / "scripts/extract_regions.py"]
    if manual is not None:
        sources.append(args.annotations)
    assemblies = list(dict.fromkeys(row["assembly_accession"] for row in metadata))
    if args.entry == "auto":
        available = {path.name for path in args.genome_dir.iterdir() if path.is_dir()
                     and ((path / "annotation.gff").is_file() or (path / "annotation.gff.gz").is_file())}
        require(available == set(assemblies), "Every genome directory with a GFF must have genome metadata, "
                f"and every indexed genome must exist (unindexed: {sorted(available - set(assemblies))}; "
                f"missing: {sorted(set(assemblies) - available)})")
    annotation_paths, fasta_paths = {}, {}
    for assembly in assemblies:
        directory = args.genome_dir / assembly
        annotation_paths[assembly] = genome_file(directory, "annotation.gff")
        fasta_paths[assembly] = genome_file(directory, "genome.fna")
        sources.extend([annotation_paths[assembly], fasta_paths[assembly]])
    hashes = [dict(file=source_name(path), sha256=sha256(path)) for path in sources]
    return metadata, manual, annotation_paths, fasta_paths, hashes


def prepare(args, binary, version, output, inputs):
    metadata, manual, annotation_paths, fasta_paths, hashes = inputs
    (output / "regions").mkdir(parents=True, exist_ok=True)
    annotations, contigs = {}, {}

    def load(assembly):
        if assembly not in annotations:
            annotations[assembly] = load_annotation(annotation_paths[assembly])
            records = list(fasta_records(fasta_paths[assembly]))
            require(len(records) == len({name for name, _ in records}), f"Duplicate FASTA contig: {assembly}")
            contigs[assembly] = dict(records)
        return annotations[assembly], contigs[assembly]

    records, regions, neighbors, genes, genomes, sequences, candidates = [], [], [], [], [], [], []
    with tempfile.TemporaryDirectory(prefix=".prepare-", dir=output) as temporary:
        work = Path(temporary)
        if args.entry == "curated":
            manual_groups = defaultdict(list)
            for row in manual:
                manual_groups[row["seq_id"]].append(row)
            source_fields = {}
            for genome in metadata:
                seq_id, assembly = genome["seq_id"], genome["assembly_accession"]
                group = manual_groups[seq_id]
                sequence = dict(seq_id=seq_id, length=max(int(row["end"]) for row in group))
                genome = dict(genome, region_length=sequence["length"])
                annotation, dna = load(assembly)
                _, target, old_target, start0, end0, _ = prepare_region(sequence, group, genome, annotation.features)
                matches = [gene for gene in annotation.genes
                           if any(child is target for child in descendants(gene, annotation.children))]
                require(len(matches) == 1, f"Curated target has no unique parent gene: {seq_id}")
                selected_ids = set()
                lookup, _ = annotation_indices(annotation.features)
                for row in group:
                    feature = lookup[(row["type"], row["feat_id"])][0]
                    parents = [gene for gene in annotation.genes if gene is feature or
                               any(child is feature for child in descendants(gene, annotation.children))]
                    require(len(parents) == 1, f"Curated feature has no unique parent gene: {seq_id}/{row['feat_id']}")
                    parent = parents[0]
                    selected_ids.add(parent.identifier)
                    source_fields[(seq_id, row["type"], row["feat_id"])] = parent.identifier
                result = extract_region(binary, output, work, seq_id, genome, annotation, dna, target,
                                        matches[0], start0, end0, f"{target.contig}_{old_target['feat_id']}",
                                        selected_ids, dict(old_target, count=len(group)))
                record, region, nearby = result
                records.append(record)
                regions.append(region)
                neighbors.extend(nearby)
                sequences.append(sequence)
                genomes.append(genome)
            nearby_index = {(row["seq_id"], row["gene_id"]): row for row in neighbors}
            for row in manual:
                original = nearby_index[(row["seq_id"], source_fields[(row["seq_id"], row["type"], row["feat_id"])])]
                enriched = dict(row)
                for key, value in original.items():
                    if key == "seq_id":
                        continue
                    enriched[key if key not in enriched else "source_" + key] = value
                enriched["is_target"] = row["note"] == "Uptake-Gene"
                genes.append(enriched)
        else:
            assembly_index = {}
            for row in metadata:
                assembly_index.setdefault(row["assembly_accession"], row)
            for assembly, genome in assembly_index.items():
                annotation, dna = load(assembly)
                targets = []
                for gene in sorted(annotation.genes, key=lambda feature: (feature.contig, feature.start,
                                                                          feature.end, feature.identifier)):
                    reason = match_target(gene, annotation.children, args.target_gene)
                    if not reason:
                        continue
                    eligible = coding_gene(gene, annotation.children)
                    nested = descendants(gene, annotation.children)
                    products = sorted({feature.attrs["product"] for feature in nested if feature.attrs.get("product")})
                    proteins = sorted({feature.attrs["protein_id"] for feature in nested if feature.attrs.get("protein_id")})
                    candidates.append(dict(assembly_accession=assembly, contig=gene.contig,
                                           gene_id=gene.identifier, locus_tag=gene.attrs.get("locus_tag", ""),
                                           gene=gene.attrs.get("gene", ""), products=";".join(products),
                                           protein_ids=";".join(proteins), start_1based=gene.start,
                                           end_1based=gene.end, strand=gene.strand,
                                           match_reason=reason, selected=eligible,
                                           exclusion_reason="" if eligible else "RNA, pseudogene or missing CDS",
                                           seq_id=""))
                    if eligible:
                        targets.append((gene, candidates[-1]))
                for rank, (target_gene, candidate) in enumerate(targets, 1):
                    seq_id = f"{assembly.replace('.', '_')}_r{rank:02d}"
                    require(len(seq_id) <= 25, f"Region sequence ID exceeds 25 characters: {seq_id}")
                    candidate["seq_id"] = seq_id
                    selected, before, after = select_neighbors(annotation.genes, target_gene,
                                                               annotation.children, args.neighbors)
                    selected_ids = {gene.identifier for gene in selected}
                    start0, end0 = min(gene.start for gene in selected) - 1, max(gene.end for gene in selected)
                    record, region, nearby = extract_region(binary, output, work, seq_id, genome,
                                                            annotation, dna, target_gene, target_gene,
                                                            start0, end0, seq_id, selected_ids)
                    region.update(target_gene_id=target_gene.identifier,
                                  requested_neighbors_per_side=args.neighbors,
                                  upstream_genes_selected=min(before, args.neighbors),
                                  downstream_genes_selected=min(after, args.neighbors),
                                  target_query=args.target_gene)
                    genus = genome["organism_name"].split()[0]
                    for row in nearby:
                        if row["selected_for_plot"]:
                            note = row["gene"] or row["products"] or "Other"
                            genes.append(dict(row, start=row["local_start_1based"], end=row["local_end_1based"],
                                              note=note, organisms=genus))
                    records.append(record)
                    regions.append(region)
                    neighbors.extend(nearby)
                    sequences.append(dict(seq_id=seq_id, length=end0 - start0))
                    genomes.append(dict(genome, seq_id=seq_id, region_length=end0 - start0))
            if candidates:
                write_table(output / "target_matches.tsv", candidates)
            else:
                write_table(output / "target_matches.tsv", [dict(assembly_accession="", contig="", gene_id="",
                            locus_tag="", gene="", products="", protein_ids="", start_1based="", end_1based="",
                            strand="", match_reason="No matching original GFF annotation", selected=False,
                            exclusion_reason="No target matches", seq_id="")])
            require(records, f"No protein-coding target matches for {args.target_gene}; "
                    f"see {output / 'target_matches.tsv'}. Any previous scientific outputs in this "
                    "directory are not results for the requested target.")

    require(len({row["seq_id"] for row in regions}) == len(regions), "Duplicate extracted region IDs")
    require(sum(row["is_target"] for row in genes) == len(regions), "Expected one target marker per region")
    current_ids = {row["seq_id"] for row in regions}
    expected_beds = {f"{identifier}{suffix}" for identifier in current_ids for suffix in (".bed", ".genes.bed")}
    for path in (output / "regions").glob("*.bed"):
        if path.name not in expected_beds:
            path.unlink()
    fasta = output / "extracted_sequences.fasta"
    fasta.write_bytes(fasta_bytes(records)[:-1])
    metrics = dict(regions=len(regions), plotted_genes=len(genes), neighbor_genes=len(neighbors),
                   target_genes=sum(row["is_target"] for row in genes),
                   clipped_left_genes=sum(row["clipped_left"] for row in neighbors),
                   clipped_right_genes=sum(row["clipped_right"] for row in neighbors),
                   negative_target_strands=sum(row["target_strand"] == "-" for row in regions),
                   total_bases=sum(len(dna) for _, dna in records), fasta_sha256=sha256(fasta))
    if args.entry == "auto":
        metrics.update(target_query=args.target_gene, requested_neighbors_per_side=args.neighbors,
                       target_candidate_genes=len(candidates))
    for name, rows in (("genomes.tsv", genomes), ("genes.tsv", genes), ("regions.tsv", regions),
                       ("neighbors.tsv", neighbors), ("input_hashes.tsv", hashes),
                       ("extraction_metrics.tsv", [dict(metric=key, value=value) for key, value in metrics.items()]),
                       ("extraction_versions.tsv", [dict(tool="bedtools", version=version),
                                                     dict(tool="python", version=sys.version.split()[0])])):
        write_table(output / name, rows)
    write_sequences(output / "sequences.csv", sequences)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entry", choices=("curated", "auto"), required=True)
    parser.add_argument("--genome-dir", type=Path, default=ROOT / "data/genomes")
    parser.add_argument("--genome-metadata", type=Path, default=ROOT / "data/genomes.tsv")
    parser.add_argument("--annotations", type=Path, default=ROOT / "data/genes.tsv")
    parser.add_argument("--target-gene", default="Nramp")
    parser.add_argument("--neighbors", type=int, default=10)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bedtools", help="BEDTools 2.31.0 executable")
    parser.add_argument("--verify-only", action="store_true", help="Regenerate independently and compare existing outputs")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    require(args.neighbors >= 0, "Neighbor count must be nonnegative")
    args.target_gene = args.target_gene.strip()
    require(bool(args.target_gene), "Target gene must not be empty")
    resolve = lambda path: (path if path.is_absolute() else ROOT / path).resolve()
    for key in ("genome_dir", "genome_metadata", "annotations"):
        setattr(args, key, resolve(getattr(args, key)))
    output = resolve(args.output) if args.output else ROOT / "results" / args.entry
    require(output != ROOT and ROOT / "data" not in (output, *output.parents)
            and args.genome_dir not in (output, *output.parents), "Output must not overwrite source inputs")
    produced = {"sequences.csv", "genomes.tsv", "genes.tsv", "regions.tsv", "neighbors.tsv",
                "extracted_sequences.fasta", "input_hashes.tsv", "extraction_metrics.tsv",
                "extraction_versions.tsv", "target_matches.tsv"}
    protected = [args.genome_metadata] + ([args.annotations] if args.entry == "curated" else [])
    require(all(source not in {output / name for name in produced} for source in protected),
            "Output would overwrite a supplied metadata or annotation file")
    binary, version = resolve_bedtools(args.bedtools)
    inputs = load_inputs(args)
    if args.verify_only:
        require(output.is_dir(), f"Missing prepared output: {output}")
        versions = table(output / "extraction_versions.tsv")
        require(any(row["tool"] == "bedtools" and row["version"] == version for row in versions),
                "BEDTools version metadata differs")
        with tempfile.TemporaryDirectory(prefix="entry-prepare-check-") as temporary:
            regenerated = Path(temporary)
            metrics = prepare(args, binary, version, regenerated, inputs)
            for expected in regenerated.rglob("*"):
                if not expected.is_file() or expected.name == "extraction_versions.tsv":
                    continue
                actual = output / expected.relative_to(regenerated)
                require(actual.is_file() and actual.read_bytes() == expected.read_bytes(),
                        f"Prepared output differs from source regeneration: {actual}")
            require({path.name for path in (output / "regions").glob("*.bed")}
                    == {path.name for path in (regenerated / "regions").glob("*.bed")},
                    "Unexpected stale region BED files")
        print(f"Verified {args.entry} preparation against source genomes and annotations.")
    else:
        metrics = prepare(args, binary, version, output, inputs)
    print(f"Prepared {metrics['regions']} {args.entry} regions, {metrics['plotted_genes']} plotted genes "
          f"and {metrics['target_genes']} targets.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"Entry preparation failed: {error}")
