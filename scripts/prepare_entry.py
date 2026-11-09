#!/usr/bin/env python3
"""Prepare curated or automatically selected genomic neighborhoods from source data."""

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
from urllib.parse import unquote

from extract_regions import (ROOT, Feature, annotation_indices, descendants,
                             fasta_bytes, fasta_records, genome_file, open_text,
                             prepare_region, require, resolve_bedtools,
                             run_bedtools, table, write_table)
from runtime import DEFAULT_THREADS, execution_context, validate_threads


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
    return repository_input(path).relative_to(ROOT).as_posix()


def repository_input(path):
    """Reject external inputs, including files reached through external symlinks."""
    path = Path(path)
    resolved = (path if path.is_absolute() else ROOT / path).resolve()
    require(resolved.is_relative_to(ROOT), f"Input must be inside this repository: {path}")
    return resolved


def protect_region_inputs(output, paths):
    """The published region directory must never contain any source input."""
    directory = (Path(output) / "regions").resolve()
    for path in paths:
        source = Path(path).resolve()
        require(not source.is_relative_to(directory),
                f"Output regions would overwrite a supplied analysis input: {source}")


def load_inputs(args):
    args.genome_dir = repository_input(args.genome_dir)
    args.genome_metadata = repository_input(args.genome_metadata)
    if args.entry == "curated":
        args.annotations = repository_input(args.annotations)
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
        require(all(re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", row["seq_id"]) for row in metadata),
                "Sequence IDs must be safe filenames using letters, digits, underscore, dot or hyphen")
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
    sources = [args.genome_metadata, Path(__file__).resolve(), ROOT / "scripts/extract_regions.py",
               ROOT / "scripts/runtime.py"]
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
        annotation_paths[assembly] = repository_input(genome_file(directory, "annotation.gff"))
        fasta_paths[assembly] = repository_input(genome_file(directory, "genome.fna"))
        sources.extend([annotation_paths[assembly], fasta_paths[assembly]])
    hashes = [dict(file=source_name(path), sha256=sha256(path)) for path in sources]
    return metadata, manual, annotation_paths, fasta_paths, hashes


def worker_environment():
    """Each BEDTools worker uses one CPU from the shared process budget."""
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[name] = "1"


def assembly_worker(task):
    """Load one assembly once and return only plain records, never GFF objects."""
    started = time.monotonic()
    assembly = task["assembly"]
    output = Path(task["output"])
    (output / "regions").mkdir(parents=True)
    annotation = load_annotation(Path(task["annotation"]))
    fasta = list(fasta_records(Path(task["fasta"])))
    require(len(fasta) == len({name for name, _ in fasta}), f"Duplicate FASTA contig: {assembly}")
    dna = dict(fasta)
    del fasta
    results, candidates, source_fields = [], [], {}
    with tempfile.TemporaryDirectory(prefix="work-", dir=output) as temporary:
        work = Path(temporary)
        if task["entry"] == "curated":
            lookup, _ = annotation_indices(annotation.features)
            for genome in task["metadata"]:
                seq_id = genome["seq_id"]
                group = task["manual_groups"][seq_id]
                sequence = dict(seq_id=seq_id, length=max(int(row["end"]) for row in group))
                genome = dict(genome, region_length=sequence["length"])
                _, target, old_target, start0, end0, _ = prepare_region(
                    sequence, group, genome, annotation.features)
                matches = [gene for gene in annotation.genes
                           if any(child is target for child in descendants(gene, annotation.children))]
                require(len(matches) == 1, f"Curated target has no unique parent gene: {seq_id}")
                selected_ids = set()
                for row in group:
                    feature = lookup[(row["type"], row["feat_id"])][0]
                    parents = [gene for gene in annotation.genes if gene is feature or
                               any(child is feature for child in descendants(gene, annotation.children))]
                    require(len(parents) == 1,
                            f"Curated feature has no unique parent gene: {seq_id}/{row['feat_id']}")
                    parent = parents[0]
                    selected_ids.add(parent.identifier)
                    source_fields[(seq_id, row["type"], row["feat_id"])] = parent.identifier
                record, region, nearby = extract_region(
                    task["binary"], output, work, seq_id, genome, annotation, dna, target,
                    matches[0], start0, end0, f"{target.contig}_{old_target['feat_id']}",
                    selected_ids, dict(old_target, count=len(group)))
                results.append(dict(record=record, region=region, neighbors=nearby,
                                    sequence=sequence, genome=genome, genes=[]))
        else:
            genome = task["metadata"][0]
            targets = []
            for gene in sorted(annotation.genes, key=lambda feature: (
                    feature.contig, feature.start, feature.end, feature.identifier)):
                reason = match_target(gene, annotation.children, task["target_gene"])
                if not reason:
                    continue
                eligible = coding_gene(gene, annotation.children)
                nested = descendants(gene, annotation.children)
                products = sorted({feature.attrs["product"] for feature in nested
                                   if feature.attrs.get("product")})
                proteins = sorted({feature.attrs["protein_id"] for feature in nested
                                   if feature.attrs.get("protein_id")})
                candidate = dict(assembly_accession=assembly, contig=gene.contig,
                                 gene_id=gene.identifier, locus_tag=gene.attrs.get("locus_tag", ""),
                                 gene=gene.attrs.get("gene", ""), products=";".join(products),
                                 protein_ids=";".join(proteins), start_1based=gene.start,
                                 end_1based=gene.end, strand=gene.strand,
                                 match_reason=reason, selected=eligible,
                                 exclusion_reason="" if eligible else "RNA, pseudogene or missing CDS",
                                 seq_id="")
                candidates.append(candidate)
                if eligible:
                    targets.append((gene, candidate))
            for rank, (target_gene, candidate) in enumerate(targets, 1):
                seq_id = f"{assembly.replace('.', '_')}_r{rank:02d}"
                require(len(seq_id) <= 25, f"Region sequence ID exceeds 25 characters: {seq_id}")
                candidate["seq_id"] = seq_id
                selected, before, after = select_neighbors(
                    annotation.genes, target_gene, annotation.children, task["neighbors"])
                selected_ids = {gene.identifier for gene in selected}
                start0 = min(gene.start for gene in selected) - 1
                end0 = max(gene.end for gene in selected)
                record, region, nearby = extract_region(
                    task["binary"], output, work, seq_id, genome, annotation, dna,
                    target_gene, target_gene, start0, end0, seq_id, selected_ids)
                region.update(target_gene_id=target_gene.identifier,
                              requested_neighbors_per_side=task["neighbors"],
                              upstream_genes_selected=min(before, task["neighbors"]),
                              downstream_genes_selected=min(after, task["neighbors"]),
                              target_query=task["target_gene"])
                genus = genome["organism_name"].split()[0]
                genes = []
                for row in nearby:
                    if row["selected_for_plot"]:
                        note = row["gene"] or row["products"] or "Other"
                        genes.append(dict(row, start=row["local_start_1based"],
                                          end=row["local_end_1based"], note=note, organisms=genus))
                results.append(dict(record=record, region=region, neighbors=nearby,
                                    sequence=dict(seq_id=seq_id, length=end0 - start0),
                                    genome=dict(genome, seq_id=seq_id, region_length=end0 - start0),
                                    genes=genes))
    return dict(assembly=assembly, results=results, candidates=candidates,
                source_fields=source_fields, elapsed_seconds=round(time.monotonic() - started, 6),
                pid=os.getpid())


def target_diagnostics(path, candidates):
    if not candidates:
        candidates = [dict(assembly_accession="", contig="", gene_id="", locus_tag="", gene="",
                           products="", protein_ids="", start_1based="", end_1based="", strand="",
                           match_reason="No matching original GFF annotation", selected=False,
                           exclusion_reason="No target matches", seq_id="")]
    write_table(path, candidates)


def publish_preparation(staged, output):
    """Publish only complete preparation results, preserving downstream products."""
    output.mkdir(parents=True, exist_ok=True)
    destination = output / "regions"
    previous = staged.parent / "previous-regions"
    if destination.exists():
        os.replace(destination, previous)
    try:
        os.replace(staged / "regions", destination)
    except OSError:
        if previous.exists():
            os.replace(previous, destination)
        raise
    for path in sorted(staged.iterdir()):
        if path.is_file():
            os.replace(path, output / path.name)
    if previous.exists():
        shutil.rmtree(previous)


def prepare(args, binary, version, output, inputs):
    """Stage independent assembly tasks and merge in the historical input order."""
    started = time.monotonic()
    threads = validate_threads(args.threads)
    metadata, manual, annotation_paths, fasta_paths, hashes = inputs
    protect_region_inputs(output, [ROOT / row["file"] for row in hashes])
    assemblies = list(dict.fromkeys(row["assembly_accession"] for row in metadata))
    workers = min(threads, len(assemblies))
    require(workers > 0, "No assemblies selected for preparation")
    grouped_metadata, manual_groups = defaultdict(list), defaultdict(list)
    for row in metadata:
        grouped_metadata[row["assembly_accession"]].append(row)
    for row in manual or []:
        manual_groups[row["seq_id"]].append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".prepare-", dir=output.parent) as temporary:
        stage = Path(temporary)
        staged = stage / "prepared"
        (staged / "regions").mkdir(parents=True)
        tasks = []
        for index, assembly in enumerate(assemblies):
            rows = grouped_metadata[assembly]
            tasks.append(dict(assembly=assembly, entry=args.entry, metadata=rows,
                              manual_groups={row["seq_id"]: manual_groups[row["seq_id"]]
                                             for row in rows} if manual is not None else {},
                              annotation=str(annotation_paths[assembly]), fasta=str(fasta_paths[assembly]),
                              output=str(stage / "workers" / f"{index:05d}"), binary=str(binary),
                              target_gene=args.target_gene, neighbors=args.neighbors))
        with ProcessPoolExecutor(max_workers=workers, initializer=worker_environment) as pool:
            # map preserves assembly order even when workers finish out of order.
            completed = list(pool.map(assembly_worker, tasks))
        candidates = [row for item in completed for row in item["candidates"]]
        collected = [row for item in completed for row in item["results"]]
        if args.entry == "curated":
            index = {row["region"]["seq_id"]: row for row in collected}
            collected = [index[row["seq_id"]] for row in metadata]
        records = [row["record"] for row in collected]
        if args.entry == "auto":
            target_diagnostics(staged / "target_matches.tsv", candidates)
            if not records:
                output.mkdir(parents=True, exist_ok=True)
                os.replace(staged / "target_matches.tsv", output / "target_matches.tsv")
                require(False, f"No protein-coding target matches for {args.target_gene}; "
                        f"see {output / 'target_matches.tsv'}. Any previous scientific outputs in this "
                        "directory are not results for the requested target.")
        regions = [row["region"] for row in collected]
        neighbors = [gene for row in collected for gene in row["neighbors"]]
        sequences = [row["sequence"] for row in collected]
        genomes = [row["genome"] for row in collected]
        genes = [gene for row in collected for gene in row["genes"]]
        if args.entry == "curated":
            source_fields = {key: value for item in completed for key, value in item["source_fields"].items()}
            nearby_index = {(row["seq_id"], row["gene_id"]): row for row in neighbors}
            for row in manual:
                original = nearby_index[(row["seq_id"], source_fields[
                    (row["seq_id"], row["type"], row["feat_id"])])]
                enriched = dict(row)
                for key, value in original.items():
                    if key != "seq_id":
                        enriched[key if key not in enriched else "source_" + key] = value
                enriched["is_target"] = row["note"] == "Uptake-Gene"
                genes.append(enriched)
        require(len({row["seq_id"] for row in regions}) == len(regions), "Duplicate extracted region IDs")
        require(sum(row["is_target"] for row in genes) == len(regions), "Expected one target marker per region")
        for task in tasks:
            for path in sorted((Path(task["output"]) / "regions").iterdir()):
                destination = staged / "regions" / path.name
                require(not destination.exists(), f"Duplicate staged BED file: {path.name}")
                os.replace(path, destination)
        fasta = staged / "extracted_sequences.fasta"
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
            write_table(staged / name, rows)
        write_sequences(staged / "sequences.csv", sequences)
        context = dict(execution_context(threads), entry=args.entry, assembly_workers=workers,
                       worker_native_threads=1, elapsed_seconds=round(time.monotonic() - started, 6),
                       assemblies=[dict(assembly_accession=item["assembly"], pid=item["pid"],
                                        elapsed_seconds=item["elapsed_seconds"]) for item in completed])
        (staged / "preparation_execution.json").write_text(
            json.dumps(context, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        publish_preparation(staged, output)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entry", choices=("curated", "auto"), required=True)
    parser.add_argument("--genome-dir", type=Path, default=ROOT / "data/genomes")
    parser.add_argument("--genome-metadata", type=Path, default=ROOT / "data/genomes.tsv")
    parser.add_argument("--annotations", type=Path, default=ROOT / "data/genes.tsv")
    parser.add_argument("--target-gene", default="Nramp")
    parser.add_argument("--neighbors", type=int, default=10)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS,
                        help="Shared CPU budget for parallel genome preparation (default: 56)")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bedtools", help="BEDTools 2.31.0 executable")
    parser.add_argument("--verify-only", action="store_true", help="Regenerate independently and compare existing outputs")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    args.threads = validate_threads(args.threads)
    require(args.neighbors >= 0, "Neighbor count must be nonnegative")
    args.target_gene = args.target_gene.strip()
    require(bool(args.target_gene), "Target gene must not be empty")
    resolve = lambda path: (path if path.is_absolute() else ROOT / path).resolve()
    for key in ("genome_dir", "genome_metadata") + (("annotations",) if args.entry == "curated" else ()):
        setattr(args, key, repository_input(getattr(args, key)))
    output = resolve(args.output) if args.output else ROOT / "results" / args.entry
    require(output != ROOT and ROOT / "data" not in (output, *output.parents)
            and args.genome_dir not in (output, *output.parents)
            and not args.genome_dir.is_relative_to(output), "Output must not overwrite source inputs")
    produced = {"sequences.csv", "genomes.tsv", "genes.tsv", "regions.tsv", "neighbors.tsv",
                "extracted_sequences.fasta", "input_hashes.tsv", "extraction_metrics.tsv",
                "extraction_versions.tsv", "target_matches.tsv", "preparation_execution.json"}
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
                if not expected.is_file() or expected.name in {"extraction_versions.tsv", "preparation_execution.json"}:
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
