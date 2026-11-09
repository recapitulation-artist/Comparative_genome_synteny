#!/usr/bin/env python3
"""Extract fixed genomic intervals and plot annotations directly from source GFFs."""

import argparse
from collections import Counter
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile

from extract_regions import (ROOT, annotation_indices, descendants, fasta_bytes,
                             fasta_records, genome_file, read_gff, require,
                             resolve_bedtools, run_bedtools, table, write_table)


MANIFEST_COLUMNS = ("seq_id", "assembly_accession", "contig", "region_start0", "region_end0",
                    "target_protein_id", "target_locus_tag")
CATEGORIES = ("Target Nramp", "Other protein-coding", "RNA gene", "Pseudogene")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_inputs():
    """Check only inputs used by this independent GFF route."""
    sources = [ROOT / "data/sequences.csv", ROOT / "data/genomes.tsv", ROOT / "data/region_manifest.tsv"]
    sequences, genomes, manifest = table(sources[0], ","), table(sources[1]), table(sources[2])
    order = [row["seq_id"] for row in sequences]
    require(len(order) == len(set(order)) == 29, "Expected 29 unique sequence IDs")
    require(len(genomes) == 29 and {row["seq_id"] for row in genomes} == set(order),
            "Genome IDs differ from sequences.csv")
    require(len(manifest) == 29 and [row["seq_id"] for row in manifest] == order,
            "Region manifest must follow sequences.csv order")
    require(tuple(manifest[0]) == MANIFEST_COLUMNS, "Unexpected region manifest columns")
    genome_index = {row["seq_id"]: row for row in genomes}
    require(len({row["assembly_accession"] for row in genomes}) == 29, "Expected 29 unique assemblies")
    for sequence, region in zip(sequences, manifest):
        identifier = sequence["seq_id"]
        require(region["assembly_accession"] == genome_index[identifier]["assembly_accession"],
                f"Manifest/genome assembly differs: {identifier}")
        start0, end0 = int(region["region_start0"]), int(region["region_end0"])
        require(0 <= start0 < end0 and end0 - start0 == int(sequence["length"])
                == int(genome_index[identifier]["region_length"]), f"Invalid region bounds: {identifier}")
        require(all(region[key] for key in MANIFEST_COLUMNS), f"Incomplete region manifest: {identifier}")
        directory = ROOT / "data/genomes" / region["assembly_accession"]
        sources.extend([genome_file(directory, "annotation.gff"), genome_file(directory, "genome.fna")])
    checksums = {}
    for line in (ROOT / "data/checksums.sha256").read_text(encoding="ascii").splitlines():
        digest, name = line.split(maxsplit=1)
        require(name not in checksums, f"Duplicate checksum entry: {name}")
        checksums[name] = digest
    hashes = []
    for source in sources:
        name = source.relative_to(ROOT).as_posix()
        require(name in checksums, f"Missing input checksum: {name}")
        actual = sha256(source)
        require(actual == checksums[name], f"Input checksum differs: {name}")
        hashes.append(dict(file=name, sha256=actual))
    return sequences, genome_index, manifest, hashes


def target_features(features, region):
    lookup, children = annotation_indices(features)
    protein, locus = region["target_protein_id"], region["target_locus_tag"]
    cds_matches, gene_matches = lookup[("CDS", protein)], lookup[("gene", locus)]
    require(len(cds_matches) == len(gene_matches) == 1,
            f"Expected one target CDS and locus: {region['assembly_accession']}/{protein}/{locus}")
    cds, gene = cds_matches[0], gene_matches[0]
    require(cds.contig == gene.contig == region["contig"] and cds.strand == gene.strand
            and gene.start <= cds.start <= cds.end <= gene.end
            and any(child is cds for child in descendants(gene, children)),
            f"Target protein/locus mapping differs: {region['seq_id']}")
    require(int(region["region_start0"]) <= gene.start - 1 < gene.end <= int(region["region_end0"]),
            f"Target gene is not fully inside its region: {region['seq_id']}")
    return cds, gene, children


def extract_one(binary, output, work, sequence, manifest, genome):
    identifier, assembly = sequence["seq_id"], manifest["assembly_accession"]
    directory = ROOT / "data/genomes" / assembly
    features = read_gff(genome_file(directory, "annotation.gff"))
    target, target_gene, children = target_features(features, manifest)
    start0, end0 = int(manifest["region_start0"]), int(manifest["region_end0"])
    contigs = [(name, dna) for name, dna in fasta_records(genome_file(directory, "genome.fna"))
               if name == manifest["contig"]]
    require(len(contigs) == 1, f"Expected one source contig: {assembly}/{manifest['contig']}")
    contig, dna = contigs[0]
    require(0 <= start0 < end0 <= len(dna), f"Region outside source contig: {identifier}")
    header = f"{contig}_{manifest['target_protein_id']}"
    bed = output / "regions" / f"{assembly}.bed"
    bed.write_text(f"{contig}\t{start0}\t{end0}\t{header}\t0\t+\n", encoding="utf-8")
    fasta, extracted = work / "contig.fna", work / "region.fasta"
    fasta.write_bytes(fasta_bytes(contigs))
    fasta.with_suffix(".fna.fai").unlink(missing_ok=True)
    run_bedtools(binary, "getfasta", "-fi", fasta, "-bed", bed, "-fo", extracted, "-nameOnly")
    records = list(fasta_records(extracted))
    require(len(records) == 1 and records[0] == (header, dna[start0:end0]),
            f"BEDTools sequence differs from the source interval: {identifier}")

    genes = [feature for feature in features if feature.kind in {"gene", "pseudogene"}]
    gene_index = {gene.identifier: gene for gene in genes}
    require(len(gene_index) == len(genes) and "" not in gene_index, f"Nonunique GFF gene IDs: {assembly}")
    genes_bed = output / "regions" / f"{assembly}.genes.bed"
    genes_bed.write_text("".join(f"{gene.contig}\t{gene.start - 1}\t{gene.end}\t{gene.identifier}\t0\t{gene.strand}\n"
                                 for gene in genes), encoding="utf-8")
    overlaps = run_bedtools(binary, "intersect", "-a", genes_bed, "-b", bed, "-wa", "-wb")
    neighbors = []
    for line in overlaps.splitlines():
        fields = line.split("\t")
        require(len(fields) == 12, f"Unexpected BEDTools intersection: {identifier}")
        gene = gene_index[fields[3]]
        nested = descendants(gene, children)
        proteins = sorted({child.attrs["protein_id"] for child in nested if child.attrs.get("protein_id")})
        products = sorted({child.attrs["product"] for child in nested if child.attrs.get("product")})
        local_start0, local_end0 = max(gene.start - 1, start0) - start0, min(gene.end, end0) - start0
        require(gene.contig == contig and 0 <= local_start0 < local_end0 <= end0 - start0,
                f"Invalid overlapping gene interval: {identifier}/{gene.identifier}")
        is_target = gene is target_gene
        require((manifest["target_protein_id"] in proteins) == is_target,
                f"Target protein maps to an unexpected gene: {identifier}/{gene.identifier}")
        neighbors.append(dict(seq_id=identifier, assembly_accession=assembly,
                              contig=contig, gene_id=gene.identifier,
                              locus_tag=gene.attrs.get("locus_tag", ""), feature_type=gene.kind,
                              gene=gene.attrs.get("gene", ""), biotype=gene.attrs.get("gene_biotype", ""),
                              start_1based=gene.start, end_1based=gene.end, strand=gene.strand,
                              local_start0=local_start0, local_end0=local_end0,
                              local_start_1based=local_start0 + 1, local_end_1based=local_end0,
                              clipped_left=gene.start - 1 < start0, clipped_right=gene.end > end0,
                              is_target=is_target, protein_ids=";".join(proteins), products=";".join(products)))
    neighbors.sort(key=lambda row: (row["start_1based"], row["end_1based"], row["gene_id"]))
    require(sum(row["is_target"] for row in neighbors) == 1, f"Expected one target gene: {identifier}")
    region = dict(seq_id=identifier, fasta_header=header, assembly_accession=assembly,
                  contig=contig, contig_length=len(dna), region_start0=start0, region_end0=end0,
                  region_start_1based=start0 + 1, region_end_1based=end0, region_strand="+", length=end0 - start0,
                  target_protein_id=manifest["target_protein_id"], target_start_1based=target.start,
                  target_end_1based=target.end, target_strand=target.strand,
                  target_locus_tag=manifest["target_locus_tag"], target_product=target.attrs.get("product", ""),
                  target_local_start_1based=target.start - start0, target_local_end_1based=target.end - start0,
                  assembly_status=genome["assembly_status"], genome_notes=genome.get("genome_notes", ""))
    return records[0], region, neighbors


def category(row):
    if row["is_target"]:
        return "Target Nramp"
    if row["feature_type"] == "pseudogene" or row["biotype"] == "pseudogene":
        return "Pseudogene"
    if row["biotype"] == "protein_coding":
        return "Other protein-coding"
    require(row["biotype"] in {"tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA"},
            f"Unrecognized original gene biotype: {row['seq_id']}/{row['gene_id']}/{row['biotype']}")
    return "RNA gene"


def extract(binary, version, output, inputs):
    sequences, genome_index, manifest, hashes = inputs
    (output / "regions").mkdir(parents=True, exist_ok=True)
    records, regions, neighbors = [], [], []
    with tempfile.TemporaryDirectory(prefix=".raw-extraction-", dir=output) as temporary:
        for sequence, interval in zip(sequences, manifest):
            record, region, nearby = extract_one(binary, output, Path(temporary), sequence,
                                                 interval, genome_index[sequence["seq_id"]])
            records.append(record)
            regions.append(region)
            neighbors.extend(nearby)
    # Retain the historical FASTA serialization without a terminal newline.
    fasta = output / "extracted_29_sequences.fasta"
    fasta.write_bytes(fasta_bytes(records)[:-1])
    genes = [dict(row, start=row["local_start_1based"], end=row["local_end_1based"], note=category(row),
                  organisms=genome_index[row["seq_id"]]["organism_name"].split()[0]) for row in neighbors]
    species = [dict(rank=rank, seq_id=row["seq_id"], organism_name=genome_index[row["seq_id"]]["organism_name"],
                    strain=genome_index[row["seq_id"]]["strain"], assembly_accession=row["assembly_accession"],
                    assembly_status=row["assembly_status"], genome_notes=row["genome_notes"], contig=row["contig"],
                    target_strand=row["target_strand"], target_protein_id=row["target_protein_id"],
                    fasta_header=row["fasta_header"], length=row["length"],
                    reverse_complemented_in_oriented_branch=row["target_strand"] == "-")
               for rank, row in enumerate(regions, 1)]
    counts = Counter(row["note"] for row in genes)
    require(set(counts) == set(CATEGORIES), "Unexpected original annotation categories")
    metrics = dict(regions=len(regions), neighbor_genes=len(neighbors), target_genes=counts["Target Nramp"],
                   other_protein_coding_genes=counts["Other protein-coding"], rna_genes=counts["RNA gene"],
                   pseudogenes=counts["Pseudogene"], clipped_left_genes=sum(row["clipped_left"] for row in neighbors),
                   clipped_right_genes=sum(row["clipped_right"] for row in neighbors),
                   total_bases=sum(len(dna) for _, dna in records), fasta_sha256=sha256(fasta))
    for expected in table(ROOT / "config/expected_raw_extraction_metrics.tsv"):
        require(str(metrics[expected["metric"]]) == expected["value"],
                f"Unexpected raw extraction metric {expected['metric']}: {metrics[expected['metric']]}")
    for name, rows in (("regions.tsv", regions), ("neighbors.tsv", neighbors), ("genes.tsv", genes),
                       ("species.tsv", species), ("extraction_input_hashes.tsv", hashes),
                       ("extraction_metrics.tsv", [dict(metric=key, value=value) for key, value in metrics.items()]),
                       ("extraction_versions.tsv", [dict(tool="bedtools", version=version),
                                                     dict(tool="python", version=sys.version.split()[0])])):
        write_table(output / name, rows)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bedtools", help="BEDTools 2.31.0 executable")
    parser.add_argument("--output", type=Path, default=ROOT / "results/raw")
    parser.add_argument("--check-inputs", action="store_true", help="Verify only the input hashes and manifest")
    parser.add_argument("--verify-only", action="store_true", help="Re-extract independently and compare existing outputs")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    require(not (args.check_inputs and args.verify_only), "Select only one check mode")
    inputs = checked_inputs()
    if args.check_inputs:
        print(f"Verified {len(inputs[3])} raw-route input hashes and 29 manifest intervals.")
        return
    binary, version = resolve_bedtools(args.bedtools)
    output = args.output.resolve()
    require(output != ROOT and ROOT / "data" not in (output, *output.parents),
            "Output must not overwrite repository inputs")
    if args.verify_only:
        require(output.is_dir(), f"Missing raw extraction output: {output}")
        versions = table(output / "extraction_versions.tsv")
        require({row["tool"] for row in versions} == {"bedtools", "python"}
                and len(versions) == 2 and all(row["version"] for row in versions)
                and next(row["version"] for row in versions if row["tool"] == "bedtools") == version,
                "Missing or incompatible raw extraction version metadata")
        with tempfile.TemporaryDirectory(prefix="raw-extraction-check-") as temporary:
            regenerated = Path(temporary)
            metrics = extract(binary, version, regenerated, inputs)
            for expected in regenerated.rglob("*"):
                if not expected.is_file() or expected.name == "extraction_versions.tsv":
                    continue
                actual = output / expected.relative_to(regenerated)
                require(actual.is_file() and actual.read_bytes() == expected.read_bytes(),
                        f"Raw extraction output differs: {actual}")
        print("Verified raw extraction outputs against source FASTA/GFF annotations.")
    else:
        metrics = extract(binary, version, output, inputs)
    print(f"Extracted {metrics['regions']} regions and {metrics['neighbor_genes']} original GFF genes "
          f"with {metrics['target_genes']} targets.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(f"Raw extraction failed: {error}")
