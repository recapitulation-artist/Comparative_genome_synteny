#!/usr/bin/env python3
"""Recreate historical regions from genome annotations and verify the FASTA."""

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


BEDTOOLS_VERSION = "bedtools v2.31.0"
ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def table(path, delimiter="\t"):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(path, rows):
    require(bool(rows), f"No rows for {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: str(value).lower() if isinstance(value, bool) else value
                         for key, value in row.items()} for row in rows)


def fasta_records(path):
    with open_text(path) as handle:
        name, parts = None, []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                require(name is not None, f"Sequence before FASTA header: {path}")
                parts.append(line)
        if name is not None:
            yield name, "".join(parts)


def fasta_bytes(records):
    return "".join(f">{name}\n{sequence}\n" for name, sequence in records).encode("ascii")


@dataclass
class Feature:
    contig: str
    kind: str
    start: int
    end: int
    strand: str
    attrs: dict

    @property
    def identifier(self):
        return self.attrs.get("ID", "")


def read_gff(path):
    features = []
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
    return features


def run_bedtools(binary, *arguments):
    command = [str(binary), *map(str, arguments)]
    result = subprocess.run(command, text=True, capture_output=True)
    require(result.returncode == 0, f"BEDTools failed ({' '.join(command)}): {result.stderr.strip()}")
    return result.stdout


def resolve_bedtools(requested):
    local = ROOT / ".tools" / "bedtools"
    candidate = requested or (str(local) if local.is_file() else shutil.which("bedtools"))
    require(candidate is not None, "BEDTools is missing; run make setup-bedtools or use --bedtools PATH")
    binary = Path(candidate).expanduser()
    if not binary.is_absolute():
        binary = ROOT / binary if (ROOT / binary).is_file() else Path(shutil.which(str(binary)) or binary)
    require(binary.is_file() and os.access(binary, os.X_OK), f"BEDTools executable not found: {candidate}")
    binary = binary.resolve()
    version = run_bedtools(binary, "--version").strip()
    require(version == BEDTOOLS_VERSION, f"Expected {BEDTOOLS_VERSION}; found {version}")
    return binary, version


def genome_file(directory, name):
    path = directory / (name + ".gz")
    if not path.is_file():
        path = directory / name
    require(path.is_file(), f"Missing genome input: {path}")
    return path


def annotation_indices(features):
    lookup, children = defaultdict(list), defaultdict(list)
    for feature in features:
        if feature.kind == "CDS" and feature.attrs.get("protein_id"):
            lookup[("CDS", feature.attrs["protein_id"])].append(feature)
        if feature.kind in {"gene", "pseudogene"} and feature.attrs.get("locus_tag"):
            lookup[("gene", feature.attrs["locus_tag"])].append(feature)
        for parent in feature.attrs.get("Parent", "").split(","):
            if parent:
                children[parent].append(feature)
    return lookup, children


def descendants(gene, children):
    found, pending, seen = [], [gene], set()
    while pending:
        feature = pending.pop()
        if id(feature) in seen:
            continue
        seen.add(id(feature))
        found.append(feature)
        pending.extend(children.get(feature.identifier, []))
    return found


def prepare_region(sequence, legacy, genome_index, features):
    seq_id, length = sequence["seq_id"], int(sequence["length"])
    require(len({(row["type"], row["feat_id"]) for row in legacy}) == len(legacy),
            f"Duplicate legacy annotations: {seq_id}")
    require(len({row["geom_id"] for row in legacy}) == 1, f"Multiple assemblies for {seq_id}")
    assembly = legacy[0]["geom_id"]
    require(genome_index["assembly_accession"] == assembly, f"Assembly index mismatch: {seq_id}")
    require(int(genome_index["region_length"]) == length, f"Region length index mismatch: {seq_id}")
    targets = [row for row in legacy if row["note"] == "Uptake-Gene"]
    require(len(targets) == 1, f"Expected one uptake target: {seq_id}")
    old_target = targets[0]
    lookup, children = annotation_indices(features)
    mapped = []
    for row in legacy:
        require(row["type"] in {"CDS", "gene"}, f"Unsupported legacy feature type: {row['type']}")
        matches = lookup[(row["type"], row["feat_id"])]
        require(len(matches) == 1, f"Expected one GFF match: {assembly}/{row['feat_id']}, found {len(matches)}")
        require(1 <= int(row["start"]) <= int(row["end"]) <= length, f"Invalid legacy interval: {row['feat_id']}")
        mapped.append((row, matches[0]))
    target = next(feature for row, feature in mapped if row is old_target)
    # Preserve the reference's historical one-base boundary shift.
    start0 = target.start - int(old_target["start"]) + 1
    end0 = start0 + length
    for row, feature in mapped:
        require(feature.contig == target.contig and feature.strand == row["strand"],
                f"Legacy contig/strand mismatch: {assembly}/{row['feat_id']}")
        require(feature.start == start0 + int(row["start"]) - 1
                and feature.end == start0 + int(row["end"]),
                f"Historical coordinate mismatch: {assembly}/{row['feat_id']}")
    require(min(feature.start for _, feature in mapped) == start0
            and max(feature.end for _, feature in mapped) == end0,
            f"Legacy annotations do not define the complete region: {seq_id}")
    require(min(int(row["start"]) for row in legacy) == 1
            and max(int(row["end"]) for row in legacy) == length,
            f"Legacy local boundaries mismatch: {seq_id}")
    return assembly, target, old_target, start0, end0, children


def extract_one(binary, output, work, sequence, legacy, genome_index):
    directory = ROOT / "data" / "genomes" / genome_index["assembly_accession"]
    features = read_gff(genome_file(directory, "annotation.gff"))
    assembly, target, old_target, start0, end0, children = prepare_region(
        sequence, legacy, genome_index, features)
    matches = [(name, dna) for name, dna in fasta_records(genome_file(directory, "genome.fna"))
               if name == target.contig]
    require(len(matches) == 1, f"Expected unique genome contig: {assembly}/{target.contig}")
    name, dna = matches[0]
    require(0 <= start0 < end0 <= len(dna), f"Region outside contig: {assembly}")
    header = f"{target.contig}_{old_target['feat_id']}"
    region_bed = output / "regions" / f"{assembly}.bed"
    region_bed.write_text(f"{name}\t{start0}\t{end0}\t{header}\t0\t+\n", encoding="utf-8")
    # BEDTools may create a .fai beside its FASTA. Index only this temporary contig.
    temporary_fasta, temporary_output = work / "contig.fna", work / "region.fasta"
    temporary_fasta.write_bytes(fasta_bytes([(name, dna)]))
    temporary_fasta.with_suffix(".fna.fai").unlink(missing_ok=True)
    run_bedtools(binary, "getfasta", "-fi", temporary_fasta, "-bed", region_bed,
                 "-fo", temporary_output, "-nameOnly")
    records = list(fasta_records(temporary_output))
    require(len(records) == 1 and records[0][0] == header and len(records[0][1]) == end0 - start0,
            f"Unexpected BEDTools FASTA output: {assembly}")
    genes = [feature for feature in features if feature.kind in {"gene", "pseudogene"}]
    gene_lookup = {gene.identifier: gene for gene in genes}
    require(len(gene_lookup) == len(genes) and "" not in gene_lookup, f"Nonunique gene IDs: {assembly}")
    genes_bed = output / "regions" / f"{assembly}.genes.bed"
    genes_bed.write_text("".join(f"{g.contig}\t{g.start - 1}\t{g.end}\t{g.identifier}\t0\t{g.strand}\n"
                                 for g in genes), encoding="utf-8")
    overlap = run_bedtools(binary, "intersect", "-a", genes_bed, "-b", region_bed, "-wa", "-wb")
    neighbors = []
    for line in overlap.splitlines():
        fields = line.split("\t")
        require(len(fields) == 12, f"Unexpected BEDTools intersection: {assembly}")
        gene = gene_lookup[fields[3]]
        nested = descendants(gene, children)
        products = sorted({f.attrs["product"] for f in nested if f.attrs.get("product")})
        proteins = sorted({f.attrs["protein_id"] for f in nested if f.attrs.get("protein_id")})
        local_start0, local_end0 = max(gene.start - 1, start0) - start0, min(gene.end, end0) - start0
        neighbors.append(dict(seq_id=sequence["seq_id"], assembly_accession=assembly,
                              contig=gene.contig, gene_id=gene.identifier,
                              locus_tag=gene.attrs.get("locus_tag", ""), feature_type=gene.kind,
                              gene=gene.attrs.get("gene", ""), biotype=gene.attrs.get("gene_biotype", ""),
                              start_1based=gene.start, end_1based=gene.end, strand=gene.strand,
                              local_start0=local_start0, local_end0=local_end0,
                              local_start_1based=local_start0 + 1, local_end_1based=local_end0,
                              clipped_left=gene.start - 1 < start0, clipped_right=gene.end > end0,
                              is_target=old_target["feat_id"] in proteins,
                              protein_ids=";".join(proteins), products=";".join(products)))
    neighbors.sort(key=lambda row: (row["start_1based"], row["end_1based"], row["gene_id"]))
    require(sum(row["is_target"] for row in neighbors) == 1, f"Expected one overlapping target gene: {assembly}")
    region = dict(seq_id=sequence["seq_id"], fasta_header=header, assembly_accession=assembly,
                  contig=target.contig, contig_length=len(dna), region_start0=start0, region_end0=end0,
                  region_start_1based=start0 + 1, region_end_1based=end0, region_strand="+", length=end0 - start0,
                  target_protein_id=old_target["feat_id"], target_start_1based=target.start,
                  target_end_1based=target.end, target_strand=target.strand,
                  target_locus_tag=target.attrs.get("locus_tag", ""), target_product=target.attrs.get("product", ""),
                  target_local_start_1based=target.start - start0, target_local_end_1based=target.end - start0,
                  legacy_target_start=int(old_target["start"]), legacy_target_end=int(old_target["end"]),
                  legacy_annotations_verified=len(legacy), assembly_status=genome_index["assembly_status"],
                  genome_notes=genome_index.get("genome_notes", ""))
    return records[0], region, neighbors


def validate_reference(records, output):
    # The reference fixture is deliberately opened only after genomic extraction.
    reference_path = ROOT / "data" / "reference_sequences.fasta"
    expected = list(fasta_records(reference_path))
    require(len(expected) == len(records), "Reference FASTA record count differs")
    validation = []
    for (header, dna), (expected_header, expected_dna) in zip(records, expected):
        same = header == expected_header and dna == expected_dna
        validation.append(dict(fasta_header=header, length=len(dna),
                               sequence_sha256=hashlib.sha256(dna.encode("ascii")).hexdigest(),
                               reference_sha256=hashlib.sha256(expected_dna.encode("ascii")).hexdigest(),
                               header_matches=header == expected_header, sequence_matches=dna == expected_dna))
        require(same, f"Extracted sequence differs from reference: {header}")
    payload = (output / "extracted_29_sequences.fasta").read_bytes()
    require(payload == reference_path.read_bytes(), "FASTA bytes differ from the reference fixture")
    write_table(output / "extraction_validation.tsv", validation)
    return hashlib.sha256(payload).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bedtools", help="BEDTools 2.31.0 executable")
    parser.add_argument("--output", type=Path, default=ROOT / "results", help="Output directory (default: results)")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    binary, version = resolve_bedtools(args.bedtools)
    sequences = table(ROOT / "data" / "sequences.csv", ",")
    legacy = table(ROOT / "data" / "genes.tsv")
    require(len(legacy) == 579, "Expected 579 legacy annotations")
    indices = table(ROOT / "data" / "genomes.tsv")
    order = [row["seq_id"] for row in sequences]
    require(len(order) == 29 and len(set(order)) == 29, "Expected 29 unique sequence IDs")
    require({row["seq_id"] for row in legacy} == set(order), "Legacy annotation sequence IDs differ")
    require(len(indices) == 29 and {row["seq_id"] for row in indices} == set(order), "Genome index IDs differ")
    require(len({row["assembly_accession"] for row in indices}) == 29, "Expected 29 unique assemblies")
    index = {row["seq_id"]: row for row in indices}
    output = args.output.resolve()
    (output / "regions").mkdir(parents=True, exist_ok=True)
    records, regions, neighbors = [], [], []
    with tempfile.TemporaryDirectory(prefix=".extraction-", dir=output) as temporary:
        for sequence in sequences:
            rows = [row for row in legacy if row["seq_id"] == sequence["seq_id"]]
            record, region, nearby = extract_one(binary, output, Path(temporary), sequence, rows, index[sequence["seq_id"]])
            records.append(record)
            regions.append(region)
            neighbors.extend(nearby)
    # The historical file uses one sequence line per record and no final newline.
    (output / "extracted_29_sequences.fasta").write_bytes(fasta_bytes(records)[:-1])
    write_table(output / "regions.tsv", regions)
    write_table(output / "neighbors.tsv", neighbors)
    digest = validate_reference(records, output)
    metrics = dict(regions=len(records), legacy_annotations_verified=len(legacy), neighbor_genes=len(neighbors),
                   target_genes=sum(row["is_target"] for row in neighbors),
                   clipped_left_genes=sum(row["clipped_left"] for row in neighbors),
                   clipped_right_genes=sum(row["clipped_right"] for row in neighbors),
                   total_bases=sum(len(dna) for _, dna in records), reference_records_verified=len(records),
                   reference_bytes_identical=True, fasta_sha256=digest)
    for expected in table(ROOT / "config" / "expected_extraction_metrics.tsv"):
        actual = metrics[expected["metric"]]
        actual = str(actual).lower() if isinstance(actual, bool) else str(actual)
        require(actual == expected["value"],
                f"Unexpected extraction metric {expected['metric']}: {actual}")
    write_table(output / "extraction_metrics.tsv", [dict(metric=key, value=value) for key, value in metrics.items()])
    write_table(output / "extraction_versions.tsv", [dict(tool="bedtools", version=version),
                                                    dict(tool="python", version=sys.version.split()[0])])
    print(f"Verified {len(records)} regions, {len(legacy)} legacy annotations and {len(neighbors)} overlapping genes.")
    print(f"FASTA is byte-identical to the reference: {digest}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"error: {error}")
