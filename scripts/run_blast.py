#!/usr/bin/env python3
"""Align the extracted regions with local BLAST+ and validate historical geometry."""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal, InvalidOperation
from pathlib import Path

from extract_regions import ROOT, fasta_bytes, fasta_records, require, table, write_table


COLUMNS = ["seq_id", "seq_id2", "pident", "length", "mismatch", "gapopen", "start", "end",
           "start2", "end2", "evalue", "bitscore"]
FIELDS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
BLAST_OPTIONS = ("task", "word_size", "reward", "penalty", "gapopen", "gapextend", "evalue",
                 "dust", "soft_masking", "strand", "max_target_seqs", "num_threads")
BLAST_ENV = {"NCBI_DONT_USE_NCBIRC": "1", "BLAST_USAGE_REPORT": "false"}


def run(command, environment):
    result = subprocess.run(list(map(str, command)), text=True, capture_output=True, env=environment)
    require(result.returncode == 0,
            f"{Path(command[0]).name} failed: {result.stderr.strip() or result.stdout.strip()}")
    return result.stdout


def resolve_tool(name, requested, version, environment):
    local = ROOT / ".tools" / "blast" / "bin" / name
    candidate = requested or (str(local) if local.is_file() else shutil.which(name))
    require(candidate, f"{name} is missing; run make setup-blast or provide --{name} PATH")
    path = Path(candidate)
    if not path.is_absolute():
        path = ROOT / path if (ROOT / path).is_file() else Path(shutil.which(str(path)) or str(path))
    require(path.is_file(), f"Executable not found: {candidate}")
    path = path.resolve()
    actual = run([path, "-version"], environment).splitlines()[0]
    require(actual == f"{name}: {version}", f"Expected {name} {version}; found {actual}")
    return path, actual


def prepare_input(output, input_dir=ROOT / "results"):
    sequences = table(ROOT / "data" / "sequences.csv", ",")
    regions = table(input_dir / "regions.tsv")
    records = list(fasta_records(input_dir / "extracted_29_sequences.fasta"))
    order = [row["seq_id"] for row in sequences]
    require(len(order) == 29 and len(set(order)) == 29, "Expected 29 unique sequence IDs")
    require(len(records) == len(regions) == len(order), "Region/FASTA/sequence counts differ")
    require([row["seq_id"] for row in regions] == order, "Region order differs from sequences.csv")
    normalized, lookup = [], []
    for rank, ((header, dna), region, sequence) in enumerate(zip(records, regions, sequences)):
        require(header == region["fasta_header"], f"FASTA header differs for {sequence['seq_id']}")
        require(len(dna) == int(region["length"]) == int(sequence["length"]),
                f"FASTA length differs for {sequence['seq_id']}")
        require(set(dna.upper()) <= set("ACGTN"), f"Unexpected nucleotide: {header}")
        normalized.append((sequence["seq_id"], dna))
        lookup.append(dict(rank=rank, seq_id=sequence["seq_id"], fasta_header=header, length=len(dna),
                           sequence_sha256=hashlib.sha256(dna.encode("ascii")).hexdigest()))
    (output / "sequences.fasta").write_bytes(fasta_bytes(normalized))
    write_table(output / "sequence_lookup.tsv", lookup)
    return {name: rank for rank, name in enumerate(order)}


def decimal(value):
    try:
        result = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"Invalid BLAST numeric value: {value}") from error
    require(result.is_finite(), f"Nonfinite BLAST value: {value}")
    return result


def compact(value):
    formatted = format(decimal(value), "f")
    return formatted.rstrip("0").rstrip(".") if "." in formatted else formatted


def sequence_id(value, ranks):
    if value in ranks:
        return value
    # makeblastdb recognizes accession-like names and adds an NCBI ID wrapper.
    candidates = [name for name in ranks if value in {f"ref|{name}|", f"lcl|{name}"}]
    require(len(candidates) == 1, f"BLAST returned an unknown sequence ID: {value}")
    return candidates[0]


def convert_hits(raw_path, destination, ranks):
    rows, pairs = [], set()
    counts = dict(raw_hsps=0, retained_hsps=0, self_hsps_removed=0, reciprocal_hsps_removed=0)
    with raw_path.open(encoding="ascii") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == len(COLUMNS), "BLAST output must have exactly 12 columns")
            query, subject = [sequence_id(value, ranks) for value in fields[:2]]
            fields[:2] = [query, subject]
            counts["raw_hsps"] += 1
            if ranks[query] == ranks[subject]:
                counts["self_hsps_removed"] += 1
                continue
            if ranks[query] > ranks[subject]:
                counts["reciprocal_hsps_removed"] += 1
                continue
            # Keep native HSP order and strand-directed coordinates.
            for value in fields[2:]:
                decimal(value)
            fields[2] = compact(fields[2])
            fields[10] = "0" if decimal(fields[10]) == 0 else fields[10].lower()
            fields[11] = compact(fields[11])
            rows.append(dict(zip(COLUMNS, fields)))
            pairs.add((query, subject))
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts.update(retained_hsps=len(rows), sequence_pairs=len(pairs))
    counts["generated_csv_sha256"] = hashlib.sha256(destination.read_bytes()).hexdigest()
    return rows, counts


def validate_reference(rows, counts, output):
    # Reference values are used only to audit the completed BLAST result.
    reference = table(ROOT / "data" / "reference_blast_hits.csv", ",")
    require(reference and list(reference[0]) == COLUMNS, "Historical BLAST columns differ")
    same_count = len(rows) == len(reference)
    validation = []
    for column in COLUMNS:
        equal = sum(actual[column] == expected[column] if column in COLUMNS[:2]
                    else decimal(actual[column]) == decimal(expected[column])
                    for actual, expected in zip(rows, reference))
        validation.append(dict(column=column, generated_rows=len(rows), reference_rows=len(reference),
                               equal_values=equal, differing_values=max(len(rows), len(reference)) - equal,
                               all_equal=same_count and equal == len(rows)))
    write_table(output / "validation.tsv", validation)
    differences = []
    for number, (actual, expected) in enumerate(zip(rows, reference), 1):
        if decimal(actual["evalue"]) != decimal(expected["evalue"]):
            differences.append(dict(row=number, **{key: actual[key] for key in
                                    ("seq_id", "seq_id2", "start", "end", "start2", "end2")},
                                    generated_evalue=actual["evalue"], reference_evalue=expected["evalue"]))
    if differences:
        write_table(output / "evalue_differences.tsv", differences)
    else:
        (output / "evalue_differences.tsv").unlink(missing_ok=True)
    counts["reference_rows"] = len(reference)
    counts["reference_evalue_differences"] = next(row["differing_values"] for row in validation
                                                  if row["column"] == "evalue")
    counts["reference_geometry_matches"] = all(row["all_equal"] for row in validation
                                               if row["column"] != "evalue")
    write_table(output / "metrics.tsv", [dict(metric=key, value=value) for key, value in counts.items()])
    require(counts["reference_geometry_matches"],
            f"BLAST geometry/order differs from the historical reference; see {output / 'validation.tsv'}")
    for expected in table(ROOT / "config" / "expected_blast_metrics.tsv"):
        actual = counts[expected["metric"]]
        actual = str(actual).lower() if isinstance(actual, bool) else str(actual)
        require(actual == expected["value"], f"Unexpected BLAST metric {expected['metric']}: {actual}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blastn", help="Local blastn executable")
    parser.add_argument("--makeblastdb", help="Local makeblastdb executable")
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "blast_hits.csv",
                        help="Output CSV (default: data/blast_hits.csv)")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "results",
                        help="Extracted FASTA and region table directory (default: results)")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "results/blast",
                        help="BLAST work and audit directory (default: results/blast)")
    parser.add_argument("--check-tools", action="store_true", help="Check tool versions and exit")
    parser.add_argument("--threads", type=int, help="Resource override (default: locked configuration)")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    config = json.loads((ROOT / "config" / "blast.json").read_text(encoding="utf-8"))
    if args.threads is not None:
        require(args.threads > 0, "Thread count must be positive")
        config["num_threads"] = args.threads
    require(config["outfmt"] == "6 " + FIELDS, "Expected the standard 12-field BLAST output")
    environment = {**os.environ, **BLAST_ENV}
    blastn, blastn_version = resolve_tool("blastn", args.blastn, config["version"], environment)
    makeblastdb, db_version = resolve_tool("makeblastdb", args.makeblastdb, config["version"], environment)
    if args.check_tools:
        print(f"{blastn_version}; {db_version}")
        return
    resolve = lambda path: (path if path.is_absolute() else ROOT / path).resolve()
    output, input_dir = resolve(args.work_dir), resolve(args.input_dir)
    require(output != ROOT and ROOT / "data" not in (output, *output.parents),
            "BLAST work directory must not overwrite repository inputs")
    output.mkdir(parents=True, exist_ok=True)
    ranks = prepare_input(output, input_dir)
    options = []
    for name in BLAST_OPTIONS:
        value = config[name]
        options.extend(["-" + name, str(value).lower() if isinstance(value, bool) else str(value)])
    portable = lambda path: str(path.relative_to(ROOT)) if ROOT in path.parents else str(path)
    templates = dict(environment=BLAST_ENV,
                     makeblastdb=["makeblastdb", "-in", portable(output / "sequences.fasta"), "-dbtype", "nucl",
                                  "-parse_seqids", "-blastdb_version", str(config["db_version"]),
                                  "-out", "<temporary>/regions"],
                     blastn=["blastn", "-query", portable(output / "sequences.fasta"), "-db", "<temporary>/regions",
                             "-out", portable(output / "all_hits.tsv"), "-outfmt", config["outfmt"], *options])
    (output / "commands.json").write_text(json.dumps(templates, indent=2) + "\n", encoding="utf-8")
    write_table(output / "versions.tsv", [dict(tool="blastn", version=blastn_version),
                                            dict(tool="makeblastdb", version=db_version),
                                            dict(tool="python", version=sys.version.split()[0])])
    with tempfile.TemporaryDirectory(prefix=".blast-", dir=output) as temporary:
        database = Path(temporary) / "regions"
        run([makeblastdb, "-in", output / "sequences.fasta", "-dbtype", "nucl", "-parse_seqids",
             "-blastdb_version", str(config["db_version"]), "-out", database], environment)
        run([blastn, "-query", output / "sequences.fasta", "-db", database, "-out", output / "all_hits.tsv",
             "-outfmt", config["outfmt"], *options], environment)
    destination = args.output if args.output.is_absolute() else ROOT / args.output
    rows, counts = convert_hits(output / "all_hits.tsv", destination, ranks)
    validate_reference(rows, counts, output)
    print(f"Generated {len(rows)} HSPs across {counts['sequence_pairs']} sequence pairs: {destination}")
    print("Historical alignment geometry and row order match in all 11 non-E-value columns.")
    print(f"Actual BLAST E-values retained; {counts['reference_evalue_differences']} differ from the historical file.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"error: {error}")
