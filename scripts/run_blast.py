#!/usr/bin/env python3
"""Compare extracted regions with local BLAST+ and retain native HSP geometry."""

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal, InvalidOperation
from pathlib import Path

from extract_regions import ROOT, fasta_bytes, fasta_records, require, table, write_table
from runtime import DEFAULT_THREADS, execution_context, thread_environment, validate_threads


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
    path = Path(candidate).expanduser()
    if not path.is_absolute():
        path = ROOT / path if (ROOT / path).is_file() else Path(shutil.which(str(path)) or str(path))
    require(path.is_file() and os.access(path, os.X_OK), f"Executable not found: {candidate}")
    path = path.resolve()
    actual = run([path, "-version"], environment).splitlines()[0]
    require(actual == f"{name}: {version}", f"Expected {name} {version}; found {actual}")
    return path, actual


def input_sources(input_dir):
    if (input_dir / "sequences.csv").is_file():
        return dict(sequences=input_dir / "sequences.csv", regions=input_dir / "regions.tsv",
                    extracted_fasta=input_dir / "extracted_sequences.fasta")
    require(input_dir in {(ROOT / "results").resolve(), (ROOT / "results/raw").resolve()},
            f"Missing entry sequence index: {input_dir / 'sequences.csv'}")
    return dict(sequences=ROOT / "data/sequences.csv", regions=input_dir / "regions.tsv",
                extracted_fasta=input_dir / "extracted_29_sequences.fasta")


def prepare_input(output, input_dir=ROOT / "results"):
    sources = input_sources(input_dir)
    sequences = table(sources["sequences"], ",")
    regions = table(sources["regions"])
    records = list(fasta_records(sources["extracted_fasta"]))
    order = [row["seq_id"] for row in sequences]
    require(len(order) >= 2 and len(set(order)) == len(order), "BLAST needs at least two unique sequence IDs")
    require(all(re.fullmatch(r"[A-Za-z0-9_.-]+", name) for name in order),
            "Sequence IDs must contain only letters, digits, underscores, dots and hyphens")
    require(len(records) == len(regions) == len(order), "Region/FASTA/sequence counts differ")
    region_index = {row["seq_id"]: row for row in regions}
    fasta_index = dict(records)
    require(len(region_index) == len(order) and set(region_index) == set(order), "Region IDs differ from sequences.csv")
    require(len(fasta_index) == len(records) and len({row["fasta_header"] for row in regions}) == len(order)
            and set(fasta_index) == {row["fasta_header"] for row in regions}, "FASTA headers differ from regions.tsv")
    normalized, lookup, lengths = [], [], {}
    for rank, sequence in enumerate(sequences):
        region = region_index[sequence["seq_id"]]
        header = region["fasta_header"]
        dna = fasta_index[header].upper()
        require(len(dna) > 0 and len(dna) == int(region["length"]) == int(sequence["length"]),
                f"FASTA length differs for {sequence['seq_id']}")
        require(set(dna) <= set("ACGTN"), f"Unexpected nucleotide: {header}")
        normalized.append((sequence["seq_id"], dna))
        lengths[sequence["seq_id"]] = len(dna)
        lookup.append(dict(rank=rank, seq_id=sequence["seq_id"], fasta_header=header, length=len(dna),
                           sequence_sha256=hashlib.sha256(dna.encode("ascii")).hexdigest()))
    (output / "sequences.fasta").write_bytes(fasta_bytes(normalized))
    write_table(output / "sequence_lookup.tsv", lookup)
    return {name: rank for rank, name in enumerate(order)}, lengths, sources


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


def validate_hsp(fields, lengths):
    values = [decimal(value) for value in fields[2:]]
    identity, length, mismatches, gaps, start, end, start2, end2, evalue, bitscore = values
    require(all(value == value.to_integral_value() for value in values[1:8]),
            "BLAST lengths, mismatch counts and coordinates must be integers")
    require(0 <= identity <= 100 and length > 0 and 0 <= mismatches <= length
            and 0 <= gaps <= length and evalue >= 0 and bitscore >= 0, "Invalid BLAST HSP statistics")
    for identifier, first, last in ((fields[0], start, end), (fields[1], start2, end2)):
        require(1 <= min(first, last) <= max(first, last) <= lengths[identifier],
                f"BLAST coordinates outside extracted sequence: {identifier}/{first}/{last}")
        require(abs(last - first) + 1 <= length, "BLAST alignment length is shorter than its coordinate span")


def convert_hits(raw_path, destination, ranks, lengths):
    rows, pairs = [], set()
    counts = dict(raw_hsps=0, retained_hsps=0, self_hsps_removed=0, reciprocal_hsps_removed=0)
    with raw_path.open(encoding="ascii") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == len(COLUMNS), "BLAST output must have exactly 12 columns")
            query, subject = [sequence_id(value, ranks) for value in fields[:2]]
            fields[:2] = [query, subject]
            validate_hsp(fields, lengths)
            counts["raw_hsps"] += 1
            if ranks[query] == ranks[subject]:
                counts["self_hsps_removed"] += 1
                continue
            if ranks[query] > ranks[subject]:
                counts["reciprocal_hsps_removed"] += 1
                continue
            # Keep native HSP order and strand-directed coordinates.
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
    parser.add_argument("--output", type=Path, help="Output CSV (default: entry directory/blast_hits.csv)")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "results",
                        help="Extracted FASTA and region table directory (default: results)")
    parser.add_argument("--work-dir", type=Path, help="BLAST work and audit directory (default: entry directory/blast)")
    parser.add_argument("--validate-reference", action=argparse.BooleanOptionalAction, default=None,
                        help="Audit against historical HSPs (enabled only for the legacy results directory by default)")
    parser.add_argument("--check-tools", action="store_true", help="Check tool versions and exit")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="CPU budget (default: 56)")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    config = json.loads((ROOT / "config" / "blast.json").read_text(encoding="utf-8"))
    config["num_threads"] = args.threads
    require(config["outfmt"] == "6 " + FIELDS, "Expected the standard 12-field BLAST output")
    environment = {**os.environ, **thread_environment(), **BLAST_ENV}
    blastn, blastn_version = resolve_tool("blastn", args.blastn, config["version"], environment)
    makeblastdb, db_version = resolve_tool("makeblastdb", args.makeblastdb, config["version"], environment)
    if args.check_tools:
        print(f"{blastn_version}; {db_version}")
        return
    config["num_threads"] = validate_threads(args.threads)
    execution = execution_context(args.threads)
    resolve = lambda path: (path if path.is_absolute() else ROOT / path).resolve()
    input_dir = resolve(args.input_dir)
    require(ROOT in input_dir.parents, "BLAST inputs must be inside the repository")
    legacy = input_dir == (ROOT / "results").resolve() and not (input_dir / "sequences.csv").is_file()
    output = resolve(args.work_dir) if args.work_dir else input_dir / "blast"
    destination = resolve(args.output) if args.output else (ROOT / "data/blast_hits.csv" if legacy
                                                           else input_dir / "blast_hits.csv")
    validate_historical = legacy if args.validate_reference is None else args.validate_reference
    require(output != ROOT and ROOT / "data" not in (output, *output.parents),
            "BLAST work directory must not overwrite repository inputs")
    output.mkdir(parents=True, exist_ok=True)
    ranks, lengths, sources = prepare_input(output, input_dir)
    require(all(ROOT in path.resolve().parents for path in sources.values()),
            "BLAST inputs must not link to files outside the repository")
    require(destination not in {path.resolve() for path in sources.values()}, "BLAST output would overwrite an input")
    sources.update(config=ROOT / "config/blast.json", script=Path(__file__).resolve(),
                   runtime_helpers=ROOT / "scripts/runtime.py")
    input_hashes = {key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in sources.items()}
    options = []
    for name in BLAST_OPTIONS:
        value = config[name]
        options.extend(["-" + name, str(value).lower() if isinstance(value, bool) else str(value)])
    portable = lambda path: str(path.relative_to(ROOT)) if ROOT in path.parents else str(path)
    templates = dict(environment={**thread_environment(), **BLAST_ENV},
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
    rows, counts = convert_hits(output / "all_hits.tsv", destination, ranks, lengths)
    for key, path in sources.items():
        require(hashlib.sha256(path.read_bytes()).hexdigest() == input_hashes[key], f"BLAST input changed during execution: {path}")
    provenance = dict(schema_version=1, threads=config["num_threads"], execution=execution,
                      inputs={key: dict(path=portable(path), sha256=input_hashes[key]) for key, path in sources.items()},
                      tools={name: dict(path=str(path), version=version,
                                        sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                             for name, path, version in (("blastn", blastn, blastn_version),
                                                         ("makeblastdb", makeblastdb, db_version))},
                      historical_validation=validate_historical,
                      output=dict(path=portable(destination), sha256=counts["generated_csv_sha256"]))
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if validate_historical:
        validate_reference(rows, counts, output)
    else:
        counts.update(input_sequences=len(ranks), coordinate_ranges_valid=True, historical_validation=False)
        write_table(output / "metrics.tsv", [dict(metric=key, value=value) for key, value in counts.items()])
    print(f"Generated {len(rows)} HSPs across {counts['sequence_pairs']} sequence pairs: {destination}")
    if validate_historical:
        print("Historical alignment geometry and row order match in all 11 non-E-value columns.")
        print(f"Actual BLAST E-values retained; {counts['reference_evalue_differences']} differ from the historical file.")
    else:
        print(f"Validated coordinates for {len(ranks)} input sequences; actual BLAST E-values retained.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"error: {error}")
