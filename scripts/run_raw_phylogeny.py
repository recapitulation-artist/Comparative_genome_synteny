#!/usr/bin/env python3
"""Infer region trees from manifest-selected genomes without curated annotations."""

import argparse
from pathlib import Path
import sys

from extract_regions import ROOT, fasta_records, require, table
from run_phylogeny import (configured, replay_settings, resolve_tools, sha256,
                           workflow, workflow_lock)


def prepare_inputs(config, input_dir):
    sources = dict(sequences=ROOT / "data/sequences.csv", genomes=ROOT / "data/genomes.tsv",
                   regions=input_dir / "regions.tsv",
                   extracted_fasta=input_dir / "extracted_29_sequences.fasta",
                   raw_workflow=Path(__file__).resolve())
    sequences, genomes = table(sources["sequences"], ","), table(sources["genomes"])
    regions = table(sources["regions"])
    records = list(fasta_records(sources["extracted_fasta"]))
    order = [row["seq_id"] for row in sequences]
    require(len(order) == config["sequences"] and len(set(order)) == len(order),
            "Expected 29 unique sequence IDs")
    require(len(records) == len(regions) == len(genomes) == len(order),
            "FASTA, region, genome and sequence counts differ")
    require([row["seq_id"] for row in regions] == order, "Region order differs from sequences.csv")
    index = {row["seq_id"]: row for row in genomes}
    require(len(index) == len(order) and set(index) == set(order), "Genome IDs differ")
    normalized, oriented, species = [], [], []
    complement = str.maketrans("ACGTN", "TGCAN")
    for rank, ((header, dna), sequence, region) in enumerate(zip(records, sequences, regions), 1):
        identifier, dna = sequence["seq_id"], dna.upper()
        genome = index[identifier]
        require(header == region["fasta_header"] == region["contig"] + "_" + region["target_protein_id"],
                f"FASTA/region header differs: {identifier}")
        require(len(dna) == int(sequence["length"]) == int(region["length"])
                == int(genome["region_length"]), f"Region length differs: {identifier}")
        require(dna and set(dna) <= set("ACGTN"), f"Unexpected DNA alphabet: {identifier}")
        require(region["assembly_accession"] == genome["assembly_accession"],
                f"Region/genome assembly differs: {identifier}")
        require(region["target_strand"] in {"+", "-"} and region["region_strand"] == "+",
                f"Unexpected source orientation: {identifier}")
        reverse = region["target_strand"] == "-"
        normalized.append((identifier, dna))
        oriented.append((identifier, dna.translate(complement)[::-1] if reverse else dna))
        species.append(dict(rank=rank, seq_id=identifier, organism_name=genome["organism_name"],
                            strain=genome["strain"], assembly_accession=region["assembly_accession"],
                            assembly_status=genome["assembly_status"], genome_notes=genome.get("genome_notes", ""),
                            contig=region["contig"], target_strand=region["target_strand"],
                            target_protein_id=region["target_protein_id"], fasta_header=header,
                            length=len(dna), reverse_complemented_in_oriented_branch=reverse))
    require(sum(row["reverse_complemented_in_oriented_branch"] for row in species)
            == config["negative_target_strands"], "Expected 12 negative-strand targets")
    return dict(default=normalized, oriented=oriented), species, {key: sha256(path) for key, path in sources.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("default", "replay"), default="default")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "results/raw")
    parser.add_argument("--output-root", type=Path, default=ROOT / "results/raw/phylogeny")
    parser.add_argument("--replay-manifest", type=Path, default=ROOT / "config/phylogeny_replay.json")
    parser.add_argument("--mafft", help="MAFFT 7.526 executable")
    parser.add_argument("--iqtree", help="IQ-TREE 3.1.4 executable")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    config = configured()
    tools, tool_info = resolve_tools(config, args.mafft, args.iqtree)
    resolve = lambda path: (path if path.is_absolute() else ROOT / path).resolve()
    inputs, species, sources = prepare_inputs(config, resolve(args.input_dir))
    replays = replay_settings(resolve(args.replay_manifest), args.mode)
    output = resolve(args.output_root)
    require(output != ROOT and ROOT / "data" not in (output, *output.parents),
            "Output must not overwrite repository inputs")
    with workflow_lock(output):
        workflow(output, inputs, species, sources, replays, args.mode, config, tools, tool_info)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"Raw phylogeny failed: {error}")
