#!/usr/bin/env python3
"""Infer a region tree from one prepared entry with local MAFFT and IQ-TREE."""

import argparse
import json
import os
from pathlib import Path
import re
import sys

from extract_regions import ROOT, fasta_bytes, fasta_records, require, table, write_table
from run_phylogeny import (CONFIG_PATH, cached_state, configured, execute_branch, payload_hash,
                           resolve_tools, save_json, scheduler_context, sha256, workflow_lock)


def prepare_inputs(input_dir):
    sources = dict(sequences=input_dir / "sequences.csv", regions=input_dir / "regions.tsv",
                   genomes=input_dir / "genomes.tsv", extracted_fasta=input_dir / "extracted_sequences.fasta")
    sequences = table(sources["sequences"], ",")
    regions, genomes = table(sources["regions"]), table(sources["genomes"])
    records = list(fasta_records(sources["extracted_fasta"]))
    order = [row["seq_id"] for row in sequences]
    require(len(order) >= 4, f"Tree inference requires at least four regions; found {len(order)}")
    require(len(set(order)) == len(order) and all(re.fullmatch(r"[A-Za-z0-9_.-]+", name) for name in order),
            "Tree sequence IDs must be unique and contain only letters, digits, underscores, dots and hyphens")
    require(len(regions) == len(genomes) == len(records) == len(order), "Sequence, region, genome and FASTA counts differ")
    region_index, genome_index = ({row["seq_id"]: row for row in rows} for rows in (regions, genomes))
    require(len(region_index) == len(genome_index) == len(order)
            and set(region_index) == set(genome_index) == set(order), "Entry table sequence IDs differ")
    fasta_index = dict(records)
    require(len(fasta_index) == len(records) and len({row["fasta_header"] for row in regions}) == len(order)
            and set(fasta_index) == {row["fasta_header"] for row in regions}, "FASTA headers differ from regions.tsv")
    normalized, species = [], []
    for rank, sequence in enumerate(sequences, 1):
        identifier = sequence["seq_id"]
        region, genome = region_index[identifier], genome_index[identifier]
        header = region["fasta_header"]
        dna = fasta_index[header].upper()
        require(dna and set(dna) <= set("ACGTN"), f"Invalid DNA sequence: {identifier}")
        require(len(dna) == int(sequence["length"]) == int(region["length"]), f"Region length differs: {identifier}")
        if genome.get("region_length"):
            require(len(dna) == int(genome["region_length"]), f"Genome metadata region length differs: {identifier}")
        require(region["target_strand"] in {"+", "-"} and region["region_strand"] == "+",
                f"Invalid target or extracted-region strand: {identifier}")
        require(region["assembly_accession"] == genome["assembly_accession"], f"Region assembly differs: {identifier}")
        require(bool(genome["organism_name"].strip()), f"Missing species name: {identifier}")
        normalized.append((identifier, dna))
        species.append(dict(rank=rank, seq_id=identifier, organism_name=genome["organism_name"],
                            strain=genome.get("strain", ""), assembly_accession=region["assembly_accession"],
                            assembly_status=genome.get("assembly_status", ""), genome_notes=genome.get("genome_notes", ""),
                            contig=region["contig"], target_strand=region["target_strand"],
                            target_protein_id=region["target_protein_id"],
                            target_locus_tag=region.get("target_locus_tag", ""),
                            target_product=region.get("target_product", ""), fasta_header=header, length=len(dna),
                            reverse_complemented=False))
    return normalized, species, sources


def load_replay(path, entry, payload):
    manifest = json.loads(path.read_text(encoding="utf-8"))
    require(manifest.get("schema_version") == 1 and manifest.get("entry") == entry
            and set(manifest.get("branches", {})) == {"default"}, "Replay manifest belongs to a different entry or schema")
    replay = manifest["branches"]["default"]
    require(isinstance(replay.get("seed"), int) and not isinstance(replay["seed"], bool)
            and replay["seed"] > 0, "Invalid entry replay seed")
    for key in ("input_sha256", "aligned_sha256", "tree_sha256"):
        require(isinstance(replay.get(key), str) and re.fullmatch(r"[0-9a-f]{64}", replay[key]) is not None,
                f"Invalid entry replay {key}")
    require(payload_hash(payload) == replay["input_sha256"], "Entry replay input differs from the recorded run")
    return replay


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entry", choices=("curated", "auto"), required=True)
    parser.add_argument("--input-dir", type=Path, help="Prepared entry directory (default: results/<entry>)")
    parser.add_argument("--output-root", type=Path, help="Tree directory (default: entry/phylogeny or entry/phylogeny_replay)")
    parser.add_argument("--mode", choices=("default", "replay"), default="default")
    parser.add_argument("--replay-manifest", type=Path, help="Recorded run (default: entry/phylogeny/replay_manifest.json)")
    parser.add_argument("--mafft", help="MAFFT 7.526 executable")
    parser.add_argument("--iqtree", help="IQ-TREE 3.1.4 executable")
    parser.add_argument("--check-inputs", action="store_true", help="Validate prepared inputs without running tools")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    resolve = lambda path: (path if path.is_absolute() else ROOT / path).resolve()
    input_dir = resolve(args.input_dir) if args.input_dir else ROOT / "results" / args.entry
    records, species, sources = prepare_inputs(input_dir)
    negative_count = sum(row["target_strand"] == "-" for row in species)
    if args.check_inputs:
        print(f"Validated {len(records)} {args.entry} regions, including {negative_count} negative-strand targets.")
        return
    output = resolve(args.output_root) if args.output_root else input_dir / (
        "phylogeny_replay" if args.mode == "replay" else "phylogeny")
    require(output != ROOT and ROOT / "data" not in (output, *output.parents), "Tree output must not overwrite repository inputs")
    payload = fasta_bytes(records)
    replay_path = resolve(args.replay_manifest) if args.replay_manifest else input_dir / "phylogeny/replay_manifest.json"
    replay = load_replay(replay_path, args.entry, payload) if args.mode == "replay" else None
    config = configured()
    config.update(sequences=len(records), negative_target_strands=negative_count)
    tools, tool_info = resolve_tools(config, args.mafft, args.iqtree)
    sources.update(entry_script=Path(__file__).resolve(), phylogeny_script=ROOT / "scripts/run_phylogeny.py",
                   input_helpers=ROOT / "scripts/extract_regions.py", config=CONFIG_PATH,
                   runtime_lock=ROOT / config["runtime_lock"])
    source_hashes = {key: sha256(path) for key, path in sources.items()}
    signature = payload_hash(json.dumps(dict(schema_version=1, entry=args.entry, branch="default", mode=args.mode,
                                         input_sha256=payload_hash(payload), source_sha256=source_hashes,
                                         tools={name: dict(version=info["version"], sha256=info["executable_sha256"])
                                                for name, info in tool_info.items()}, replay=replay), sort_keys=True).encode("utf-8"))
    with workflow_lock(output):
        directory = output / "default"
        state = cached_state(directory, signature, payload, records, replay, 0, config)
        needs_compute = state is None or state.get("status") != "completed"
        if needs_compute:
            require(os.environ.get("SLURM_NTASKS", "1") == "1", "Tree inference requires one Slurm task")
        execution = scheduler_context(config) if needs_compute else None
        species_path = output / "species.tsv"
        if species_path.exists():
            expected = [{key: str(value).lower() if isinstance(value, bool) else str(value)
                         for key, value in row.items()} for row in species]
            require(table(species_path) == expected, "Cached species metadata differs from current inputs")
        else:
            write_table(species_path, species)
        execute_branch(directory, "default", records, payload, signature, state, tools, tool_info,
                       config, args.mode, replay, execution, source_hashes)
        for key, path in sources.items():
            require(sha256(path) == source_hashes[key], f"Tree input changed during execution: {path}")
        state = json.loads((directory / "run.json").read_text(encoding="utf-8"))
        state.update(entry=args.entry, input_sequence_count=len(records), negative_target_strands=negative_count)
        save_json(directory / "run.json", state)
        provenance = dict(schema_version=1, entry=args.entry, mode=args.mode, input_sequences=len(records),
                          negative_target_strands=negative_count, threads=config["threads"],
                          inputs={key: dict(path=str(path.relative_to(ROOT)) if ROOT in path.parents else str(path),
                                            sha256=source_hashes[key]) for key, path in sources.items()},
                          tools=tool_info)
        save_json(output / "provenance.json", provenance)
        actual = {key: state[key] for key in ("seed", "input_sha256", "aligned_sha256", "tree_sha256")}
        save_json(output / "replay_manifest.json", dict(schema_version=1, entry=args.entry, branches=dict(default=actual)))
    print(f"{args.entry} tree: {directory / 'tree.treefile'}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"Entry phylogeny failed: {error}")
