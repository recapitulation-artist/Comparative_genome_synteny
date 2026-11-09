#!/usr/bin/env python3
"""Submit a configurable single-node workflow to Slurm, or show its command."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from runtime import DEFAULT_THREADS


ROOT = Path(__file__).resolve().parents[1]
GENE_PALETTES = ("main-muted", "automatic", "nature", "cell", "science")
LINK_PALETTES = ("blue", "green", "purple", "orange", "original")
TARGETS = ("verify", "run", "plot-main", "plot-curated", "plot-auto",
           "verify-main", "verify-curated", "verify-auto", "prepare-curated", "prepare-auto",
           "blast-curated", "blast-auto", "phylogeny-curated", "phylogeny-auto")


def validate_link_colors(value):
    if value == "":
        return
    colors = [color.strip() for color in value.split(",")]
    if (len(colors) != 3 or any(re.fullmatch(r"#[0-9A-Fa-f]{6}", color) is None for color in colors)
            or len({color.upper() for color in colors}) != 3):
        raise ValueError("Provide three distinct comma-separated #RRGGBB link colors, "
                         "ordered from high to low identity.")


def positive_threads(value):
    """Validate a future allocation without inspecting the login node's CPUs."""
    if isinstance(value, bool) or not str(value).isdigit() or int(value) < 1:
        raise ValueError("THREADS must be a positive integer")
    return int(value)


def make_assignment(name, value):
    # GNU Make expands command-line variables; protect literal dollar signs.
    return name + "=" + str(value).replace("$", "$$")


def argument_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("default", "replay"), default="default")
    parser.add_argument("--target", choices=TARGETS, default="verify")
    parser.add_argument("--threads", default=DEFAULT_THREADS,
                        help=f"CPUs per task and workflow CPU budget (default: {DEFAULT_THREADS})")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--target-gene", default="Nramp")
    parser.add_argument("--neighbors", type=int, default=10)
    parser.add_argument("--genome-dir", default="data/genomes")
    parser.add_argument("--genome-metadata", default="data/genomes.tsv")
    parser.add_argument("--annotations", default="data/genes.tsv")
    parser.add_argument("--gene-palette", choices=GENE_PALETTES, default="main-muted")
    parser.add_argument("--link-palette", choices=LINK_PALETTES, default="blue")
    parser.add_argument("--link-colors", default="",
                        help="Three distinct #RRGGBB colors (high, medium, low identity); overrides the preset")
    parser.add_argument("--rscript", default="", help="Optional Rscript override; otherwise use repository/PATH lookup")
    for name in ("bedtools", "blastn", "makeblastdb", "mafft", "iqtree"):
        parser.add_argument("--" + name, default="")
    for name in ("partition", "node", "time", "memory"):
        parser.add_argument("--" + name, default="", help="Optional Slurm " + name)
    parser.add_argument("--dry-run", action="store_true", help="Print JSON argv without submitting or writing files")
    return parser


def build_command(args, root=ROOT):
    threads = positive_threads(args.threads)
    if args.neighbors < 0 or not args.target_gene.strip():
        raise ValueError("Provide a nonempty target gene and a nonnegative neighbor count.")
    validate_link_colors(args.link_colors)
    root = Path(root).resolve()
    command = ["sbatch", "--parsable", "--wait", "--nodes=1", "--ntasks=1",
               f"--cpus-per-task={threads}", "--job-name=cd-uptake-synteny",
               "--chdir=" + str(root), "--output=" + str(root / "results/slurm/%j.log")]
    for attribute, option in (("partition", "partition"), ("node", "nodelist"),
                              ("time", "time"), ("memory", "mem")):
        value = getattr(args, attribute)
        if value:
            command.append(f"--{option}={value}")
    command.extend([str(root / "scripts/run_slurm.sh"), str(threads), args.python, args.target])
    options = dict(MODE=args.mode, THREADS=threads, PYTHON=args.python,
                   BEDTOOLS=args.bedtools, BLASTN=args.blastn, MAKEBLASTDB=args.makeblastdb,
                   MAFFT=args.mafft, IQTREE=args.iqtree, TARGET_GENE=args.target_gene,
                   NEIGHBORS=args.neighbors, GENOME_DIR=args.genome_dir,
                   GENOME_METADATA=args.genome_metadata, GENE_ANNOTATIONS=args.annotations,
                   GENE_PALETTE=args.gene_palette, LINK_PALETTE=args.link_palette,
                   LINK_COLORS=args.link_colors)
    if args.rscript:
        options["RSCRIPT"] = args.rscript
    command.extend(make_assignment(name, value) for name, value in options.items())
    return command


def submission_environment():
    # An sbatch invoked inside another job must not inherit its resource defaults.
    return {name: value for name, value in os.environ.items() if not name.startswith("SBATCH_")}


def main():
    args = argument_parser().parse_args()
    command = build_command(args)
    if args.dry_run:
        print(json.dumps(dict(argv=command, cwd=str(ROOT), threads=positive_threads(args.threads)), indent=2))
        return
    (ROOT / "results/slurm").mkdir(parents=True, exist_ok=True)
    print(f"Submitting one Slurm task with {positive_threads(args.threads)} CPUs. "
          "Logs: results/slurm/<job_id>.log", flush=True)
    process = subprocess.run(command, cwd=ROOT, env=submission_environment())
    if process.returncode:
        raise ValueError(f"Slurm workflow failed (exit {process.returncode}); inspect results/slurm/.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"error: {error}")
