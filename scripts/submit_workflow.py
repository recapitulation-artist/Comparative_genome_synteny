#!/usr/bin/env python3
"""Submit one sequential 96-CPU workflow to node2 and wait for its result."""

import argparse
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
GENE_PALETTES = ("main-muted", "automatic", "nature", "cell", "science")
LINK_PALETTES = ("blue", "green", "purple", "orange", "original")


def validate_link_colors(value):
    if value == "":
        return
    colors = [color.strip() for color in value.split(",")]
    if (len(colors) != 3 or any(re.fullmatch(r"#[0-9A-Fa-f]{6}", color) is None for color in colors)
            or len({color.upper() for color in colors}) != 3):
        raise ValueError("Provide three distinct comma-separated #RRGGBB link colors, "
                         "ordered from high to low identity.")


def check_node():
    if (not os.environ.get("SLURM_JOB_ID") or socket.gethostname().split(".")[0] != "node2"
            or os.environ.get("SLURM_CPUS_PER_TASK") != "96"
            or os.environ.get("SLURM_NTASKS", "1") != "1"):
        raise ValueError("Run inside a node2 Slurm allocation with one task and 96 CPUs; "
                         "use make submit TARGET=verify, plot-main, plot-curated or plot-auto.")


def rscript_path(value):
    candidate = shutil.which(value)
    if candidate:
        return str(Path(candidate).resolve())
    if value != "Rscript":
        raise ValueError(f"Rscript executable not found: {value}")
    for prefix in (ROOT / ".tools/r", Path.home() / "miniforge3/envs/r_env"):
        candidate = prefix / "bin/Rscript"
        if candidate.is_file():
            return str(candidate.resolve())
    raise ValueError("Rscript is missing; activate the locked R environment or set RSCRIPT=/path/to/Rscript")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-node", action="store_true")
    parser.add_argument("--mode", choices=("default", "replay"), default="default")
    parser.add_argument("--target", choices=("verify", "run", "plot-main", "plot-curated", "plot-auto",
                                           "verify-main", "verify-curated", "verify-auto", "prepare-curated", "prepare-auto",
                                           "blast-curated", "blast-auto", "phylogeny-curated", "phylogeny-auto",
                                           "phylogeny", "plot-reference", "plot-inferred", "extract-raw", "blast-raw",
                                           "phylogeny-raw", "plot-raw-inferred", "verify-raw"),
                        default="verify")
    parser.add_argument("--target-gene", default="Nramp")
    parser.add_argument("--neighbors", type=int, default=10)
    parser.add_argument("--genome-dir", default="data/genomes")
    parser.add_argument("--genome-metadata", default="data/genomes.tsv")
    parser.add_argument("--annotations", default="data/genes.tsv")
    parser.add_argument("--gene-palette", choices=GENE_PALETTES, default="main-muted")
    parser.add_argument("--link-palette", choices=LINK_PALETTES, default="blue")
    parser.add_argument("--link-colors", default="",
                        help="Three distinct #RRGGBB colors (high, medium, low identity); overrides the preset")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--bedtools", default=".tools/bedtools")
    for name in ("blastn", "makeblastdb", "mafft", "iqtree"):
        parser.add_argument("--" + name, default="")
    args = parser.parse_args()
    if args.check_node:
        check_node()
        print("Resource check passed: node2, one task, 96 CPUs.")
        return
    if args.neighbors < 0 or not args.target_gene.strip():
        raise ValueError("Provide a nonempty target gene and a nonnegative neighbor count.")
    validate_link_colors(args.link_colors)
    output = ROOT / "results/slurm"
    output.mkdir(parents=True, exist_ok=True)
    command = ["sbatch", "--parsable", "--wait", "--partition=normal", "--nodelist=node2",
               "--nodes=1", "--ntasks=1", "--cpus-per-task=96", "--job-name=cd-uptake-synteny",
               "--chdir=" + str(ROOT), "--output=" + str(output / "%j.log"),
               str(ROOT / "scripts/run_node2.sh"), args.target, args.mode, rscript_path(args.rscript),
               args.bedtools, args.blastn, args.makeblastdb, args.mafft, args.iqtree,
               args.target_gene, str(args.neighbors), args.genome_dir, args.genome_metadata, args.annotations,
               args.gene_palette, args.link_palette, args.link_colors]
    print("Submitting to node2 (96 CPUs). Logs: results/slurm/<job_id>.log", flush=True)
    process = subprocess.run(command, cwd=ROOT)
    if process.returncode:
        raise ValueError(f"Slurm workflow failed (exit {process.returncode}); inspect results/slurm/.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"error: {error}")
