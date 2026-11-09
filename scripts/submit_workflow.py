#!/usr/bin/env python3
"""Submit one sequential 96-CPU workflow to node2 and wait for its result."""

import argparse
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def check_node():
    if (not os.environ.get("SLURM_JOB_ID") or socket.gethostname().split(".")[0] != "node2"
            or os.environ.get("SLURM_CPUS_PER_TASK") != "96"
            or os.environ.get("SLURM_NTASKS", "1") != "1"):
        raise ValueError("Run inside a node2 Slurm allocation with one task and 96 CPUs; "
                         "use make submit TARGET=verify (or plot-reference / plot-inferred).")


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
    parser.add_argument("--target", choices=("verify", "run", "phylogeny", "plot-reference", "plot-inferred"),
                        default="verify")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--bedtools", default=".tools/bedtools")
    for name in ("blastn", "makeblastdb", "mafft", "iqtree"):
        parser.add_argument("--" + name, default="")
    args = parser.parse_args()
    if args.check_node:
        check_node()
        print("Resource check passed: node2, one task, 96 CPUs.")
        return
    if args.mode == "replay" and not (ROOT / "config/phylogeny_replay.json").is_file():
        raise ValueError("Recorded seeds are missing; complete the initial MODE=default run first.")
    output = ROOT / "results/slurm"
    output.mkdir(parents=True, exist_ok=True)
    command = ["sbatch", "--parsable", "--wait", "--partition=normal", "--nodelist=node2",
               "--nodes=1", "--ntasks=1", "--cpus-per-task=96", "--job-name=cd-uptake-synteny",
               "--chdir=" + str(ROOT), "--output=" + str(output / "%j.log"),
               str(ROOT / "scripts/run_node2.sh"), args.target, args.mode, rscript_path(args.rscript),
               args.bedtools, args.blastn, args.makeblastdb, args.mafft, args.iqtree]
    print("Submitting to node2 (96 CPUs). Logs: results/slurm/<job_id>.log", flush=True)
    process = subprocess.run(command, cwd=ROOT)
    if process.returncode:
        raise ValueError(f"Slurm workflow failed (exit {process.returncode}); inspect results/slurm/.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"error: {error}")
