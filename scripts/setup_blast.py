#!/usr/bin/env python3
"""Use matching BLAST executables or restore the locked Conda runtime."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blastn", help="Explicit blastn executable")
    parser.add_argument("--makeblastdb", help="Explicit makeblastdb executable")
    parser.add_argument("--check", action="store_true", help="Check without installing")
    args = parser.parse_args()
    check = [sys.executable, str(ROOT / "scripts/run_blast.py"), "--check-tools"]
    for name in ("blastn", "makeblastdb"):
        if getattr(args, name):
            check.extend(["--" + name, getattr(args, name)])
    result = subprocess.run(check, capture_output=True, text=True)
    if result.returncode == 0:
        print(result.stdout.strip())
        return
    if args.check or args.blastn or args.makeblastdb:
        raise ValueError(result.stderr.strip() or result.stdout.strip())
    conda = shutil.which("conda")
    if conda is None:
        raise ValueError("BLAST+ 2.17.0+ is required. Install Conda to restore its locked runtime.")
    prefix = ROOT / ".tools/blast"
    if prefix.exists():
        raise ValueError(f"Existing BLAST environment failed validation: {result.stderr.strip()}")
    subprocess.run([conda, "create", "--yes", "--prefix", str(prefix),
                    "--file", str(ROOT / "environment/blast_linux_64.lock")], check=True)
    subprocess.run(check, check=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.exit(f"BLAST setup failed: {error}")
