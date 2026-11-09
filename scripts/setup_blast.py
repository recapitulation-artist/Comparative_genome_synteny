#!/usr/bin/env python3
"""Use a matching local BLAST installation or restore the locked Conda runtime."""

from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def main():
    check = [sys.executable, str(ROOT / "scripts/run_blast.py"), "--check-tools"]
    result = subprocess.run(check, capture_output=True, text=True)
    if result.returncode == 0:
        print(result.stdout.strip())
        return
    conda = shutil.which("conda")
    if conda is None:
        sys.exit("BLAST+ 2.17.0+ is required. Install Conda to restore its locked runtime.")
    prefix = ROOT / ".tools/blast"
    if prefix.exists():
        sys.exit(f"Existing BLAST environment failed validation: {result.stderr.strip()}")
    subprocess.run([conda, "create", "--yes", "--prefix", str(prefix),
                    "--file", str(ROOT / "environment/blast_linux_64.lock")], check=True)
    subprocess.run(check, check=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, subprocess.CalledProcessError) as error:
        sys.exit(f"BLAST setup failed: {error}")
