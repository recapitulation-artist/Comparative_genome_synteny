#!/usr/bin/env python3
"""Find matching local phylogeny tools or restore their small locked runtime."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

from run_phylogeny import ROOT, configured, resolve_tools


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check local tools without installation")
    parser.add_argument("--prefix", type=Path, default=ROOT / ".tools/phylogeny")
    parser.add_argument("--mafft", help="Explicit MAFFT executable")
    parser.add_argument("--iqtree", help="Explicit IQ-TREE executable")
    args = parser.parse_args()
    config = configured()
    prefix = args.prefix if args.prefix.is_absolute() else ROOT / args.prefix
    mafft = args.mafft or (str(prefix / "bin/mafft") if (prefix / "bin/mafft").is_file() else None)
    iqtree = args.iqtree or (str(prefix / "bin/iqtree3") if (prefix / "bin/iqtree3").is_file() else None)
    try:
        _, info = resolve_tools(config, mafft, iqtree)
    except ValueError:
        if args.check or args.mafft or args.iqtree:
            raise
        if prefix.exists():
            raise ValueError(f"Existing phylogeny environment failed validation: {prefix}")
        conda = shutil.which("conda")
        if conda is None:
            candidate = Path.home() / "miniforge3/bin/conda"
            conda = str(candidate) if candidate.is_file() else None
        if conda is None:
            raise ValueError("Install Conda to restore MAFFT 7.526 and IQ-TREE 3.1.4")
        subprocess.run([conda, "create", "--yes", "--prefix", str(prefix),
                        "--file", str(ROOT / config["runtime_lock"])], check=True)
        _, info = resolve_tools(config, str(prefix / "bin/mafft"), str(prefix / "bin/iqtree3"))
    for name, record in info.items():
        print(f"{name} {record['version']}: {record['executable']}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(f"Phylogeny setup failed: {error}")
