#!/usr/bin/env python3
"""Install or verify the locked official BEDTools binary without system changes."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
LOCK = json.loads((ROOT / "environment/bedtools.lock.json").read_text())


def verify(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != LOCK["sha256"]:
        raise ValueError(f"BEDTools SHA256 mismatch: {path}")
    version = subprocess.check_output([str(path), "--version"], text=True).strip()
    if version != f"bedtools v{LOCK['version']}":
        raise ValueError(f"Unexpected BEDTools version: {version}")
    print(f"{version}: SHA256 verified")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify without downloading")
    parser.add_argument("--bedtools", type=Path, default=ROOT / ".tools/bedtools")
    args = parser.parse_args()
    if sys.version_info < (3, 9):
        raise ValueError("Python >=3.9 is required")
    path = args.bedtools.resolve()
    if args.check or path.exists():
        verify(path)
        return
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "AMD64"):
        raise ValueError("The locked BEDTools binary requires Linux x86_64")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent) as temporary:
        downloaded = Path(temporary) / "bedtools"
        subprocess.run([
            "curl", "--fail", "--location", "--retry", "2",
            "--connect-timeout", "15", "--max-time", "300",
            "--output", str(downloaded), LOCK["url"],
        ], check=True)
        # Check the bytes before executing the downloaded program.
        if hashlib.sha256(downloaded.read_bytes()).hexdigest() != LOCK["sha256"]:
            raise ValueError("Downloaded BEDTools SHA256 mismatch")
        downloaded.chmod(0o755)
        verify(downloaded)
        os.replace(downloaded, path)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.exit(f"BEDTools setup failed: {error}")
