#!/usr/bin/env python3
"""Verify an installed BEDTools or download the checksum-locked official binary."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
LOCK = json.loads((ROOT / "environment/bedtools.lock.json").read_text())


def verify(path, official=False):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if official and digest != LOCK["sha256"]:
        raise ValueError(f"Official BEDTools SHA256 mismatch: {path}")
    version = subprocess.check_output([str(path), "--version"], text=True).strip()
    if version != f"bedtools v{LOCK['version']}":
        raise ValueError(f"Unexpected BEDTools version: {version}")
    print(f"{version}: {path}; SHA256={digest}" + (" (official checksum verified)" if official else ""))
    return dict(version=version, executable=str(path), executable_sha256=digest,
                official_checksum_verified=official)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify without downloading")
    parser.add_argument("--bedtools", help="Explicit executable; otherwise .tools/bedtools then PATH")
    args = parser.parse_args()
    if sys.version_info < (3, 9):
        raise ValueError("Python >=3.9 is required")
    local = ROOT / ".tools/bedtools"
    candidate = args.bedtools or (str(local) if local.is_file() else shutil.which("bedtools"))
    if candidate:
        path = Path(candidate).expanduser()
        if not path.is_absolute():
            path = ROOT / path if (ROOT / path).is_file() else Path(shutil.which(str(path)) or path)
        path = path.resolve()
        if not path.is_file() or not os.access(path, os.X_OK):
            raise ValueError(f"BEDTools is not executable: {candidate}")
        # A user-installed build can differ in bytes while using the locked version.
        provenance = ROOT / ".tools/bedtools.provenance.json"
        previous = json.loads(provenance.read_text()) if provenance.is_file() else {}
        official = previous.get("official_checksum_verified", False) and previous.get("executable") == str(path)
        record = verify(path, official=official)
    else:
        if args.check:
            raise ValueError("BEDTools is missing; run make setup-bedtools or specify --bedtools PATH")
        if platform.system() != "Linux" or platform.machine() not in ("x86_64", "AMD64"):
            raise ValueError("The locked BEDTools binary requires Linux x86_64")
        local.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=local.parent) as temporary:
            downloaded = Path(temporary) / "bedtools"
            subprocess.run([
                "curl", "--fail", "--location", "--retry", "2",
                "--connect-timeout", "15", "--max-time", "300",
                "--output", str(downloaded), LOCK["url"],
            ], check=True)
            # Verify the downloaded bytes before executing the program.
            if hashlib.sha256(downloaded.read_bytes()).hexdigest() != LOCK["sha256"]:
                raise ValueError("Downloaded BEDTools SHA256 mismatch")
            downloaded.chmod(0o755)
            record = verify(downloaded, official=True)
            os.replace(downloaded, local)
        record["executable"] = str(local)
    if not args.check:
        (ROOT / ".tools").mkdir(exist_ok=True)
        (ROOT / ".tools/bedtools.provenance.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.exit(f"BEDTools setup failed: {error}")
