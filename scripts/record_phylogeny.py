#!/usr/bin/env python3
"""Record first-run seeds and hashes, or prepare recorded trees for replay checks."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

from extract_regions import ROOT, require


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phylogeny-dir", type=Path, default=ROOT / "results/phylogeny")
    parser.add_argument("--manifest", type=Path, default=ROOT / "config/phylogeny_replay.json")
    parser.add_argument("--prepare-replay", action="store_true")
    args = parser.parse_args()
    directory = args.phylogeny_dir.resolve()
    manifest = args.manifest.resolve()
    if args.prepare_replay:
        data = json.loads(manifest.read_text())
        require(data["schema_version"] == 1, "Unsupported replay manifest")
        for branch in ("default", "oriented"):
            record = data["branches"][branch]
            payload = record["tree_newick"].encode("ascii")
            require(digest(payload) == record["tree_sha256"], "Recorded tree hash differs: " + branch)
            (directory / ("recorded_" + branch + ".nwk")).write_bytes(payload)
        return
    data = dict(schema_version=1, branches={})
    for branch in ("default", "oriented"):
        source = directory / branch
        run = json.loads((source / "run.json").read_text())
        require(run["status"] == "completed" and run["mode"] == "default" and run["threads"] == 96,
                "Record only completed initial default-mode runs")
        for filename, key in (("input.fasta", "input_sha256"), ("aligned.fasta", "aligned_sha256"),
                              ("tree.treefile", "tree_sha256")):
            require(digest((source / filename).read_bytes()) == run[key], "Run output hash differs: " + filename)
        data["branches"][branch] = {key: run[key] for key in
                                   ("seed", "input_sha256", "aligned_sha256", "tree_sha256", "selected_model")}
        data["branches"][branch]["tree_newick"] = (source / "tree.treefile").read_text(encoding="ascii")
    if manifest.exists():
        require(json.loads(manifest.read_text()) == data,
                "Existing recorded baseline differs; select a new --manifest path")
    else:
        manifest.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print("Recorded first-run seeds and reference tree hashes: " + str(manifest))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"error: {error}")
