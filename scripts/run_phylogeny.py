#!/usr/bin/env python3
"""Shared MAFFT/IQ-TREE execution, resumable state and provenance helpers."""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from extract_regions import ROOT, fasta_records, require, write_table
from runtime import DEFAULT_THREADS, execution_context, thread_environment


CONFIG_PATH = ROOT / "config/phylogeny.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload_hash(payload):
    return hashlib.sha256(payload).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def save_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def configured(threads=None):
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    require(config["schema_version"] == 1, "Unsupported phylogeny configuration schema")
    value = config.get("threads", DEFAULT_THREADS) if threads is None else threads
    require(not isinstance(value, bool) and str(value).isdigit() and int(value) > 0,
            "Phylogeny threads must be a positive integer")
    config["threads"] = int(value)
    # Preserve the scientific defaults; these flags specify resources only.
    config["mafft_options"] = ["--thread", str(value), "--threadtb", str(value)]
    config["iqtree_options"] = ["-T", str(value)]
    return config


def resolve_tool(name, requested, version):
    candidates = [requested] if requested else [
        str(ROOT / ".tools/phylogeny/bin" / name),
        shutil.which(name),
    ]
    problems = []
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if not path.is_absolute():
            path = ROOT / path if (ROOT / path).is_file() else Path(shutil.which(str(path)) or path)
        if not path.is_file() or not os.access(path, os.X_OK):
            problems.append(f"Not executable: {candidate}")
            continue
        path = path.resolve()
        flag = "--version"
        result = subprocess.run([str(path), flag], capture_output=True, text=True)
        text = (result.stdout + result.stderr).strip()
        pattern = (rf"^v{re.escape(version)}(?:\s|$)" if name == "mafft"
                   else rf"^IQ-TREE version {re.escape(version)}(?:\s|$)")
        if result.returncode == 0 and re.search(pattern, text, re.MULTILINE):
            return path, dict(version=version, version_output=text, executable_sha256=sha256(path),
                              executable=str(path))
        problems.append(f"Unexpected {name} version at {path}: {text.splitlines()[:1]}")
    raise ValueError(f"{name} {version} is required; run make setup-phylogeny. " + "; ".join(problems))


def resolve_tools(config, mafft=None, iqtree=None):
    alignment, alignment_info = resolve_tool("mafft", mafft, config["mafft_version"])
    inference, inference_info = resolve_tool("iqtree3", iqtree, config["iqtree_version"])
    return dict(mafft=alignment, iqtree=inference), dict(mafft=alignment_info, iqtree=inference_info)


def scheduler_context(config):
    """Return the same execution metadata for local and allocated runs."""
    return execution_context(config["threads"])


def portable_commands(config, replay):
    inference = ["iqtree3", "-s", "aligned.fasta", *config["iqtree_options"], "--prefix", "tree"]
    if replay is not None:
        inference.extend(["--seed", str(replay["seed"])])
    return dict(mafft=["mafft", *config["mafft_options"], "input.fasta"],
                mafft_stdout="aligned.fasta", mafft_stderr="mafft.log", iqtree=inference,
                iqtree_stdout="iqtree.stdout.log", iqtree_stderr="iqtree.stderr.log",
                environment=dict(LC_ALL="C", **thread_environment(config["threads"])),
                unset_environment=["MAFFT_BINARIES"])


def alignment_metrics(directory, records, reversed_count):
    aligned = list(fasta_records(directory / "aligned.fasta"))
    require([name for name, _ in aligned] == [name for name, _ in records],
            "Alignment must preserve all sequence IDs and their input order")
    lengths = {len(dna) for _, dna in aligned}
    require(len(lengths) == 1 and next(iter(lengths)) > 0, "Alignment rows have unequal or empty lengths")
    for (name, dna), (_, original) in zip(aligned, records):
        require(set(dna.upper()) <= set("ACGTN-"), f"Unexpected aligned DNA alphabet: {name}")
        require(dna.replace("-", "").upper() == original.upper(), f"Alignment changed input DNA: {name}")
    return dict(input_sequences=len(records), input_bases=sum(len(dna) for _, dna in records),
                aligned_sequences=len(aligned), alignment_columns=next(iter(lengths)),
                gap_characters=sum(dna.count("-") for _, dna in aligned),
                ungapped_input_matches=True, reverse_complemented_sequences=reversed_count,
                input_sha256=sha256(directory / "input.fasta"), aligned_sha256=sha256(directory / "aligned.fasta"))


def logged_threads(directory, expected):
    text = (directory / "mafft.log").read_text(encoding="utf-8", errors="replace")
    mafft = {}
    for name, count in re.findall(r"\b(nthread(?:pair|tb)?)\s*=\s*(\d+)", text):
        mafft.setdefault(name, set()).add(int(count))
    require(expected in mafft.get("nthread", set()) | mafft.get("nthreadpair", set()),
            f"MAFFT log does not confirm the {expected}-thread request")
    require(expected in mafft.get("nthreadtb", set()), f"MAFFT log does not confirm {expected} progressive-stage threads")
    texts = []
    for name in ("tree.log", "iqtree.stdout.log"):
        path = directory / name
        if path.is_file():
            texts.append(path.read_text(encoding="utf-8", errors="replace"))
    inference = {int(value) for value in re.findall(r"^Kernel:.*?\b(\d+)\s+threads?\b", "\n".join(texts), re.MULTILINE)}
    require(inference == {expected}, f"IQ-TREE log does not confirm exactly {expected} threads: {sorted(inference)}")
    return dict(mafft={key: sorted(values) for key, values in mafft.items()}, iqtree=sorted(inference))


def inference_details(directory, order):
    tree = (directory / "tree.treefile").read_text(encoding="utf-8").strip()
    require(tree.endswith(";"), "Incomplete IQ-TREE Newick output")
    tips = re.findall(r"(?:\(|,)\s*([A-Za-z0-9_.-]+)\s*(?=:|,|\))", tree)
    require(len(tips) == len(order) and len(set(tips)) == len(order) and set(tips) == set(order),
            "IQ-TREE tips differ from the input sequence IDs")
    report = (directory / "tree.iqtree").read_text(encoding="utf-8", errors="replace")
    seeds = set(re.findall(r"^Random seed number:\s*(\d+)", report, re.MULTILINE))
    if not seeds:
        for name in ("tree.log", "iqtree.stdout.log"):
            path = directory / name
            if path.is_file():
                seeds.update(re.findall(r"^Seed:\s*(\d+)", path.read_text(errors="replace"), re.MULTILINE))
    require(len(seeds) == 1, f"Could not determine one actual IQ-TREE random seed: {sorted(seeds)}")
    models = re.findall(r"^Model of substitution:\s*(\S+)", report, re.MULTILINE)
    if not models:
        models = re.findall(r"^Best-fit model according to BIC:\s*(\S+)", report, re.MULTILINE)
    require(len(set(models)) == 1, "Could not determine the selected substitution model")
    return dict(seed=int(next(iter(seeds))), selected_model=models[0],
                tree_sha256=sha256(directory / "tree.treefile"))


def completed_report_exists(directory):
    report = directory / "tree.iqtree"
    stdout = directory / "iqtree.stdout.log"
    if not report.is_file() or not stdout.is_file() or not (directory / "tree.treefile").is_file():
        return False
    text = stdout.read_text(errors="replace")
    return "Analysis results written" in text and "Total wall-clock time used" in report.read_text(errors="replace")


def expected_replay(replay, metrics, inference):
    if replay is None:
        return {}
    require(inference["seed"] == replay["seed"], "IQ-TREE did not use the recorded replay seed")
    require(metrics["aligned_sha256"] == replay["aligned_sha256"], "Replayed alignment SHA256 differs")
    if replay.get("input_sha256"):
        require(metrics["input_sha256"] == replay["input_sha256"], "Replay input SHA256 differs")
    return dict(exact_tree_match=inference["tree_sha256"] == replay["tree_sha256"],
                expected_tree_sha256=replay["tree_sha256"], expected_aligned_sha256=replay["aligned_sha256"])


def cached_state(directory, signature, payload, records, replay, reversed_count, config):
    path = directory / "run.json"
    if not path.is_file():
        require(not directory.exists() or not any(directory.iterdir()),
                f"Unrecognized existing output in {directory}; choose a new --output-root")
        return None
    state = json.loads(path.read_text(encoding="utf-8"))
    require(state.get("cache_signature") == signature,
            f"Input/tool/settings mismatch in {directory}; choose a new --output-root")
    require((directory / "input.fasta").is_file()
            and (directory / "input.fasta").read_bytes() == payload, f"Cached input changed: {directory}")
    require(state.get("phase") in {"prepared", "aligning", "aligned", "inferring", "inferred", "completed"},
            f"Unrecognized incomplete run phase in {directory}")
    if state["phase"] in {"aligned", "inferring", "inferred", "completed"}:
        require(sha256(directory / "aligned.fasta") == state.get("aligned_sha256"),
                f"Cached alignment changed: {directory}")
        alignment_metrics(directory, records, reversed_count)
    if state.get("status") == "completed":
        require(state["phase"] == "completed", f"Inconsistent completion state: {directory}")
        for name, digest in state["output_sha256"].items():
            require((directory / name).is_file() and sha256(directory / name) == digest,
                    f"Cached output changed: {directory / name}")
        metrics = alignment_metrics(directory, records, reversed_count)
        inference = inference_details(directory, [name for name, _ in records])
        require(inference["seed"] == state["seed"] and inference["selected_model"] == state["selected_model"],
                f"Cached inference provenance differs: {directory}")
        logged_threads(directory, config["threads"])
        comparison = expected_replay(replay, metrics, inference)
        for key, value in comparison.items():
            require(state.get(key) == value, f"Cached replay comparison differs: {directory}")
    return state


def run_command(command, directory, stdout, stderr, environment, append=False):
    with (directory / stdout).open("a" if append else "w", encoding="utf-8") as out, \
            (directory / stderr).open("a" if append else "w", encoding="utf-8") as err:
        result = subprocess.run(list(map(str, command)), cwd=directory, env=environment, stdout=out, stderr=err)
    require(result.returncode == 0, f"{Path(command[0]).name} exited {result.returncode}; see {directory / stderr}")


def execute_branch(directory, branch, records, payload, signature, state, tools, tool_info,
                   config, mode, replay, execution, sources):
    require(branch == "default", "Only the prepared entry orientation is supported")
    reverse_count = 0
    if state is not None and state.get("status") == "completed":
        print(f"Validated cached {branch} phylogeny: {directory}", flush=True)
        return
    directory.mkdir(parents=True, exist_ok=True)
    commands = portable_commands(config, replay)
    if state is None:
        (directory / "input.fasta").write_bytes(payload)
        save_json(directory / "commands.json", commands)
        state = dict(schema_version=1, branch=branch, mode=mode, phase="prepared", status="running",
                     threads=config["threads"], cache_signature=signature, input_sha256=payload_hash(payload),
                     source_sha256=sources, tools=tool_info, execution_attempts=[], created_at=now())
    state.update(status="running", execution=execution, updated_at=now())
    state["execution_attempts"].append(execution)
    save_json(directory / "run.json", state)
    environment = {**os.environ, **commands["environment"]}
    environment.pop("MAFFT_BINARIES", None)
    environment["PATH"] = os.pathsep.join([str(tools["mafft"].parent), str(tools["iqtree"].parent),
                                            environment.get("PATH", "")])
    try:
        if state["phase"] in {"prepared", "aligning"}:
            state["phase"] = "aligning"
            save_json(directory / "run.json", state)
            print(f"Aligning {branch} with {config['threads']} requested threads", flush=True)
            run_command([tools["mafft"], *commands["mafft"][1:]], directory, "aligned.fasta", "mafft.log", environment)
            metrics = alignment_metrics(directory, records, reverse_count)
            write_table(directory / "alignment_metrics.tsv", [dict(metric=key, value=value) for key, value in metrics.items()])
            state.update(phase="aligned", aligned_sha256=metrics["aligned_sha256"])
            save_json(directory / "run.json", state)
        if state["phase"] == "inferring" and completed_report_exists(directory):
            state["phase"] = "inferred"
        if state["phase"] in {"aligned", "inferring"}:
            resuming = state["phase"] == "inferring"
            if resuming:
                artifacts = list(directory.glob("tree.*"))
                require(not artifacts or (directory / "tree.ckp.gz").is_file(),
                        f"Incomplete IQ-TREE output lacks its checkpoint: {directory}")
            state["phase"] = "inferring"
            save_json(directory / "run.json", state)
            print(f"Inferring {branch} with IQ-TREE defaults and {config['threads']} threads", flush=True)
            run_command([tools["iqtree"], *commands["iqtree"][1:]], directory,
                        "iqtree.stdout.log", "iqtree.stderr.log", environment, append=resuming)
            state["phase"] = "inferred"
            save_json(directory / "run.json", state)
        metrics = alignment_metrics(directory, records, reverse_count)
        inference = inference_details(directory, [name for name, _ in records])
        threads = logged_threads(directory, config["threads"])
        comparison = expected_replay(replay, metrics, inference)
        require((directory / "commands.json").read_text() == json.dumps(commands, indent=2, sort_keys=True) + "\n",
                f"Recorded commands changed: {directory}")
        write_table(directory / "alignment_metrics.tsv", [dict(metric=key, value=value) for key, value in metrics.items()])
        outputs = ["input.fasta", "aligned.fasta", "alignment_metrics.tsv", "commands.json", "mafft.log",
                   "tree.treefile", "tree.iqtree", "tree.log", "iqtree.stdout.log", "iqtree.stderr.log"]
        state.update(**inference, **comparison, aligned_sha256=metrics["aligned_sha256"], logged_threads=threads,
                     status="completed", phase="completed", completed_at=now(), updated_at=now(),
                     output_sha256={name: sha256(directory / name) for name in outputs})
        state.pop("error", None)
        save_json(directory / "run.json", state)
        print(f"Completed {branch}: model={state['selected_model']}; seed={state['seed']}", flush=True)
        if comparison.get("exact_tree_match") is False:
            print(f"Warning: {branch} replay tree bytes differ; retain this result for topology comparison.",
                  file=sys.stderr, flush=True)
    except (Exception, KeyboardInterrupt) as error:
        state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                     error=str(error), updated_at=now())
        save_json(directory / "run.json", state)
        raise


@contextmanager
def workflow_lock(output):
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".workflow.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError(f"Another phylogeny process is using {output}") from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mafft", help="MAFFT 7.526 executable")
    parser.add_argument("--iqtree", help="IQ-TREE 3.1.4 executable")
    parser.add_argument("--check-tools", action="store_true", help="Check versions without computing")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    require(args.check_tools, "Use scripts/infer_entry.py --entry curated|auto for tree inference")
    _, tool_info = resolve_tools(configured(), args.mafft, args.iqtree)
    for name, info in tool_info.items():
        print(f"{name} {info['version']}: {info['executable']}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(f"Phylogeny failed: {error}")
