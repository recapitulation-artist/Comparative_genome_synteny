#!/usr/bin/env python3
"""Align both region orientations and infer trees on a 96-CPU node2 Slurm task."""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import socket
import subprocess
import sys

from extract_regions import ROOT, fasta_bytes, fasta_records, require, table, write_table


CONFIG_PATH = ROOT / "config/phylogeny.json"
BRANCHES = ("default", "oriented")


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


def configured():
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    require(config["schema_version"] == 1 and config["threads"] == 96,
            "Expected the 96-thread phylogeny configuration")
    require(tuple(config["branches"]) == BRANCHES, "Expected default and oriented branches")
    require(config["mafft_options"] == ["--thread", "96", "--threadtb", "96"]
            and config["iqtree_options"] == ["-T", "96"],
            "Only the documented resource options are supported")
    return config


def resolve_tool(name, requested, version):
    candidates = [requested] if requested else [
        str(ROOT / ".tools/phylogeny/bin" / name),
        str(Path.home() / "miniforge3/envs/Megagenomic_profile/bin" / name),
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
    host = socket.gethostname().split(".", 1)[0]
    job = os.environ.get("SLURM_JOB_ID", "")
    task_cpus = os.environ.get("SLURM_CPUS_PER_TASK", "")
    node_cpus = os.environ.get("SLURM_CPUS_ON_NODE", "")
    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    require(host == config["required_node"], "Phylogeny computation is restricted to node2")
    require(job.isdigit(), "Phylogeny computation requires an active Slurm job")
    require(task_cpus.isdigit() and int(task_cpus) == config["threads"],
            "Slurm must allocate exactly 96 CPUs per task")
    require(node_cpus.isdigit() and int(node_cpus) >= config["threads"],
            "Slurm reports fewer than 96 CPUs on the allocated node")
    require(len(affinity) >= config["threads"], "Process CPU affinity exposes fewer than 96 CPUs")
    nodes = os.environ.get("SLURM_JOB_NUM_NODES", os.environ.get("SLURM_NNODES", ""))
    require(nodes == "1", "The phylogeny task requires a single-node Slurm allocation")
    return dict(node=host, job_id=job, cpus_per_task=int(task_cpus), cpus_on_node=int(node_cpus),
                cpu_affinity=affinity, slurm_job_nodelist=os.environ.get("SLURM_JOB_NODELIST", ""),
                slurm_step_id=os.environ.get("SLURM_STEP_ID", ""), platform=platform.platform(),
                python_version=sys.version.split()[0], started_at=now())


def prepare_inputs(config):
    sources = dict(sequences=ROOT / "data/sequences.csv", genomes=ROOT / "data/genomes.tsv",
                   regions=ROOT / "results/regions.tsv",
                   extracted_fasta=ROOT / "results/extracted_29_sequences.fasta")
    sequences = table(sources["sequences"], ",")
    regions, genomes = table(sources["regions"]), table(sources["genomes"])
    records = list(fasta_records(sources["extracted_fasta"]))
    order = [row["seq_id"] for row in sequences]
    require(len(order) == config["sequences"] and len(set(order)) == len(order),
            "Expected exactly 29 unique sequence IDs")
    require(len(records) == len(regions) == len(genomes) == len(order),
            "FASTA, region, genome and sequence counts differ")
    require([row["seq_id"] for row in regions] == order, "Region order differs from sequences.csv")
    genome_index = {row["seq_id"]: row for row in genomes}
    require(len(genome_index) == len(order) and set(genome_index) == set(order), "Genome IDs differ")
    normalized, oriented, species = [], [], []
    complement = str.maketrans("ACGTN", "TGCAN")
    for rank, ((header, dna), sequence, region) in enumerate(zip(records, sequences, regions), 1):
        identifier = sequence["seq_id"]
        genome = genome_index[identifier]
        dna = dna.upper()
        require(header == region["fasta_header"], f"FASTA/region header differs: {identifier}")
        require(len(dna) == int(sequence["length"]) == int(region["length"])
                == int(genome["region_length"]), f"Region length differs: {identifier}")
        require(dna and set(dna) <= set("ACGTN"), f"Unexpected DNA alphabet: {identifier}")
        require(region["assembly_accession"] == genome["assembly_accession"],
                f"Region/genome assembly differs: {identifier}")
        require(region["target_strand"] in {"+", "-"} and region["region_strand"] == "+",
                f"Unexpected source orientation: {identifier}")
        require(header == region["contig"] + "_" + region["target_protein_id"],
                f"FASTA header is not the mapped contig/target pair: {identifier}")
        reverse = region["target_strand"] == "-"
        normalized.append((identifier, dna))
        oriented.append((identifier, dna.translate(complement)[::-1] if reverse else dna))
        species.append(dict(rank=rank, seq_id=identifier, organism_name=genome["organism_name"],
                            strain=genome["strain"], assembly_accession=region["assembly_accession"],
                            assembly_status=genome["assembly_status"], genome_notes=genome.get("genome_notes", ""),
                            contig=region["contig"], target_strand=region["target_strand"],
                            target_protein_id=region["target_protein_id"], fasta_header=header,
                            length=len(dna), reverse_complemented_in_oriented_branch=reverse))
    require(sum(row["reverse_complemented_in_oriented_branch"] for row in species)
            == config["negative_target_strands"], "Expected 12 negative-strand targets")
    return dict(default=normalized, oriented=oriented), species, {key: sha256(path) for key, path in sources.items()}


def replay_settings(path, mode):
    if mode == "default":
        return {branch: None for branch in BRANCHES}
    manifest = json.loads(path.read_text(encoding="utf-8"))
    require(manifest.get("schema_version") == 1 and set(manifest["branches"]) == set(BRANCHES),
            "Replay manifest must define both branches with schema_version=1")
    branches = manifest["branches"]
    for branch, settings in branches.items():
        require(isinstance(settings["seed"], int) and not isinstance(settings["seed"], bool)
                and settings["seed"] > 0, f"Invalid replay seed: {branch}")
        for key in ("aligned_sha256", "tree_sha256"):
            require(re.fullmatch(r"[0-9a-f]{64}", settings[key]) is not None,
                    f"Invalid replay {key}: {branch}")
    return branches


def portable_commands(config, replay):
    inference = ["iqtree3", "-s", "aligned.fasta", *config["iqtree_options"], "--prefix", "tree"]
    if replay is not None:
        inference.extend(["--seed", str(replay["seed"])])
    return dict(mafft=["mafft", *config["mafft_options"], "input.fasta"],
                mafft_stdout="aligned.fasta", mafft_stderr="mafft.log", iqtree=inference,
                iqtree_stdout="iqtree.stdout.log", iqtree_stderr="iqtree.stderr.log",
                environment=dict(LC_ALL="C", OMP_NUM_THREADS=str(config["threads"]), OMP_DYNAMIC="FALSE"),
                unset_environment=["MAFFT_BINARIES"])


def alignment_metrics(directory, records, reversed_count):
    aligned = list(fasta_records(directory / "aligned.fasta"))
    require([name for name, _ in aligned] == [name for name, _ in records],
            "Alignment must preserve all 29 sequence IDs and their input order")
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
            "MAFFT log does not confirm the 96-thread request")
    require(expected in mafft.get("nthreadtb", set()), "MAFFT log does not confirm 96 progressive-stage threads")
    texts = []
    for name in ("tree.log", "iqtree.stdout.log"):
        path = directory / name
        if path.is_file():
            texts.append(path.read_text(encoding="utf-8", errors="replace"))
    inference = {int(value) for value in re.findall(r"^Kernel:.*?\b(\d+)\s+threads?\b", "\n".join(texts), re.MULTILINE)}
    require(inference == {expected}, f"IQ-TREE log does not confirm exactly 96 threads: {sorted(inference)}")
    return dict(mafft={key: sorted(values) for key, values in mafft.items()}, iqtree=sorted(inference))


def inference_details(directory, order):
    tree = (directory / "tree.treefile").read_text(encoding="utf-8").strip()
    require(tree.endswith(";"), "Incomplete IQ-TREE Newick output")
    tips = re.findall(r"(?:\(|,)\s*([A-Za-z0-9_.-]+)\s*(?=:|,|\))", tree)
    require(len(tips) == len(order) and len(set(tips)) == len(order) and set(tips) == set(order),
            "IQ-TREE tips differ from the 29 input sequence IDs")
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
    reverse_count = config["negative_target_strands"] if branch == "oriented" else 0
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
            print(f"Aligning {branch} with 96 requested threads", flush=True)
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
            print(f"Inferring {branch} with IQ-TREE defaults and 96 threads", flush=True)
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


def workflow(output, inputs, species, sources, replays, mode, config, tools, tool_info):
    plans = []
    for branch in BRANCHES:
        records, replay = inputs[branch], replays[branch]
        payload = fasta_bytes(records)
        if replay and replay.get("input_sha256"):
            require(payload_hash(payload) == replay["input_sha256"], f"Replay input differs: {branch}")
        signature = payload_hash(json.dumps(dict(schema_version=1, branch=branch, mode=mode,
                                     input_sha256=payload_hash(payload), source_sha256=sources,
                                     config_sha256=sha256(CONFIG_PATH), script_sha256=sha256(Path(__file__)),
                                     runtime_lock_sha256=sha256(ROOT / config["runtime_lock"]),
                                     tool_sha256={name: info["executable_sha256"] for name, info in tool_info.items()},
                                     replay=replay), sort_keys=True).encode("utf-8"))
        directory = output / branch
        state = cached_state(directory, signature, payload, records, replay,
                             config["negative_target_strands"] if branch == "oriented" else 0, config)
        plans.append((directory, branch, records, payload, signature, state, replay))
    needs_compute = any(state is None or state.get("status") != "completed"
                        for _, _, _, _, _, state, _ in plans)
    execution = scheduler_context(config) if needs_compute else None
    output.mkdir(parents=True, exist_ok=True)
    species_path = output / "species.tsv"
    if species_path.exists():
        require(table(species_path) == [{key: str(value).lower() if isinstance(value, bool) else str(value)
                                        for key, value in row.items()} for row in species],
                "Existing species.tsv differs from current inputs")
    else:
        write_table(species_path, species)
    for directory, branch, records, payload, signature, state, replay in plans:
        execute_branch(directory, branch, records, payload, signature, state, tools, tool_info,
                       config, mode, replay, execution, sources)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("default", "replay"), default="default")
    parser.add_argument("--mafft", help="MAFFT 7.526 executable")
    parser.add_argument("--iqtree", help="IQ-TREE 3.1.4 executable")
    parser.add_argument("--check-tools", action="store_true", help="Check versions without scheduler requirements")
    parser.add_argument("--output-root", type=Path, default=ROOT / "results/phylogeny")
    parser.add_argument("--replay-manifest", type=Path, default=ROOT / "config/phylogeny_replay.json")
    args = parser.parse_args()
    require(sys.version_info >= (3, 9), "Python 3.9 or newer is required")
    config = configured()
    tools, tool_info = resolve_tools(config, args.mafft, args.iqtree)
    if args.check_tools:
        for name, info in tool_info.items():
            print(f"{name} {info['version']}: {info['executable']}")
        return
    inputs, species, sources = prepare_inputs(config)
    manifest = args.replay_manifest if args.replay_manifest.is_absolute() else ROOT / args.replay_manifest
    replays = replay_settings(manifest, args.mode)
    output = args.output_root if args.output_root.is_absolute() else ROOT / args.output_root
    output = output.resolve()
    require(output != ROOT and ROOT / "data" not in (output, *output.parents),
            "Output must not overwrite repository inputs")
    with workflow_lock(output):
        workflow(output, inputs, species, sources, replays, args.mode, config, tools, tool_info)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(f"Phylogeny failed: {error}")
