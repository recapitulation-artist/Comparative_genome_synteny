#!/usr/bin/env python3
"""Shared CPU budget and execution metadata for local and Slurm workflows."""

from datetime import datetime, timezone
import os
import platform
import socket
import sys


DEFAULT_THREADS = 56


def available_cpus():
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def validate_threads(value):
    if isinstance(value, bool) or not str(value).isdigit() or int(value) < 1:
        raise ValueError("THREADS must be a positive integer")
    threads = int(value)
    available = available_cpus()
    if threads > available:
        raise ValueError(f"Requested THREADS={threads}, but CPU affinity exposes {available} CPUs; set THREADS accordingly")
    if os.environ.get("SLURM_JOB_ID"):
        allocated = os.environ.get("SLURM_CPUS_PER_TASK", "")
        if not allocated.isdigit() or int(allocated) < threads:
            raise ValueError(f"Slurm must allocate at least {threads} CPUs per task; request --cpus-per-task={threads} or reduce THREADS")
        if os.environ.get("SLURM_NTASKS", "1") != "1":
            raise ValueError("The workflow requires one Slurm task with multiple CPUs")
    return threads


def execution_context(threads):
    threads = validate_threads(threads)
    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    return dict(mode="slurm" if os.environ.get("SLURM_JOB_ID") else "local",
                threads=threads, node=socket.gethostname(),
                job_id=os.environ.get("SLURM_JOB_ID", ""),
                cpus_per_task=int(os.environ["SLURM_CPUS_PER_TASK"]) if os.environ.get("SLURM_JOB_ID") else None,
                cpus_on_node=os.environ.get("SLURM_CPUS_ON_NODE", ""),
                available_cpus=available_cpus(), cpu_affinity=affinity,
                slurm_job_nodelist=os.environ.get("SLURM_JOB_NODELIST", ""),
                slurm_step_id=os.environ.get("SLURM_STEP_ID", ""),
                platform=platform.platform(), python_version=sys.version.split()[0],
                started_at=datetime.now(timezone.utc).isoformat(timespec="seconds"))


def thread_environment(threads=1):
    """Limit nested libraries while the native tool manages its CPU budget."""
    return dict(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1",
                MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
                NUMEXPR_NUM_THREADS="1")
