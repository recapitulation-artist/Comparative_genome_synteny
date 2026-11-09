#!/usr/bin/env bash
# Run one sequential workflow within its requested Slurm allocation.
set -euo pipefail
threads="$1"
python="$2"
target="$3"
shift 3
: "${SLURM_JOB_ID:?Submit this script with scripts/submit_workflow.py}"
export OMP_NUM_THREADS="$threads"
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export SLURM_EXPORT_ENV=ALL
printf 'Job %s on %s; CPUs per task: %s\n' "$SLURM_JOB_ID" "$(hostname)" "${SLURM_CPUS_PER_TASK:-unset}"
# Validate inside srun, where affinity reflects the actual task allocation.
srun --nodes=1 --ntasks=1 --cpus-per-task="$threads" --cpu-bind=threads \
  "$python" -c 'import os, sys; sys.path.insert(0, "scripts"); from runtime import validate_threads; validate_threads(sys.argv[1]); os.execvp("make", ["make", "--no-print-directory", "-j1", *sys.argv[2:]])' \
  "$threads" "$target" "$@"
