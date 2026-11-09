#!/usr/bin/env bash
# Execute the selected workflow within a single node2 allocation.
set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:?Submit through make submit from the repository}"
export OMP_NUM_THREADS=96
export OPENBLAS_NUM_THREADS=96
export MKL_NUM_THREADS=96
export PYTHONUNBUFFERED=1
python3 scripts/submit_workflow.py --check-node
printf 'Job %s on %s; CPUs per task: %s\n' "$SLURM_JOB_ID" "$(hostname)" "$SLURM_CPUS_PER_TASK"
srun --nodes=1 --ntasks=1 --cpus-per-task=96 --cpu-bind=threads \
  make --no-print-directory -j1 "$1" MODE="$2" RSCRIPT="$3" BEDTOOLS="$4" \
  BLASTN="$5" MAKEBLASTDB="$6" MAFFT="$7" IQTREE="$8" RAW_PHYLOGENY_DIR="$9" THREADS=96
