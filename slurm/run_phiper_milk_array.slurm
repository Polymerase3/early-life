#!/bin/bash
#SBATCH --job-name=phiper_mmilk
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=45G
#SBATCH --output=/lisc/data/scratch/ccr/MFK/early_life_milk/logs/%x-%A_%a.out
#SBATCH --error=/lisc/data/scratch/ccr/MFK/early_life_milk/logs/%x-%A_%a.err
#SBATCH --chdir=/lisc/data/scratch/ccr/MFK/early_life_milk
#SBATCH --no-requeue

set -euo pipefail
KEY_FILE="scripts/mom_milk_cmp_keys.txt"
TASK_ID="${SLURM_ARRAY_TASK_ID:?missing array id}"
CMP_KEY="$(sed -n "${TASK_ID}p" "${KEY_FILE}")"
[ -n "${CMP_KEY}" ] || { echo "No CMP_KEY for task ${TASK_ID}"; exit 1; }

# R/03-analysis.R auto-reads SLURM_ARRAY_TASK_ID as CMP_INDEX.
# We run by CMP_KEY here, so remove it to avoid a key/index conflict.
unset SLURM_ARRAY_TASK_ID

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

apptainer exec \
  --bind /lisc/data/scratch/ccr/MFK/early_life_milk:/work \
  /lisc/data/scratch/ccr/MFK/early_life_milk/phiper_0.2.6.sif \
  Rscript /work/R/03-analysis.R \
    N_CORES=${SLURM_CPUS_PER_TASK} \
    MAX_GB=45 \
    "CMP_KEY=${CMP_KEY}"
