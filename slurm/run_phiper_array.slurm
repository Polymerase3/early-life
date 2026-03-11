#!/bin/bash
#SBATCH --job-name=phiper
#SBATCH --array=1-57%10
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=45G
#SBATCH --output=/lisc/data/scratch/ccr/MFK/early_life_final/logs/%x-%A_%a.out
#SBATCH --error=/lisc/data/scratch/ccr/MFK/early_life_final/logs/%x-%A_%a.err
#SBATCH --chdir=/lisc/data/scratch/ccr/MFK/early_life_final
#SBATCH --no-requeue

set -euo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

apptainer exec \
  --bind /lisc/data/scratch/ccr/MFK/early_life_final:/work \
  /lisc/data/scratch/ccr/MFK/early_life_final/phiper_0.2.6.sif \
  Rscript /work/R/03-analysis.R \
    N_CORES=${SLURM_CPUS_PER_TASK} \
    MAX_GB=45
