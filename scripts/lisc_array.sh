#!/bin/bash
#
#SBATCH --job-name=phiper_cmp
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=100G
#SBATCH --array=1-39
#SBATCH --output=%x-%A_%a-%N.out
#SBATCH --error=%x-%A_%a-%N.err
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE_NAME="phiper-analysis:latest"

export TMPDIR=$(mktemp -d)
export XDG_RUNTIME_DIR="$TMPDIR"

WORK_DIR="$TMPDIR/project"
mkdir -p "$WORK_DIR"

rsync -a --exclude '.git' --exclude 'results' "$PROJECT_DIR/" "$WORK_DIR/"

pushd "$WORK_DIR" > /dev/null

docker run --rm \
  -v "$WORK_DIR:/work" \
  -w /work \
  --cpus 30 \
  --memory 100g \
  "$IMAGE_NAME" \
  R/02-analysis.R N_CORES=30 MAX_GB=100 LOG=TRUE CMP_INDEX="$SLURM_ARRAY_TASK_ID"

mkdir -p "$PROJECT_DIR/results"
rsync -a "$WORK_DIR/results/" "$PROJECT_DIR/results/"

popd > /dev/null
rm -rf "$TMPDIR"
