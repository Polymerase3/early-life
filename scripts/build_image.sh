#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE_NAME="phiper-analysis:latest"

export TMPDIR=$(mktemp -d)

docker build -t "$IMAGE_NAME" "$PROJECT_DIR"

rm -rf "$TMPDIR"
