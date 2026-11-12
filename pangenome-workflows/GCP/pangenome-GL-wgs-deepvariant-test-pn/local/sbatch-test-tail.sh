#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-tail-%x-%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-tail-%x-%j.err

# Slurm batch script to run the full local/test-tail.cwl
# (phenofrommetadata + exclude-eh-regions + slivar -> vep -> exomiser
#  -> acmg-extract-candidates -> vep-candidates -> acmg-agent with Unimatrix Qwen)
# on IBEX.
set -euo pipefail

COHORT="${1:-${COHORT:-TS}}"
WORK_ROOT="${NELLY_LOCAL_ROOT:-/ibex/scratch/projects/c2014/nelly-local}"
WF_DIR="$WORK_ROOT/workflow"
OUTDIR="${2:-$WORK_ROOT/runs/tail-${COHORT}-full}"
INPUT_JSON="${3:-local/test-tail.${COHORT}.json}"

export HOST="ibex"
export NELLY_LOCAL_ROOT="$WORK_ROOT"

cd "$WF_DIR"
bash local/run.sh \
  "$OUTDIR" \
  local/test-tail.cwl \
  "$INPUT_JSON"
