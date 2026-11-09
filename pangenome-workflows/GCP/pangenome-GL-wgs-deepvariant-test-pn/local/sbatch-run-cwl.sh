#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=24
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-cwl-%x-%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-cwl-%x-%j.err

# Generic Slurm wrapper to run any CWL workflow via local/run.sh on IBEX:
#   sbatch --job-name=<name> local/sbatch-run-cwl.sh <outdir> <workflow.cwl> <inputs.json>
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <outdir> <workflow.cwl> <inputs.json> [extra cwltool args...]" >&2
  exit 2
fi

OUTDIR="$1"
WORKFLOW="$2"
INPUTS_JSON="$3"
shift 3

WORK_ROOT="${NELLY_LOCAL_ROOT:-/ibex/scratch/projects/c2014/nelly-local}"
WF_DIR="$WORK_ROOT/workflow"
export HOST="ibex"
export NELLY_LOCAL_ROOT="$WORK_ROOT"
export NELLY_STEP_CORES="${SLURM_CPUS_PER_TASK:-24}"
export TMPDIR="$WORK_ROOT/tmp/job-${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR" "$OUTDIR" "$WORK_ROOT/logs"

# Ensure nelly-acmg-agent:v3.sif and vep_cache are staged (idempotent, lock-guarded)
STAGE_ONLY=1 bash "$WF_DIR/local/sbatch-stage-and-run-tail-ibex.sh"

cd "$WF_DIR"
bash local/run.sh "$OUTDIR" "$WORKFLOW" "$INPUTS_JSON" "$@"
