#!/usr/bin/env bash
#SBATCH --job-name=tail-p2-v4
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-tail-p2-v4-%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-tail-p2-v4-%j.err

set -euo pipefail

WORK_ROOT="/ibex/scratch/projects/c2014/nelly-local"
SIF_DIR="$WORK_ROOT/sif_cache"
export TMPDIR="$WORK_ROOT/tmp"
export SINGULARITY_TMPDIR="$WORK_ROOT/tmp"
export SINGULARITY_CACHEDIR="$WORK_ROOT/singularity_build_cache"
mkdir -p "$TMPDIR" "$SINGULARITY_CACHEDIR" "$WORK_ROOT/deliverables/TS" "$WORK_ROOT/deliverables/HD"

module load singularity 2>/dev/null || true

if [[ ! -s "$SIF_DIR/nelly-acmg-agent:v4.sif" ]]; then
  echo "[$(date -Is)] Building nelly-acmg-agent:v4.sif from $SIF_DIR/nelly-acmg-agent-v4.tar.gz ..."
  singularity build "$SIF_DIR/nelly-acmg-agent:v4.sif.tmp" "docker-archive://$SIF_DIR/nelly-acmg-agent-v4.tar.gz"
  mv "$SIF_DIR/nelly-acmg-agent:v4.sif.tmp" "$SIF_DIR/nelly-acmg-agent:v4.sif"
  ln -sf "nelly-acmg-agent:v4.sif" "$SIF_DIR/nelly-acmg-agent_v4.sif"
  echo "[$(date -Is)] Built $SIF_DIR/nelly-acmg-agent:v4.sif ($(stat -c %s "$SIF_DIR/nelly-acmg-agent:v4.sif") bytes)"
fi

cd "$WORK_ROOT/workflow"

echo "[$(date -Is)] === Running tail-phase2-TS-v4 (sequentially, 1/2) ==="
rm -rf "$WORK_ROOT/runs/tail-phase2-TS-v4"
bash local/run.sh "$WORK_ROOT/runs/tail-phase2-TS-v4" local/test-tail.cwl local/test-tail.phase2.TS.json

cp -f "$WORK_ROOT/runs/tail-phase2-TS-v4/acmg-report.pdf" "$WORK_ROOT/deliverables/TS/acmg-report.pdf"
cp -f "$WORK_ROOT/runs/tail-phase2-TS-v4/acmg-report.md" "$WORK_ROOT/deliverables/TS/acmg-report.md"
cp -f "$WORK_ROOT/runs/tail-phase2-TS-v4/acmg-assessment.json" "$WORK_ROOT/deliverables/TS/acmg-assessment.json"
echo "[$(date -Is)] === Staged TS deliverables to $WORK_ROOT/deliverables/TS/ ==="
ls -lh "$WORK_ROOT/deliverables/TS/"

echo "[$(date -Is)] === Running tail-phase2-HD-v4 (sequentially, 2/2) ==="
rm -rf "$WORK_ROOT/runs/tail-phase2-HD-v4"
bash local/run.sh "$WORK_ROOT/runs/tail-phase2-HD-v4" local/test-tail.cwl local/test-tail.phase2.HD.json

cp -f "$WORK_ROOT/runs/tail-phase2-HD-v4/acmg-report.pdf" "$WORK_ROOT/deliverables/HD/acmg-report.pdf"
cp -f "$WORK_ROOT/runs/tail-phase2-HD-v4/acmg-report.md" "$WORK_ROOT/deliverables/HD/acmg-report.md"
cp -f "$WORK_ROOT/runs/tail-phase2-HD-v4/acmg-assessment.json" "$WORK_ROOT/deliverables/HD/acmg-assessment.json"
echo "[$(date -Is)] === Staged HD deliverables to $WORK_ROOT/deliverables/HD/ ==="
ls -lh "$WORK_ROOT/deliverables/HD/"

echo "[$(date -Is)] === ALL DONE ==="
