#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --time=12:00:00

# Slurm batch wrapper for Phase 3 GIAB benchmark runs (and main-vg runs) on IBEX:
# 1. Atomically stages any required GIAB / legacy inputs from Arvados Keep into
#    /ibex/scratch/projects/c2014/nelly-local/inputs/ with lock-dir coordination
#    and exact byte-size verification so concurrent jobs share staged inputs safely.
# 2. Pre-pulls legacy Singularity images (vg 1.54.0, DeepVariant 1.6.0) when running
#    benchmark-giab-legacy.cwl to verify clean .sif conversion up-front.
# 3. Invokes local/run.sh <outdir> <workflow.cwl> <inputs.json>.
#
# Usage:
#   sbatch --job-name=<name> --output=<outdir>/slurm-%j.out \
#     local/sbatch-benchmark-giab.sh <outdir> <workflow.cwl> <inputs.json>
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
REPO_DIR="$WORK_ROOT/workflow"
export PATH="$WORK_ROOT/venv/bin:$PATH"
export TMPDIR="$WORK_ROOT/tmp"
mkdir -p "$TMPDIR" "$OUTDIR"

if [ -f "$HOME/.config/arvados/settings.conf" ]; then
  while IFS='=' read -r k v; do
    k="$(echo "$k" | tr -d '[:space:]')"
    v="$(echo "$v" | tr -d '[:space:]')"
    if [ -n "$k" ] && [[ ! "$k" =~ ^# ]]; then
      export "$k=$v"
    fi
  done < "$HOME/.config/arvados/settings.conf"
fi

if ! command -v singularity >/dev/null 2>&1; then
  if [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  module load singularity
fi
mkdir -p "$WORK_ROOT/sif_cache/.cache"
export SINGULARITY_CACHEDIR="$WORK_ROOT/sif_cache/.cache"
export CWL_SINGULARITY_CACHE="$WORK_ROOT/sif_cache"
export SINGULARITY_TMPDIR="$TMPDIR"

stage_keep_file() {
  local keep_src="$1"
  local dest="$2"
  local expected_bytes="$3"
  local dest_dir
  dest_dir="$(dirname "$dest")"
  mkdir -p "$dest_dir"

  if [ -f "$dest" ]; then
    local actual_bytes
    actual_bytes="$(stat -c%s "$dest" 2>/dev/null || echo 0)"
    if [ "$actual_bytes" = "$expected_bytes" ]; then
      echo "[$(date -Is)] Verified existing $dest ($actual_bytes bytes)"
      return 0
    fi
  fi

  local lock_dir="${dest}.lock"
  while ! mkdir "$lock_dir" 2>/dev/null; do
    echo "[$(date -Is)] Waiting for concurrent job staging $dest ..."
    sleep 10
    if [ -f "$dest" ]; then
      local actual_bytes
      actual_bytes="$(stat -c%s "$dest" 2>/dev/null || echo 0)"
      if [ "$actual_bytes" = "$expected_bytes" ]; then
        echo "[$(date -Is)] Concurrent job finished staging $dest ($actual_bytes bytes)"
        return 0
      fi
    fi
  done

  # Re-check after acquiring lock
  if [ -f "$dest" ]; then
    local actual_bytes
    actual_bytes="$(stat -c%s "$dest" 2>/dev/null || echo 0)"
    if [ "$actual_bytes" = "$expected_bytes" ]; then
      rmdir "$lock_dir" 2>/dev/null || true
      echo "[$(date -Is)] Verified existing $dest ($actual_bytes bytes)"
      return 0
    fi
  fi

  local tmp_file="${dest}.part.$$"
  rm -f "$tmp_file"
  echo "[$(date -Is)] Staging $keep_src -> $dest ($expected_bytes bytes) ..."
  if ! "$WORK_ROOT/venv/bin/arv-get" -f --threads 16 "$keep_src" "$tmp_file"; then
    rm -f "$tmp_file"
    rmdir "$lock_dir" 2>/dev/null || true
    echo "[$(date -Is)] ERROR: arv-get failed for $keep_src" >&2
    return 1
  fi

  local dl_bytes
  dl_bytes="$(stat -c%s "$tmp_file" 2>/dev/null || echo 0)"
  if [ "$dl_bytes" != "$expected_bytes" ]; then
    rm -f "$tmp_file"
    rmdir "$lock_dir" 2>/dev/null || true
    echo "[$(date -Is)] ERROR: size mismatch for $dest (got $dl_bytes, expected $expected_bytes)" >&2
    return 1
  fi

  mv -f "$tmp_file" "$dest"
  rmdir "$lock_dir" 2>/dev/null || true
  echo "[$(date -Is)] Successfully staged $dest ($dl_bytes bytes)"
}

cd "$REPO_DIR"

# Stage GIAB inputs if referenced by INPUTS_JSON
if grep -q "giab/hg003" "$INPUTS_JSON"; then
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "$WORK_ROOT/inputs/giab/truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" 147212428
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
    "$WORK_ROOT/inputs/giab/truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" 1651307
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "$WORK_ROOT/inputs/giab/truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" 13219856
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/GRCh38_alldifficultregions.bed.gz" \
    "$WORK_ROOT/inputs/giab/truth/GRCh38_alldifficultregions.bed.gz" 31493051
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/GRCh38_notinalldifficultregions.bed.gz" \
    "$WORK_ROOT/inputs/giab/truth/GRCh38_notinalldifficultregions.bed.gz" 32740449
  stage_keep_file "4bc83e4703035a3eeb142eb17fcc1ce2+44786/HG003_R1.fastq.gz" \
    "$WORK_ROOT/inputs/giab/hg003/HG003_R1.fastq.gz" 31969557485
  stage_keep_file "4bc83e4703035a3eeb142eb17fcc1ce2+44786/HG003_R2.fastq.gz" \
    "$WORK_ROOT/inputs/giab/hg003/HG003_R2.fastq.gz" 33451687186
fi

if grep -q "giab/hg002" "$INPUTS_JSON"; then
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "$WORK_ROOT/inputs/giab/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" 156252944
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
    "$WORK_ROOT/inputs/giab/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" 1657747
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "$WORK_ROOT/inputs/giab/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" 11494021
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG002.novaseq.pcr-free.30x.R1.fastq.gz" \
    "$WORK_ROOT/inputs/giab/hg002/HG002.novaseq.pcr-free.30x.R1.fastq.gz" 25664092287
  stage_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970/HG002.novaseq.pcr-free.30x.R2.fastq.gz" \
    "$WORK_ROOT/inputs/giab/hg002/HG002.novaseq.pcr-free.30x.R2.fastq.gz" 26634286830
fi

if grep -q "ksa_hg38.dist" "$INPUTS_JSON"; then
  stage_keep_file "1d5ad2aefb6e809dddaffea3a0558415+5657/ksa_hg38.dist" \
    "$WORK_ROOT/inputs/graph/ksa_hg38.dist" 1567026720
  stage_keep_file "56e94ed8e0d67d0c8e2616959d23ef2e+23423/ksa_hg38.min" \
    "$WORK_ROOT/inputs/graph/ksa_hg38.min" 37348221952
fi

if [[ "$WORKFLOW" == *"legacy"* ]]; then
  echo "[$(date -Is)] Pre-pulling and verifying legacy Singularity images ..."
  "$WORK_ROOT/venv/bin/python" - <<'PY'
import os, sys
from cwltool.singularity import SingularityCommandLineJob
for img in ["quay.io/biocontainers/vg:1.54.0--h9ee0642_0", "google/deepvariant:1.6.0"]:
    req = {"dockerPull": img}
    ok = SingularityCommandLineJob.get_image(req, True, os.environ["TMPDIR"] + "/out-")
    print(f"[legacy-sif] {img} -> ok={ok} path={req.get('dockerImageId')}", flush=True)
    if not ok:
        sys.exit(1)
PY
fi

exec bash "$REPO_DIR/local/run.sh" "$OUTDIR" "$WORKFLOW" "$INPUTS_JSON" "$@"
