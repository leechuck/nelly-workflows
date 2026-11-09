#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-align-%x-%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-align-%x-%j.err

# Slurm batch script to stage inputs/SIFs (on compute node, never login node)
# and run align-and-call.cwl via local/run.sh on IBEX.
#
# Submit examples:
#   sbatch --job-name=align-TS-standard local/sbatch-align-and-call.sh TS standard
#   sbatch --job-name=align-HD-standard local/sbatch-align-and-call.sh HD standard
set -euo pipefail

COHORT="${1:-${COHORT:-TS}}"
DV_MODE="${2:-${DV_MODE:-standard}}"
COHORT_LC="$(echo "$COHORT" | tr '[:upper:]' '[:lower:]')"

WORK_ROOT="${NELLY_LOCAL_ROOT:-/ibex/scratch/projects/c2014/nelly-local}"
WF_DIR="$WORK_ROOT/workflow"
INPUTS="$WORK_ROOT/inputs"
SIF_CACHE="$WORK_ROOT/sif_cache"
export TMPDIR="$WORK_ROOT/tmp/job-${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR" "$INPUTS/graph" "$INPUTS/patients/$COHORT_LC" "$SIF_CACHE/.cache" "$WORK_ROOT/runs" "$WORK_ROOT/logs"

if [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module load singularity

export SINGULARITY_CACHEDIR="$SIF_CACHE/.cache"
export CWL_SINGULARITY_CACHE="$SIF_CACHE"
export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_BINDPATH="/ibex/scratch/projects/c2014,/encrypted/e3008,$TMPDIR"

ARV_GET="$WORK_ROOT/venv/bin/arv-get"
PYTHON="$WORK_ROOT/venv/bin/python"

echo "[$(date -Is)] Node=$(hostname) Job=${SLURM_JOB_ID:-local} Cohort=$COHORT DV_MODE=$DV_MODE"

# Helper: atomic lock-guarded SIF staging
stage_sif_copy() {
  local sif_name="$1"
  local src_path="$2"
  local dst="$SIF_CACHE/$sif_name"
  local lock="$SIF_CACHE/.lock-${sif_name}"
  if [ -s "$dst" ]; then return 0; fi
  if mkdir "$lock" 2>/dev/null; then
    echo "[$(date -Is)] Copying pre-built SIF $src_path -> $dst"
    cp -f "$src_path" "$dst.tmp.$$"
    mv -f "$dst.tmp.$$" "$dst"
    rmdir "$lock" 2>/dev/null || true
  else
    echo "[$(date -Is)] Waiting for another job to finish staging $dst..."
    while [ ! -s "$dst" ]; do sleep 5; done
  fi
}

stage_sif_pull() {
  local docker_uri="$1"
  local sif_name
  sif_name="$("$PYTHON" -c "from cwl_utils.image_puller import SingularityImagePuller; print(SingularityImagePuller('$docker_uri', '$SIF_CACHE', 'singularity', False).get_image_name())")"
  local dst="$SIF_CACHE/$sif_name"
  local lock="$SIF_CACHE/.lock-${sif_name}"
  if [ -s "$dst" ]; then return 0; fi
  if mkdir "$lock" 2>/dev/null; then
    echo "[$(date -Is)] Pulling docker://$docker_uri -> $dst"
    singularity pull --force "$dst.tmp.$$" "docker://$docker_uri"
    mv -f "$dst.tmp.$$" "$dst"
    rmdir "$lock" 2>/dev/null || true
  else
    echo "[$(date -Is)] Waiting for another job to finish pulling $dst..."
    while [ ! -s "$dst" ]; do sleep 5; done
  fi
}

# 1. Stage SIF images needed for align-and-call.cwl
EXISTING_TOOLS="/encrypted/e3008/Rob/yang-pangenome-full-cohort-20261012/modernisation-20261108/tools"
if [ -s "$EXISTING_TOOLS/dv110.sif" ]; then
  stage_sif_copy "google_s_deepvariant:1.10.0.sif" "$EXISTING_TOOLS/dv110.sif"
else
  stage_sif_pull "google/deepvariant:1.10.0"
fi
if [ -s "$EXISTING_TOOLS/padv110.sif" ]; then
  stage_sif_copy "google_s_deepvariant:pangenome___aware___deepvariant-1.10.0.sif" "$EXISTING_TOOLS/padv110.sif"
fi

stage_sif_pull "quay.io/vgteam/vg:v1.78.0"
stage_sif_pull "coolmaksat/kmc"
stage_sif_pull "coolmaksat/minimap2-samtools"
stage_sif_pull "dunglam2k/expansionhunter-grch38:v5.0.0"

# 2. Stage patient FASTQ + metadata for this job's cohort first (so TS and HD download in parallel!)
if [ "$COHORT" = "TS" ]; then
  PATIENT_PDH="5c88241c12f7dc7553ffb3f2ac1ba54b+39418"
elif [ "$COHORT" = "HD" ]; then
  PATIENT_PDH="d03401bdd2443e6340903acb552fd0c3+39418"
else
  echo "Unknown COHORT=$COHORT" >&2
  exit 2
fi

for f in metadata.yaml reads1.fastq.gz reads2.fastq.gz; do
  dst="$INPUTS/patients/$COHORT_LC/$f"
  if [ ! -s "$dst" ]; then
    echo "[$(date -Is)] Downloading $COHORT $f from Keep ($PATIENT_PDH/$f) -> $dst"
    "$ARV_GET" -f --threads 16 "$PATIENT_PDH/$f" "$dst.partial.$$"
    mv -f "$dst.partial.$$" "$dst"
  fi
done

# 3. Stage shared pangenome graph + reference (atomic lock per file so both jobs can help or wait cleanly)
stage_keep_file() {
  local pdh_path="$1"
  local dst="$2"
  local base
  base="$(basename "$dst")"
  local lock="$INPUTS/graph/.lock-${base}"
  if [ -s "$dst" ]; then return 0; fi
  if mkdir "$lock" 2>/dev/null; then
    echo "[$(date -Is)] Downloading $pdh_path -> $dst"
    "$ARV_GET" -f --threads 16 "$pdh_path" "$dst.partial.$$"
    mv -f "$dst.partial.$$" "$dst"
    rmdir "$lock" 2>/dev/null || true
  else
    echo "[$(date -Is)] Waiting for another job to finish downloading $dst..."
    while [ ! -s "$dst" ]; do sleep 5; done
  fi
}

GRAPH_PDH="1d5ad2aefb6e809dddaffea3a0558415+5657"
VG170_IDX_PDH="d704d713f7039bbc4026d55fafec7c33+2743"

stage_keep_file "$GRAPH_PDH/ksa_hg38.gbz" "$INPUTS/graph/ksa_hg38.gbz"
stage_keep_file "$VG170_IDX_PDH/ksa_hg38.hapl" "$INPUTS/graph/ksa_hg38.hapl"
stage_keep_file "$GRAPH_PDH/hg38.fa" "$INPUTS/graph/hg38.fa"
stage_keep_file "$GRAPH_PDH/hg38.fa.fai" "$INPUTS/graph/hg38.fa.fai"

# 4. Run align-and-call.cwl via local/run.sh
cd "$WF_DIR"
OUTDIR="$WORK_ROOT/runs/align-${COHORT}-${DV_MODE}"
export HOST="ibex"
export NELLY_LOCAL_ROOT="$WORK_ROOT"

CORES="${SLURM_CPUS_PER_TASK:-16}"
HALF_CORES=$(( CORES > 16 ? 16 : CORES ))
OVERRIDES_FILE="$TMPDIR/cwltool-overrides.yml"
cat > "$OVERRIDES_FILE" <<EOF
cwltool:overrides:
  vg-giraffe-bam.cwl:
    requirements:
      ResourceRequirement:
        coresMin: ${CORES}
  kmc.cwl:
    requirements:
      ResourceRequirement:
        coresMin: ${HALF_CORES}
  samtools-sort.cwl:
    requirements:
      ResourceRequirement:
        coresMin: ${HALF_CORES}
  samtools-markdup.cwl:
    requirements:
      ResourceRequirement:
        coresMin: ${HALF_CORES}
  deepvariant.cwl:
    requirements:
      ResourceRequirement:
        coresMin: ${CORES}
EOF

bash local/run.sh \
  "$OUTDIR" \
  align-and-call.cwl \
  "local/align-and-call.${COHORT}.${DV_MODE}.json" \
  --overrides "$OVERRIDES_FILE"
