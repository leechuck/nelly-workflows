#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --job-name=tail-acmg-ibex
#SBATCH --output=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-tail-acmg-%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-tail-acmg-%j.err

# Stage all tail SIFs + inputs in parallel with any active staging job, then run
# the full local/test-tail.cwl (including acmg-extract-candidates -> vep-candidates -> acmg-agent
# with Unimatrix Qwen) on TS and HD.
set -euo pipefail

WORK_ROOT="${NELLY_LOCAL_ROOT:-/ibex/scratch/projects/c2014/nelly-local}"
WF_DIR="$WORK_ROOT/workflow"
INPUTS="$WORK_ROOT/inputs"
SIF_CACHE="$WORK_ROOT/sif_cache"
IMAGES="$WORK_ROOT/images"
export TMPDIR="$WORK_ROOT/tmp/tail-acmg-${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR" "$INPUTS"/{graph,slivar_vep_custom,vep_fasta,vep_cache,exomiser,acmg,tail/ts,tail/hd,patients/ts,patients/hd} "$SIF_CACHE/.cache" "$IMAGES"

if [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module load singularity

export SINGULARITY_CACHEDIR="$SIF_CACHE/.cache"
export CWL_SINGULARITY_CACHE="$SIF_CACHE"
export SINGULARITY_TMPDIR="$TMPDIR"

ARV_GET="$WORK_ROOT/venv/bin/arv-get"
PYTHON="$WORK_ROOT/venv/bin/python"

stage_sif_pull_wait() {
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
  fi
}

stage_sif_from_keep_tar() {
  local sif_name="$1"
  local keep_tar_uri="$2"
  local dst="$SIF_CACHE/$sif_name"
  local lock="$SIF_CACHE/.lock-${sif_name}"
  if [ -s "$dst" ]; then return 0; fi
  if mkdir "$lock" 2>/dev/null; then
    local safe_base
    safe_base="$(echo "${sif_name%.sif}" | tr ':' '-')"
    local tar_path="$IMAGES/${safe_base}.tar"
    if [ ! -s "$tar_path" ]; then
      echo "[$(date -Is)] Downloading Docker tar $keep_tar_uri -> $tar_path"
      "$ARV_GET" -f --threads 16 "$keep_tar_uri" "$tar_path.partial.$$"
      mv -f "$tar_path.partial.$$" "$tar_path"
    fi
    echo "[$(date -Is)] Building $dst from docker-archive://$tar_path"
    singularity build --force "$dst.tmp.$$" "docker-archive://$tar_path"
    mv -f "$dst.tmp.$$" "$dst"
    rmdir "$lock" 2>/dev/null || true
  fi
}

stage_keep_file() {
  local pdh_path="$1"
  local dst="$2"
  local dir
  dir="$(dirname "$dst")"
  local base
  base="$(basename "$dst")"
  mkdir -p "$dir"
  local lock="$dir/.lock-${base}"
  if [ -s "$dst" ]; then return 0; fi
  if mkdir "$lock" 2>/dev/null; then
    echo "[$(date -Is)] Downloading $pdh_path -> $dst"
    "$ARV_GET" -f --threads 16 "$pdh_path" "$dst.partial.$$"
    mv -f "$dst.partial.$$" "$dst"
    rmdir "$lock" 2>/dev/null || true
  fi
}

# 1. Stage tail SIFs first
stage_sif_from_keep_tar "nelly-exomiser:v4.sif" "7712a129b6c89fb29f6f39eb5842ad90+429/sha256:1e800dce4081b13665db09bcb3e650f8fe38d241ab2eea2f2031a205ecda6dcd.tar"
stage_sif_from_keep_tar "nelly-acmg-agent:v3.sif" "49f9dd239b70c0fa5a67dc51edea9e99+1731/sha256:dda1610e81cdbac256e71f1eb78d1bc6b5308da8f332c440e897c8261d7a1877.tar"
ln -sf "nelly-acmg-agent:v3.sif" "$SIF_CACHE/nelly-acmg-agent:v2.sif"
ln -sf "nelly-acmg-agent:v3.sif" "$SIF_CACHE/nelly-acmg-agent:v1.sif"

for uri in \
  "ensemblorg/ensembl-vep:release_116.3" \
  "brentp/slivar:v0.3.0" \
  "dunglam2k/phenofrommetadata:v1.01" \
  "staphb/bcftools:1.22" \
  "dunglam2k/expansionhunter-grch38:v5.0.0" \
  "coolmaksat/minimap2-samtools" \
  "coolmaksat/kmc"; do
  stage_sif_pull_wait "$uri"
done

# 2. Stage tail inputs first
GRAPH_PDH="1d5ad2aefb6e809dddaffea3a0558415+5657"
stage_keep_file "$GRAPH_PDH/hg38.fa.fai" "$INPUTS/graph/hg38.fa.fai"

for f in metadata.yaml; do
  stage_keep_file "5c88241c12f7dc7553ffb3f2ac1ba54b+39418/$f" "$INPUTS/patients/ts/$f"
  stage_keep_file "d03401bdd2443e6340903acb552fd0c3+39418/$f" "$INPUTS/patients/hd/$f"
done

SLIVAR_PDH="148bcb74e5d8268ac486d315e657f70b+3240"
for f in \
  gnomad.hg38.genomes.v3.fix.zip \
  sorted.pavs.plugin.38.vcf.gz \
  sorted.pavs.plugin.38.vcf.gz.tbi \
  sorted.plugin.go.bed.gz \
  sorted.plugin.go.bed.gz.tbi \
  sorted.plugin.pheno.bed.gz \
  sorted.plugin.pheno.bed.gz.tbi \
  sorted.plugin.ppi.bed.gz \
  sorted.plugin.ppi.bed.gz.tbi; do
  stage_keep_file "$SLIVAR_PDH/$f" "$INPUTS/slivar_vep_custom/$f"
done

VEP_FASTA_PDH="00528ccfa8b07fd08907077b1f4fe3e9+927"
for f in fastafile fastafile.fai fastafile.gzi; do
  stage_keep_file "$VEP_FASTA_PDH/$f" "$INPUTS/vep_fasta/$f"
done
touch "$INPUTS/vep_fasta/fastafile.fai" "$INPUTS/vep_fasta/fastafile.gzi"

for f in normalized.vcf.gz normalized.vcf.gz.tbi eh_repeats.vcf.gz eh_repeats.vcf.gz.tbi eh_repeats.json; do
  stage_keep_file "b2c8c2b940de5bf4d7881b652e12c24f+18026/$f" "$INPUTS/tail/ts/$f"
  stage_keep_file "6d04e84a0ab500c5036ec6713756783b+17943/$f" "$INPUTS/tail/hd/$f"
done

ACMG_PDH="28a24f20fb9659e2ba4e4b9a17f7bec0+778"
for f in clinvar.vcf.gz clinvar.vcf.gz.tbi gnomad.v4.1.constraint_metrics.tsv variant_summary.txt.gz; do
  stage_keep_file "$ACMG_PDH/$f" "$INPUTS/acmg/$f"
done
for f in gene_table.tsv STRchive-loci.json; do
  stage_keep_file "95d0b6027bd22592ac4ebde6412c0d63+100/$f" "$INPUTS/acmg/$f"
done

EXISTING_EXOMISER="/encrypted/e3008/Rob/yang-pangenome-full-cohort-20261012/tools/exomiser/data"
if [ -d "$EXISTING_EXOMISER/2512_hg38" ] && [ ! -d "$INPUTS/exomiser/data/2512_hg38" ]; then
  if mkdir "$INPUTS/exomiser/.lock-data" 2>/dev/null; then
    echo "[$(date -Is)] Copying Exomiser 2512 data from $EXISTING_EXOMISER -> $INPUTS/exomiser/data"
    cp -rf "$EXISTING_EXOMISER" "$INPUTS/exomiser/data.tmp.$$"
    mv -f "$INPUTS/exomiser/data.tmp.$$" "$INPUTS/exomiser/data"
    rmdir "$INPUTS/exomiser/.lock-data" 2>/dev/null || true
  fi
fi

if [ ! -d "$INPUTS/vep_cache/homo_sapiens/116_GRCh38/chr1" ]; then
  rmdir "$INPUTS/vep_cache/.lock-vep" 2>/dev/null || true
  if mkdir "$INPUTS/vep_cache/.lock-vep" 2>/dev/null; then
    echo "[$(date -Is)] Downloading VEP 116 cache from Keep (8afec8cdb796367144fa2bcec3a3dc59+624229) with 24 parallel workers -> $INPUTS/vep_cache/"
    "$PYTHON" - "$INPUTS/vep_cache" <<'PY'
import concurrent.futures, os, sys
import arvados, arvados.collection

dst_root = sys.argv[1]
cr = arvados.collection.CollectionReader("8afec8cdb796367144fa2bcec3a3dc59+624229")

items = []
def walk_col(col, rel=""):
    for k, v in col.items():
        p = f"{rel}/{k}" if rel else k
        if isinstance(v, arvados.collection.Subcollection):
            walk_col(v, p)
        else:
            items.append((p, v.size()))

walk_col(cr)
print(f"Total VEP cache files to check/stage: {len(items)}", flush=True)

def fetch_one(item):
    rel_path, expected_size = item
    out_path = os.path.join(dst_root, rel_path)
    if os.path.exists(out_path) and os.path.getsize(out_path) == expected_size:
        return 0
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    local_cr = arvados.collection.CollectionReader("8afec8cdb796367144fa2bcec3a3dc59+624229")
    tmp_path = f"{out_path}.part.{os.getpid()}"
    with local_cr.open(rel_path, "rb") as src, open(tmp_path, "wb") as dst:
        while True:
            buf = src.read(4 * 1024 * 1024)
            if not buf:
                break
            dst.write(buf)
    os.replace(tmp_path, out_path)
    return expected_size

with concurrent.futures.ThreadPoolExecutor(max_workers=24) as ex:
    total = sum(ex.map(fetch_one, items))
print(f"Staged {total / (1024**3):.2f} GiB into {dst_root}", flush=True)
PY
    rmdir "$INPUTS/vep_cache/.lock-vep" 2>/dev/null || true
  fi
fi

# Wait for any locked tail SIFs/inputs to finish if another job was staging them
for req_sif in \
  "nelly-exomiser:v4.sif" \
  "nelly-acmg-agent:v3.sif" \
  "ensemblorg_s_ensembl-vep:release___116.3.sif" \
  "brentp_s_slivar:v0.3.0.sif" \
  "dunglam2k_s_phenofrommetadata:v1.01.sif" \
  "staphb_s_bcftools:1.22.sif"; do
  while [ ! -s "$SIF_CACHE/$req_sif" ]; do
    echo "[$(date -Is)] Waiting for $SIF_CACHE/$req_sif..."
    sleep 5
  done
done
while [ ! -d "$INPUTS/exomiser/data/2512_hg38" ] || [ ! -d "$INPUTS/vep_cache/homo_sapiens/116_GRCh38" ]; do
  echo "[$(date -Is)] Waiting for exomiser/vep_cache..."
  sleep 5
done

if [ "${STAGE_ONLY:-0}" = "1" ]; then
  echo "[$(date -Is)] STAGE_ONLY=1 completed successfully!"
  exit 0
fi

# 3. Run full test-tail.cwl (including acmg-extract-candidates -> vep-candidates -> acmg-agent with Qwen) on TS and HD!
cd "$WF_DIR"
export HOST="ibex"
export NELLY_LOCAL_ROOT="$WORK_ROOT"

echo "[$(date -Is)] Starting full test-tail.cwl (with ACMG agent + Qwen) on TS..."
bash local/run.sh "$WORK_ROOT/runs/tail-TS-full" local/test-tail.cwl local/test-tail.TS.json

echo "[$(date -Is)] Starting full test-tail.cwl (with ACMG agent + Qwen) on HD..."
bash local/run.sh "$WORK_ROOT/runs/tail-HD-full" local/test-tail.cwl local/test-tail.HD.json

echo "[$(date -Is)] Both tail-TS-full and tail-HD-full completed successfully!"
