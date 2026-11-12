#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --job-name=nelly-stage-ibex
#SBATCH --output=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-stage-%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/nelly-local/logs/slurm-stage-%j.err

# Stage all Singularity .sif images (public + custom Keep docker-archive tars)
# and all workflow inputs onto /ibex/scratch/projects/c2014/nelly-local/ inside Slurm.
set -euo pipefail

WORK_ROOT="${NELLY_LOCAL_ROOT:-/ibex/scratch/projects/c2014/nelly-local}"
INPUTS="$WORK_ROOT/inputs"
SIF_CACHE="$WORK_ROOT/sif_cache"
IMAGES="$WORK_ROOT/images"
export TMPDIR="$WORK_ROOT/tmp/stage-${SLURM_JOB_ID:-$$}"
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

# 1. Stage SIF images
EXISTING_TOOLS="/encrypted/e3008/Rob/yang-pangenome-full-cohort-20261012/modernisation-20261108/tools"
if [ -s "$EXISTING_TOOLS/dv110.sif" ]; then
  stage_sif_copy "google_s_deepvariant:1.10.0.sif" "$EXISTING_TOOLS/dv110.sif"
fi
if [ -s "$EXISTING_TOOLS/padv110.sif" ]; then
  stage_sif_copy "google_s_deepvariant:pangenome___aware___deepvariant-1.10.0.sif" "$EXISTING_TOOLS/padv110.sif"
fi

for uri in \
  "quay.io/vgteam/vg:v1.78.0" \
  "coolmaksat/kmc" \
  "coolmaksat/minimap2-samtools" \
  "dunglam2k/expansionhunter-grch38:v5.0.0" \
  "staphb/bcftools:1.22" \
  "dunglam2k/phenofrommetadata:v1.01" \
  "brentp/slivar:v0.3.0" \
  "ensemblorg/ensembl-vep:release_116.3" \
  "jmcdani20/hap.py:v0.3.12"; do
  stage_sif_pull "$uri"
done

stage_sif_from_keep_tar "nelly-exomiser:v4.sif" "7712a129b6c89fb29f6f39eb5842ad90+429/sha256:1e800dce4081b13665db09bcb3e650f8fe38d241ab2eea2f2031a205ecda6dcd.tar"
stage_sif_from_keep_tar "nelly-acmg-agent:v2.sif" "c15a4a032ac6d0bcb2c10e8bd9653a05+1731/sha256:8ed50b5e81a5b87c4d05aca26b7dccf818c76c2f0afe125a74680e4f667fe341.tar"
ln -sf "nelly-acmg-agent:v2.sif" "$SIF_CACHE/nelly-acmg-agent:v1.sif"

# 2. Stage graph + reference + patient FASTQs
GRAPH_PDH="1d5ad2aefb6e809dddaffea3a0558415+5657"
VG170_IDX_PDH="d704d713f7039bbc4026d55fafec7c33+2743"
stage_keep_file "$GRAPH_PDH/ksa_hg38.gbz" "$INPUTS/graph/ksa_hg38.gbz"
stage_keep_file "$VG170_IDX_PDH/ksa_hg38.hapl" "$INPUTS/graph/ksa_hg38.hapl"
stage_keep_file "$GRAPH_PDH/hg38.fa" "$INPUTS/graph/hg38.fa"
stage_keep_file "$GRAPH_PDH/hg38.fa.fai" "$INPUTS/graph/hg38.fa.fai"

for f in metadata.yaml reads1.fastq.gz reads2.fastq.gz; do
  stage_keep_file "5c88241c12f7dc7553ffb3f2ac1ba54b+39418/$f" "$INPUTS/patients/ts/$f"
  stage_keep_file "d03401bdd2443e6340903acb552fd0c3+39418/$f" "$INPUTS/patients/hd/$f"
done

# 3. Stage tail inputs (slivar, vep_fasta, exomiser, acmg, tail VCFs, vep_cache)
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
  echo "[$(date -Is)] Copying Exomiser 2512 data from $EXISTING_EXOMISER -> $INPUTS/exomiser/data"
  mkdir -p "$INPUTS/exomiser"
  cp -rf "$EXISTING_EXOMISER" "$INPUTS/exomiser/data.tmp.$$"
  mv -f "$INPUTS/exomiser/data.tmp.$$" "$INPUTS/exomiser/data"
fi

if [ ! -d "$INPUTS/vep_cache/homo_sapiens/116_GRCh38" ]; then
  echo "[$(date -Is)] Downloading VEP 116 cache from Keep (8afec8cdb796367144fa2bcec3a3dc59+624229) -> $INPUTS/vep_cache/"
  "$ARV_GET" -f --threads 16 "8afec8cdb796367144fa2bcec3a3dc59+624229/homo_sapiens" "$INPUTS/vep_cache/homo_sapiens"
fi

echo "[$(date -Is)] IBEX staging complete!"
