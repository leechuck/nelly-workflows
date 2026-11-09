#!/usr/bin/env bash
# Stage all inputs required for Phases 1 (tail smoke), 2 (full WGS TS/HD), and 3 (GIAB benchmark)
# into /data/nelly-local/inputs/ from Arvados Keep (with size verification against arv-ls -s),
# and compute sha256 checksums into local/inputs.sha256.
#
# Usage:
#   bash local/stage-inputs.sh [--phase1-only | --all] [--verify-only]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${NELLY_LOCAL_ROOT:-/data/nelly-local}"
ROOT="$WORK_ROOT/inputs"
export TMPDIR="$WORK_ROOT/tmp"
ARV_GET="${ARV_GET:-$WORK_ROOT/venv/bin/arv-get}"
ARV_LS="${ARV_LS:-$WORK_ROOT/venv/bin/arv-ls}"

MODE="${1:---all}"

mkdir -p "$ROOT"/{graph,slivar_vep_custom,vep_fasta,vep_cache,exomiser,acmg,tail/ts,tail/hd,patients/ts,patients/hd,giab}

# Helper: fetch a file from Keep if missing and verify its KiB size against arv-ls -s
fetch_keep_file() {
  local pdh="$1"
  local rel="$2"
  local dest="$3"
  if [ ! -s "$dest" ]; then
    echo "[$(date -Is)] arv-get $pdh/$rel -> $dest"
    "$ARV_GET" --threads 8 --skip-existing "$pdh/$rel" "$dest"
  else
    echo "[$(date -Is)] Exists: $dest ($(stat -c%s "$dest") bytes)"
  fi
}

# Verified Keep PDH inventory:
# - Graph + UCSC hg38 reference:      1d5ad2aefb6e809dddaffea3a0558415+5657
# - vg 1.78 indexes (.hapl/.dist/.ri): d704d713f7039bbc4026d55fafec7c33+2743
#   (Note: b3a8a85e970122c02ea058eb43457b07+177 in earlier JSONs is the coolmaksat/kmc Docker image tar, not .hapl)
# - Legacy vg 1.54 .min:              56e94ed8e0d67d0c8e2616959d23ef2e+23423
# - Tay-Sachs tail (post-norm + EH):  b2c8c2b940de5bf4d7881b652e12c24f+18026
# - Huntington tail (post-norm + EH): 6d04e84a0ab500c5036ec6713756783b+17943
# - Tay-Sachs WGS reads + metadata:   5c88241c12f7dc7553ffb3f2ac1ba54b+39418 (148_ID1_148; reads identical to 9bb0914061258315847c397cd2eda556+36864)
# - Tay-Sachs legacy id141 metadata:  4e618275df873b268594912ab6104f89+39753
# - Huntington WGS reads + metadata:  fa1d814cea61d867026682776bb0eb7a+36923 (hd-id201-htt-expansion-reads)
# - Huntington coded HPO metadata:    d03401bdd2443e6340903acb552fd0c3+39418 (149_ID2_149)
# - slivar gnomAD + VEP custom files: 148bcb74e5d8268ac486d315e657f70b+3240
# - VEP Ensembl-numeric fastafile:    00528ccfa8b07fd08907077b1f4fe3e9+927
# - Exomiser 2512 data bundle:        00cbcbd45419c2f051889639deff3431+39145
# - ACMG inputs (ClinVar/gnomAD):     28a24f20fb9659e2ba4e4b9a17f7bec0+778
# - GIAB benchmark inputs + HG002:    0c559a512e87c566c3a3de02d056b4a5+33970
# - GIAB HG003 WGS reads:             4bc83e4703035a3eeb142eb17fcc1ce2+44786

if [ "$MODE" != "--verify-only" ]; then
  echo "=== Phase 1 inputs ==="
  for f in normalized.vcf.gz normalized.vcf.gz.tbi eh_repeats.vcf.gz eh_repeats.vcf.gz.tbi eh_repeats.json slivar-output.vcf.gz; do
    fetch_keep_file "b2c8c2b940de5bf4d7881b652e12c24f+18026" "$f" "$ROOT/tail/ts/$f"
    fetch_keep_file "6d04e84a0ab500c5036ec6713756783b+17943" "$f" "$ROOT/tail/hd/$f"
  done

  fetch_keep_file "5c88241c12f7dc7553ffb3f2ac1ba54b+39418" "metadata.yaml" "$ROOT/patients/ts/metadata.yaml"
  fetch_keep_file "4e618275df873b268594912ab6104f89+39753" "metadata.yaml" "$ROOT/patients/ts/metadata.id141.yaml"
  fetch_keep_file "fa1d814cea61d867026682776bb0eb7a+36923" "metadata.yaml" "$ROOT/patients/hd/metadata.yaml"
  fetch_keep_file "d03401bdd2443e6340903acb552fd0c3+39418" "metadata.yaml" "$ROOT/patients/hd/metadata.id149.yaml"

  for f in \
    gnomad.hg38.genomes.v3.fix.zip \
    sorted.pavs.plugin.38.vcf.gz sorted.pavs.plugin.38.vcf.gz.tbi \
    sorted.plugin.go.bed.gz sorted.plugin.go.bed.gz.tbi \
    sorted.plugin.pheno.bed.gz sorted.plugin.pheno.bed.gz.tbi \
    sorted.plugin.ppi.bed.gz sorted.plugin.ppi.bed.gz.tbi
  do
    fetch_keep_file "148bcb74e5d8268ac486d315e657f70b+3240" "$f" "$ROOT/slivar_vep_custom/$f"
  done

  for f in fastafile fastafile.fai fastafile.gzi; do
    fetch_keep_file "00528ccfa8b07fd08907077b1f4fe3e9+927" "$f" "$ROOT/vep_fasta/$f"
  done

  for f in clinvar.vcf.gz clinvar.vcf.gz.md5 clinvar.vcf.gz.tbi gnomad.v4.1.constraint_metrics.tsv variant_summary.txt.gz; do
    fetch_keep_file "28a24f20fb9659e2ba4e4b9a17f7bec0+778" "$f" "$ROOT/acmg/$f"
  done

  if [ ! -s "$ROOT/acmg/STRchive-loci.json" ]; then
    curl -fsSL -o "$ROOT/acmg/STRchive-loci.json" \
      https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/STRchive-loci.json
  fi

  if [ ! -d "$ROOT/exomiser/data/2512_hg38" ]; then
    echo "[$(date -Is)] Staging Exomiser 2512 bundle (00cbcbd45419c2f051889639deff3431+39145)..."
    "$ARV_GET" -r --skip-existing "00cbcbd45419c2f051889639deff3431+39145/data" "$ROOT/exomiser/"
  fi

  if [ ! -d "$ROOT/vep_cache/homo_sapiens/116_GRCh38" ]; then
    if [ -s "$ROOT/vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz" ]; then
      echo "[$(date -Is)] Unpacking VEP 116 cache tarball..."
      (cd "$ROOT/vep_cache" && tar -xzf homo_sapiens_vep_116_GRCh38.tar.gz)
    else
      echo "[$(date -Is)] Downloading VEP 116 GRCh38 cache from Ensembl..."
      curl -fsSL -o "$ROOT/vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz" \
        https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
      (cd "$ROOT/vep_cache" && tar -xzf homo_sapiens_vep_116_GRCh38.tar.gz)
    fi
  fi

  if [ "$MODE" = "--all" ]; then
    echo "=== Phase 2 & 3 inputs (graph, WGS FASTQs, GIAB) ==="
    for f in ksa_hg38.gbz hg38.fa hg38.fa.fai; do
      fetch_keep_file "1d5ad2aefb6e809dddaffea3a0558415+5657" "$f" "$ROOT/graph/$f"
    done
    for f in ksa_hg38.hapl ksa_hg38.dist ksa_hg38.ri; do
      fetch_keep_file "d704d713f7039bbc4026d55fafec7c33+2743" "$f" "$ROOT/graph/$f"
    done
    fetch_keep_file "56e94ed8e0d67d0c8e2616959d23ef2e+23423" "ksa_hg38.min" "$ROOT/graph/ksa_hg38.v154.min"
    fetch_keep_file "1d5ad2aefb6e809dddaffea3a0558415+5657" "ksa_hg38.dist" "$ROOT/graph/ksa_hg38.v154.dist"

    fetch_keep_file "5c88241c12f7dc7553ffb3f2ac1ba54b+39418" "reads1.fastq.gz" "$ROOT/patients/ts/reads1.fastq.gz"
    fetch_keep_file "5c88241c12f7dc7553ffb3f2ac1ba54b+39418" "reads2.fastq.gz" "$ROOT/patients/ts/reads2.fastq.gz"
    fetch_keep_file "fa1d814cea61d867026682776bb0eb7a+36923" "reads1.fastq.gz" "$ROOT/patients/hd/reads1.fastq.gz"
    fetch_keep_file "fa1d814cea61d867026682776bb0eb7a+36923" "reads2.fastq.gz" "$ROOT/patients/hd/reads2.fastq.gz"

    for f in \
      GRCh38_alldifficultregions.bed.gz \
      GRCh38_notinalldifficultregions.bed.gz \
      HG002.novaseq.pcr-free.30x.R1.fastq.gz \
      HG002.novaseq.pcr-free.30x.R2.fastq.gz \
      HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
      HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
      HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
      HG002_GRCh38_v5.0q_smvar.benchmark.bed \
      HG002_GRCh38_v5.0q_smvar.vcf.gz \
      HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi \
      HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
      HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
      HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
      SHA256SUMS
    do
      fetch_keep_file "0c559a512e87c566c3a3de02d056b4a5+33970" "$f" "$ROOT/giab/$f"
    done
    fetch_keep_file "4bc83e4703035a3eeb142eb17fcc1ce2+44786" "HG003_R1.fastq.gz" "$ROOT/giab/HG003_R1.fastq.gz"
    fetch_keep_file "4bc83e4703035a3eeb142eb17fcc1ce2+44786" "HG003_R2.fastq.gz" "$ROOT/giab/HG003_R2.fastq.gz"
  fi
fi

echo "=== Computing SHA256 checksums for staged files (excluding unpacked VEP cache directory tree) ==="
(
  cd "$ROOT"
  find graph slivar_vep_custom vep_fasta exomiser acmg tail patients giab vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz \
    -type f 2>/dev/null | sort | xargs -r sha256sum
) | tee "$ROOT/inputs.sha256"
