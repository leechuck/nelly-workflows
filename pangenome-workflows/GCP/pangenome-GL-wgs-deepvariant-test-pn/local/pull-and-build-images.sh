#!/usr/bin/env bash
# Pull all public Docker images required by the modernised nelly-workflows pipeline
# and build the lab wrapper images (nelly-exomiser:v4 and nelly-acmg-agent:v3).
# Usage: bash local/pull-and-build-images.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
WORK_ROOT="${NELLY_LOCAL_ROOT:-/data/nelly-local}"
mkdir -p "$WORK_ROOT"/{images,tmp,logs}
export TMPDIR="$WORK_ROOT/tmp"

PUBLIC_IMAGES=(
  "staphb/bcftools:1.22"
  "brentp/slivar:v0.3.0"
  "ensemblorg/ensembl-vep:release_116.3"
  "dunglam2k/phenofrommetadata:v1.01"
  "dunglam2k/expansionhunter-grch38:v5.0.0"
  "coolmaksat/kmc"
  "coolmaksat/minimap2-samtools"
  "quay.io/vgteam/vg:v1.78.0"
  "google/deepvariant:1.10.0"
  "google/deepvariant:pangenome_aware_deepvariant-1.10.0"
  "jmcdani20/hap.py:v0.3.12"
)

echo "[$(date -Is)] Building nelly-exomiser:v4 from $PIPELINE_DIR/exomiser-wrapper/..."
DOCKER_BUILDKIT=0 docker build -t nelly-exomiser:v4 "$PIPELINE_DIR/exomiser-wrapper"

if [ -f "$PIPELINE_DIR/acmg-agent/Dockerfile" ]; then
  echo "[$(date -Is)] Building nelly-acmg-agent:v3 from $PIPELINE_DIR/acmg-agent/..."
  if docker image inspect nelly-acmg-agent:v2 >/dev/null 2>&1; then
    DOCKER_BUILDKIT=0 docker build -t nelly-acmg-agent:v3 -f - "$PIPELINE_DIR/acmg-agent" <<'EOF'
FROM nelly-acmg-agent:v2
COPY acmg_agent /app/acmg_agent
ENV PYTHONPATH=/app ACMG_MODEL_PATH=/model/model.gguf ACMG_LLM_BACKEND=llama
CMD ["python"]
EOF
  else
    DOCKER_BUILDKIT=0 docker build -t nelly-acmg-agent:v3 "$PIPELINE_DIR/acmg-agent"
  fi
fi

for img in "${PUBLIC_IMAGES[@]}"; do
  echo "[$(date -Is)] Pulling $img..."
  docker pull "$img"
done

echo "[$(date -Is)] Recording docker image digests to $WORK_ROOT/images/docker-images-digests.txt..."
docker images --digests | tee "$WORK_ROOT/images/docker-images-digests.txt"
