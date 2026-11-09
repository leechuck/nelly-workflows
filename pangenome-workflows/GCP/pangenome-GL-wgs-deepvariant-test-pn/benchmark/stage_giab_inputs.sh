#!/usr/bin/env bash
# stage_giab_inputs.sh
# Staging script for GIAB benchmark datasets (HG002, HG003, truth sets, stratifications).
# Runs on ws, downloads into /mnt/data1/nelly-sim/giab/
set -euo pipefail

GIAB_DIR="${1:-/mnt/data1/nelly-sim/giab}"
HG002_DIR="${GIAB_DIR}/hg002"
TRUTH_DIR="${GIAB_DIR}/truth"

mkdir -p "${HG002_DIR}" "${TRUTH_DIR}"

echo "=== GIAB Benchmark Data Staging ==="
echo "Target directory: ${GIAB_DIR}"

# Helper function to download idempotently and verify file size against Content-Length
download_and_verify() {
    local url="$1"
    local dest="$2"
    local expected_size="${3:-}"

    echo "==> Target: $(basename "${dest}")"

    # Query server content-length if not explicitly provided
    if [[ -z "${expected_size}" ]]; then
        expected_size=$(curl -sIL "${url}" | grep -i '^content-length:' | tail -n1 | tr -d '\r' | awk '{print $2}')
    fi

    if [[ -z "${expected_size}" ]]; then
        echo "ERROR: Failed to retrieve Content-Length from ${url}" >&2
        return 1
    fi

    # Check if another process is currently writing to dest
    if [[ -f "${dest}" ]] && command -v fuser >/dev/null 2>&1; then
        if fuser "${dest}" >/dev/null 2>&1; then
            echo "    [WAIT] Another process is actively writing to ${dest}. Waiting for it to finish..."
            while fuser "${dest}" >/dev/null 2>&1; do
                sleep 5
            done
        fi
    fi

    # Check existing file
    if [[ -f "${dest}" ]]; then
        local current_size
        current_size=$(stat -c %s "${dest}")
        if [[ "${current_size}" -eq "${expected_size}" ]]; then
            echo "    [SKIP] Exists and size matches (${current_size} bytes)."
            return 0
        elif [[ "${current_size}" -eq 996 ]]; then
            echo "    [WARN] File is 996 bytes (known 404 HTML error page). Removing."
            rm -f "${dest}"
        elif [[ "${current_size}" -gt "${expected_size}" ]]; then
            echo "    [WARN] File size (${current_size}) > expected (${expected_size}). Removing."
            rm -f "${dest}"
        else
            echo "    [RESUME] Partial file (${current_size} / ${expected_size} bytes). Resuming..."
        fi
    fi

    # Download with curl retry and resume
    curl -sSL --retry 5 -C - -o "${dest}" "${url}"

    # Verify downloaded size
    local final_size
    final_size=$(stat -c %s "${dest}")
    if [[ "${final_size}" -ne "${expected_size}" ]]; then
        echo "ERROR: Size mismatch for ${dest}! Expected: ${expected_size}, got: ${final_size}" >&2
        return 1
    fi
    echo "    [OK] Verified size (${final_size} bytes)."
}

echo ""
echo "--- 1. HG002 30x NovaSeq PCR-free reads ---"
download_and_verify \
    "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R1.fastq.gz" \
    "${HG002_DIR}/HG002.novaseq.pcr-free.30x.R1.fastq.gz" \
    "25664092287"

download_and_verify \
    "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R2.fastq.gz" \
    "${HG002_DIR}/HG002.novaseq.pcr-free.30x.R2.fastq.gz" \
    "26634286830"

echo ""
echo "--- 2. HG002 GIAB v4.2.1 benchmark truth set (GRCh38) ---"
download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "${TRUTH_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "156252944"

download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
    "${TRUTH_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
    "1657747"

download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "${TRUTH_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "11494021"

echo ""
echo "--- 3. HG002 GIAB v5.0q benchmark truth set (GRCh38) ---"
download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz" \
    "${TRUTH_DIR}/HG002_GRCh38_v5.0q_smvar.vcf.gz" \
    "42067644"

download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi" \
    "${TRUTH_DIR}/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi" \
    "1547524"

download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed" \
    "${TRUTH_DIR}/HG002_GRCh38_v5.0q_smvar.benchmark.bed" \
    "685544"

echo ""
echo "--- 4. HG003 GIAB v4.2.1 benchmark truth set (GRCh38) ---"
download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "${TRUTH_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "147212428"

download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
    "${TRUTH_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
    "1651307"

download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "${TRUTH_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "13219856"

echo ""
echo "--- 5. GIAB v3.5 Stratification BEDs (GRCh38 Union) ---"
download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz" \
    "${TRUTH_DIR}/GRCh38_alldifficultregions.bed.gz" \
    "31493051"

download_and_verify \
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_notinalldifficultregions.bed.gz" \
    "${TRUTH_DIR}/GRCh38_notinalldifficultregions.bed.gz" \
    "32740449"

echo ""
echo "--- 6. Computing SHA256 checksums ---"
(
    cd "${GIAB_DIR}"
    sha256sum \
        hg002/HG002.novaseq.pcr-free.30x.R1.fastq.gz \
        hg002/HG002.novaseq.pcr-free.30x.R2.fastq.gz \
        truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
        truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        truth/HG002_GRCh38_v5.0q_smvar.vcf.gz \
        truth/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi \
        truth/HG002_GRCh38_v5.0q_smvar.benchmark.bed \
        truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
        truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        truth/GRCh38_alldifficultregions.bed.gz \
        truth/GRCh38_notinalldifficultregions.bed.gz \
        > "${GIAB_DIR}/SHA256SUMS"
)
echo "Wrote ${GIAB_DIR}/SHA256SUMS."

echo ""
echo "=== Staging Complete ==="
echo "To upload staged datasets to Keep, run:"
echo "arv-put --name giab-benchmark-inputs-2026-11 --project-uuid <left blank for supervisor> ${HG002_DIR}/* ${TRUTH_DIR}/* ${GIAB_DIR}/SHA256SUMS"
