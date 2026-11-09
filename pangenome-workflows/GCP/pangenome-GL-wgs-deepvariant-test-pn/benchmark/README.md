# GIAB Benchmark Inputs Staging

This directory contains assets and scripts for benchmarking the modernized pangenome alignment and variant calling workflow against Genome in a Bottle (GIAB) reference materials (HG002 and HG003).

## Data Sources and Verified URLs

All URLs and file sizes verified as of **2026-11-09**.

### 1. HG002 30x NovaSeq PCR-Free FASTQ Reads
- **Source**: Google Brain Genomics DeepVariant training dataset
- **Files**:
  - `HG002.novaseq.pcr-free.30x.R1.fastq.gz`
    - URL: `https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R1.fastq.gz`
    - Size: 25,664,092,287 bytes (~25.7 GB)
  - `HG002.novaseq.pcr-free.30x.R2.fastq.gz`
    - URL: `https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R2.fastq.gz`
    - Size: 26,634,286,830 bytes (~26.6 GB)

### 2. HG002 GIAB Benchmark v4.2.1 (GRCh38)
- **Source**: NCBI GIAB FTP Release
- **Files**:
  - `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
    - Size: 156,252,944 bytes (~156 MB)
  - `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi`
    - Size: 1,657,747 bytes (~1.7 MB)
  - `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`
    - Size: 11,494,021 bytes (~11.5 MB)

### 3. HG002 GIAB Benchmark v5.0q (GRCh38)
- **Source**: NCBI GIAB FTP Release
- **Note on URL Structure**: Files are located directly under `/v5.0q/` (not `/v5.0q/GRCh38/`).
- **Files**:
  - `HG002_GRCh38_v5.0q_smvar.vcf.gz`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz`
    - Size: 42,067,644 bytes (~42 MB)
  - `HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi`
    - Size: 1,547,524 bytes (~1.5 MB)
  - `HG002_GRCh38_v5.0q_smvar.benchmark.bed`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed`
    - Size: 685,544 bytes (~686 KB)

### 4. HG003 GIAB Benchmark v4.2.1 (GRCh38)
- **Source**: NCBI GIAB FTP Release
- **Files**:
  - `HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
    - Size: 147,212,428 bytes (~147 MB)
  - `HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi`
    - Size: 1,651,307 bytes (~1.7 MB)
  - `HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`
    - Size: 13,219,856 bytes (~13.2 MB)

### 5. GIAB v3.5 Stratification BEDs (GRCh38 Union)
- **Source**: NCBI GIAB FTP Release
- **Note on URL Structure**: Files are located under `/genome-stratifications/v3.5/GRCh38@all/Union/`.
- **Files**:
  - `GRCh38_alldifficultregions.bed.gz`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz`
    - Size: 31,493,051 bytes (~31.5 MB)
  - `GRCh38_notinalldifficultregions.bed.gz`
    - URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_notinalldifficultregions.bed.gz`
    - Size: 32,740,449 bytes (~32.7 MB)

## HG003 Existing Keep Assets
HG003 reads and graph assets already reside in Keep collections:
- HG003 reads: `keep:4bc83e4703035a3eeb142eb17fcc1ce2+44786/HG003_R1.fastq.gz` and `HG003_R2.fastq.gz`
- Graph assets (KSA): `keep:1d5ad2aefb6e809dddaffea3a0558415+5657/{ksa_hg38.gbz,ksa_hg38.hapl,hg38.fa,hg38.fa.fai}`
- Legacy graph assets:
  - `dist`: `keep:1d5ad2aefb6e809dddaffea3a0558415+5657/ksa_hg38.dist`
  - `min`: `keep:56e94ed8e0d67d0c8e2616959d23ef2e+23423/ksa_hg38.min`

## Usage

To stage the datasets on workstation `ws`:
```bash
./benchmark/stage_giab_inputs.sh /mnt/data1/nelly-sim/giab
```

The script is idempotent: it validates existing files against server `Content-Length`, removes corrupted/error files (e.g. 996-byte 404 pages), and writes `SHA256SUMS`.

To upload staged assets to Arvados Keep:
```bash
arv-put --name giab-benchmark-inputs-2026-11 --project-uuid <PROJECT_UUID> /mnt/data1/nelly-sim/giab/hg002/* /mnt/data1/nelly-sim/giab/truth/* /mnt/data1/nelly-sim/giab/SHA256SUMS
```
