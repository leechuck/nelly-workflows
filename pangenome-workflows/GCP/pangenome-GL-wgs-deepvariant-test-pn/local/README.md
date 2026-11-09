# Local Execution & Cloud Deployment Guide (`local/`)

This directory contains everything needed to run the modernised `nelly-workflows` pangenome pipeline (`main-vg.cwl`, `local/test-tail.cwl`, and `benchmark/benchmark-giab.cwl`) on a single host with `cwltool` + Docker (exercised on `cbontsr01.kaust.edu.sa` / ontolinator, 256 threads, 1 TB RAM), or to provision an AWS / Oracle Cloud deployment.

## 1. Host Prerequisites & Directory Layout

All local paths live under `/data/nelly-local/` (override via `NELLY_LOCAL_ROOT`), with `TMPDIR=/data/nelly-local/tmp` (never `/tmp`):

```text
/data/nelly-local/
├── venv/                  # Python 3.11 venv (cwltool 3.3+, arvados-python-client)
├── inputs/                # Staged reference, graph, caches, patient & GIAB inputs
│   ├── graph/             # ksa_hg38.gbz, ksa_hg38.hapl (v1.78), ksa_hg38.dist, ksa_hg38.ri, hg38.fa, hg38.fa.fai
│   ├── slivar_vep_custom/ # gnomad.hg38.genomes.v3.fix.zip + 4 custom VEP VCF/BED plugins (+ .tbi)
│   ├── vep_fasta/         # Ensembl-numeric fastafile (+ .fai, .gzi)
│   ├── vep_cache/         # homo_sapiens/116_GRCh38/
│   ├── exomiser/          # data/2512_hg38/ + data/2512_phenotype/
│   ├── acmg/              # clinvar.vcf.gz(+.tbi), gnomad.v4.1.constraint_metrics.tsv, variant_summary.txt.gz, gene_table.tsv, STRchive-loci.json
│   ├── tail/{ts,hd}/      # Pre-computed normalized.vcf.gz(+.tbi), eh_repeats.vcf.gz(+.tbi), eh_repeats.json
│   ├── patients/{ts,hd}/  # reads1.fastq.gz, reads2.fastq.gz, metadata.yaml
│   ├── giab/              # HG002 & HG003 FASTQs + GIAB v4.2.1 / v5.0q truth VCFs/BEDs + stratifications
│   └── inputs.sha256      # SHA-256 manifest of all staged files
├── images/                # Docker build context & recorded image digests
├── runs/                  # Per-run output directories (incl. step-resources.tsv/.md)
├── logs/                  # Background job logs
└── tmp/                   # TMPDIR for cwltool and container scratch
```

### Credentials (mode `0600`, never committed)
1. **Arvados Keep** (for staging from `cborg.projectnelly.com`): `~/.config/arvados/settings.conf`
   ```ini
   ARVADOS_API_HOST=cborg.projectnelly.com
   ARVADOS_API_TOKEN=<token>
   ```
2. **Unimatrix Qwen vLLM endpoint** (for `acmg-agent`): `~/.config/borg-llm/api-key`
   Contains the raw Bearer API key on a single line. Committed JSON templates set `"llm_api_key": "FROM-ENV"`, and `local/run.sh` injects the key into a mode-600 runtime-only JSON in `$TMPDIR` that is deleted on exit.

## 2. Step-by-Step Setup on a Fresh Box

```bash
cd pangenome-workflows/GCP/pangenome-GL-wgs-deepvariant-test-pn

# 1. Create /data/nelly-local/venv (Python 3.11 + cwltool 3.x + arvados-python-client)
bash local/setup-venv.sh

# 2. Pull the 11 public Docker images and build nelly-exomiser:v4 (+ nelly-acmg-agent:v1)
setsid nohup bash local/pull-and-build-images.sh > /data/nelly-local/logs/stage-images.log 2>&1 < /dev/null &

# 3. Stage all inputs from Arvados Keep and verify checksums
setsid nohup bash local/stage-inputs.sh --all > /data/nelly-local/logs/stage-inputs.log 2>&1 < /dev/null &
```

## 3. Verified Arvados Keep Collections (Single Source of Truth)

| Bundle | Verified Keep PDH | Key Files | Local Destination |
|---|---|---|---|
| Graph + UCSC GRCh38 ref | `1d5ad2aefb6e809dddaffea3a0558415+5657` | `ksa_hg38.gbz`, `hg38.fa`, `hg38.fa.fai` | `inputs/graph/` |
| vg 1.78 indexes | `d704d713f7039bbc4026d55fafec7c33+2743` | `ksa_hg38.hapl`, `ksa_hg38.dist`, `ksa_hg38.ri` | `inputs/graph/` |
| Legacy vg 1.54 `.min` | `56e94ed8e0d67d0c8e2616959d23ef2e+23423` | `ksa_hg38.min` | `inputs/graph/ksa_hg38.v154.min` |
| slivar gnomAD + VEP custom | `148bcb74e5d8268ac486d315e657f70b+3240` | `gnomad.hg38.genomes.v3.fix.zip`, 4 custom VCF/BED + `.tbi` | `inputs/slivar_vep_custom/` |
| VEP Ensembl-numeric FASTA | `00528ccfa8b07fd08907077b1f4fe3e9+927` | `fastafile`, `fastafile.fai`, `fastafile.gzi` | `inputs/vep_fasta/` |
| VEP 116 GRCh38 cache | `vep-116-GRCh38-cache` (Ensembl 116) | `homo_sapiens/116_GRCh38/` | `inputs/vep_cache/` |
| Exomiser 2512 bundle | `00cbcbd45419c2f051889639deff3431+39145` | `data/2512_hg38/`, `data/2512_phenotype/` | `inputs/exomiser/` |
| ACMG databases | `28a24f20fb9659e2ba4e4b9a17f7bec0+778` | `clinvar.vcf.gz`(+.tbi), `gnomad.v4.1.constraint_metrics.tsv`, `variant_summary.txt.gz` | `inputs/acmg/` |
| Tay-Sachs tail inputs | `b2c8c2b940de5bf4d7881b652e12c24f+18026` | `normalized.vcf.gz`(+.tbi), `eh_repeats.vcf.gz`(+.tbi), `eh_repeats.json` | `inputs/tail/ts/` |
| Huntington tail inputs | `6d04e84a0ab500c5036ec6713756783b+17943` | `normalized.vcf.gz`(+.tbi), `eh_repeats.vcf.gz`(+.tbi), `eh_repeats.json` | `inputs/tail/hd/` |
| Tay-Sachs WGS reads + metadata | `5c88241c12f7dc7553ffb3f2ac1ba54b+39418` (`9bb0914061258315847c397cd2eda556+36864`) | `reads1.fastq.gz`, `reads2.fastq.gz`, `metadata.yaml` | `inputs/patients/ts/` |
| Huntington WGS reads + metadata | `fa1d814cea61d867026682776bb0eb7a+36923` (`d03401bdd2443e6340903acb552fd0c3+39418`) | `reads1.fastq.gz`, `reads2.fastq.gz`, `metadata.yaml` | `inputs/patients/hd/` |
| GIAB truth + HG002 reads | `0c559a512e87c566c3a3de02d056b4a5+33970` | `HG002.novaseq.pcr-free.30x.R{1,2}.fastq.gz`, v4.2.1 & v5.0q truths, stratifications | `inputs/giab/` |
| GIAB HG003 reads | `4bc83e4703035a3eeb142eb17fcc1ce2+44786` | `HG003_R1.fastq.gz`, `HG003_R2.fastq.gz` | `inputs/giab/` |

> **Note on Keep PDHs:** Earlier templates referenced `b3a8a85e970122c02ea058eb43457b07+177` for `ksa_hg38.hapl`; that PDH is actually the `coolmaksat/kmc` Docker image tar. The real `ksa_hg38-vg1.78-indexes` collection in Keep is `d704d713f7039bbc4026d55fafec7c33+2743`.

## 4. Running Phases 1, 2, and 3

All runs go through `local/run.sh`, which sets `TMPDIR=/data/nelly-local/tmp`, injects `llm_api_key` into a temporary mode-600 input JSON, launches `local/sample-resources.py` to record per-step container wall time and peak RSS (`<outdir>/step-resources.tsv` and `<outdir>/step-resources.md`), and invokes `cwltool --parallel --timestamps`.

### Phase 1: Tail-Only Smoke (`local/test-tail.cwl`)
```bash
# Tay-Sachs tail
setsid nohup bash local/run.sh \
  /data/nelly-local/runs/tail-TS \
  local/test-tail.cwl \
  local/test-tail.TS.json \
  > /data/nelly-local/logs/tail-TS.log 2>&1 < /dev/null &

# Huntington tail
setsid nohup bash local/run.sh \
  /data/nelly-local/runs/tail-HD \
  local/test-tail.cwl \
  local/test-tail.HD.json \
  > /data/nelly-local/logs/tail-HD.log 2>&1 < /dev/null &
```

### Phase 2: Full WGS End-to-End (`main-vg.cwl`)
```bash
# Tay-Sachs (standard DeepVariant 1.10)
setsid nohup bash local/run.sh \
  /data/nelly-local/runs/wgs-TS-standard \
  main-vg.cwl \
  local/main-vg.TS.standard.json \
  > /data/nelly-local/logs/wgs-TS-standard.log 2>&1 < /dev/null &

# Huntington (standard DeepVariant 1.10)
setsid nohup bash local/run.sh \
  /data/nelly-local/runs/wgs-HD-standard \
  main-vg.cwl \
  local/main-vg.HD.standard.json \
  > /data/nelly-local/logs/wgs-HD-standard.log 2>&1 < /dev/null &
```

### Phase 3: GIAB Benchmark (`benchmark/benchmark-giab.cwl`)
```bash
for cfg in HG003.standard HG003.pangenome_aware HG002.standard HG002.pangenome_aware; do
  setsid nohup bash local/run.sh \
    "/data/nelly-local/runs/giab-$cfg" \
    benchmark/benchmark-giab.cwl \
    "local/benchmark-giab.$cfg.json" \
    > "/data/nelly-local/logs/giab-$cfg.log" 2>&1 < /dev/null &
done
```
