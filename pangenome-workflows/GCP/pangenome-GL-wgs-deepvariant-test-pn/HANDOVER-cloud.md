# Cloud Deployment & Migration Handover — Modernised Nelly Workflow

**Audience:** Dung Lam (Workflow Owner & Arvados Cluster Operator)  
**Date:** 2026-11-09  
**Repository:** `nelly-workflows` (`pangenome-workflows/GCP/pangenome-GL-wgs-deepvariant-test-pn/`)  
**Target Environments:** AWS / Oracle Cloud Infrastructure (OCI) / On-premise Arvados clusters  
**Accompanying Machine-Readable Tables:** `handover/` (`data-bundles.tsv`, `container-images.tsv`, `workflow-inputs-main-vg.tsv`, `giab-sha256sums.txt`)

---

## 1. What Changed

Between 2026-11-08 and 2026-11-09, the entire alignment, variant calling, annotation, and interpretation stack of `nelly-workflows` was modernised to match the November 2026 Yang et al. Saudi Pangenome benchmark specification and replace the ungrounded `genome-linter` LLM narrative with a verified, deterministic ACMG/AMP 2015 evaluation agent.

```
OLD STACK (October 2026 cborg runs):
  FASTQ -> vg Giraffe 1.54.0 (37 GB .min + .dist) -> DeepVariant 1.6.0 -> ExpansionHunter 5.0.0
        -> bcftools 1.19 -> slivar 0.3.0 -> VEP 102.0 -> Exomiser 15.0.0 (data 2512)
        -> genome-linter-llm (OpenRouter / local llama-cpp blurb.py)

NEW STACK (November 2026 Modernisation):
  FASTQ -> KMC 29-mer counts (.kff) -> vg Giraffe 1.78.0 (Haplotype Sampling, .hapl)
        -> samtools 1.22 sort/markdup -> DeepVariant 1.10.0 OR PADV 1.10.0 (dv_mode switch)
        -> ExpansionHunter 5.0.0 -> bcftools 1.22 concat/norm -> exclude-eh-regions
        -> slivar 0.3.0 -> VEP 116.3 (Pass 1: slivar VCF) -> Exomiser 15.1.0 (data 2512)
        -> acmg-extract-candidates -> VEP 116.3 (Pass 2: candidates mini-VCF)
        -> acmg-agent (Deterministic 28-criteria ClinGen SVI engine + Bounded LLM Reviewer)
```

### 1.1 Summary of Core Changes

| Component | Pre-Modernisation | Modernised | Key Rationale & Technical Notes |
|---|---|---|---|
| **Graph Aligner** | `vg:1.54.0` | `quay.io/vgteam/vg:v1.78.0` | Switched to **haplotype sampling** (`--haplotype-name ksa_hg38.hapl --kff-name sample.kff --index-basename personal`). Drops the monolithic 37.3 GB `.min` minimizer index; builds personalized graph indexes in memory per sample. |
| **K-mer Counter** | None | `coolmaksat/kmc` (latest) | Generates exact 29-mer frequency counts (`sample.kff`, 590 MB) from raw FASTQ reads in <6 s for haplotype sampling. |
| **Graph Indexes** | `ksa_hg38.dist`, `ksa_hg38.min` (v1.54 format) | `ksa_hg38.hapl`, `ksa_hg38.dist`, `ksa_hg38.ri` (v1.78 format) | **Critical:** vg 1.54 `.hapl`/`.dist` files in Keep are version 2 format and crash vg 1.78 with `Expected version 5 to 6, got version 2`. All indexes were rebuilt with vg 1.78 on `ws` and staged to `ksa_hg38-vg1.78-indexes` (`d704d713f7039bbc4026d55fafec7c33+2743`). |
| **Variant Caller** | `google/deepvariant:1.6.0` | `google/deepvariant:1.10.0` OR `pangenome_aware_deepvariant-1.10.0` | Configurable via `dv_mode: standard \| pangenome_aware`. Standard mode calls against GRCh38; PADV calls directly against the GBZ pangenome graph. |
| **PADV /dev/shm Handling** | N/A | Modified stub invocation in `pangenome-aware-deepvariant.cwl` | Arvados `crunch-run` and standard Docker allocate 64 MB `/dev/shm`. Default PADV attempts to load `ksa_hg38.gbz` (~467 MB) into `/dev/shm/GBZ_SHARED_MEMORY` and crashes with `Bus error` (`SIGBUS`). The CWL stubs `load_gbz_into_shared_memory_command` and passes `use_loaded_gbz_shared_memory=false`. Each shard loads the graph into heap via `GbzReader` (12.7 GiB/shard); `num_shards` is set to 4 (40.94 GiB peak) to stay well within 90 GiB. |
| **VCF Utilities** | `staphb/bcftools:1.19` | `staphb/bcftools:1.22` | Upgraded `bcftools-concat.cwl`, `bcftools-norm.cwl`, and `exclude-eh-regions.cwl`. |
| **VEP Annotation** | `ensembl-vep:release_102.0` | `ensembl-vep:release_116.3` | Pinned `--cache_version 116`; added `--mane --protein --pick_allele_gene --af_gnomade --af_gnomadg --max_af`. Dropped `--af_gnomad` (superseded by `--af_gnomade`; VEP 116 fails if both are passed). Cache staged to `vep-116-GRCh38-cache` (`8afec8cdb796367144fa2bcec3a3dc59+624229`, 26 GB). |
| **Two-Pass VEP** | Single pass on slivar output | Two passes: Pass 1 (slivar) + Pass 2 (`vep-candidates`) | Pass 1 annotates slivar-filtered variants for Exomiser. `acmg-extract-candidates` extracts Exomiser top candidates and repeat expansions into a mini-VCF (`candidates.vcf.gz`, ~20–50 records), and Pass 2 annotates it with full VEP 116 CSQ fields to guarantee transcript and frequency annotations for every candidate. |
| **Exomiser** | `nelly-exomiser:v3` (15.0.0) | `nelly-exomiser:v4` (15.1.0) | Distribution downloaded at image build time from GitHub release (no pre-unpacked build context needed). Base image updated to `eclipse-temurin:21-jre`. Output remains `results/exomiser.jsonl`. Exomiser data bundle remains release **2512**. |
| **Candidate Interpretation** | `genome-linter-llm.cwl` (`blurb.py` + local/remote LLM) | `acmg-agent.cwl` (`nelly-acmg-agent:v2`) | Deterministic 28-criteria ACMG/AMP 2015 + ClinGen SVI engine coupled with STRchive repeat expansion tiers. Includes a bounded LLM reviewer (OpenAI-compatible vLLM endpoint or local CPU GGUF fallback). Fail-closed contract enforced in code. |
| **Sub-Workflow Architecture** | Monolithic `main-vg.cwl` | `align-and-call.cwl` + `main-vg.cwl` + `benchmark-giab.cwl` | Core alignment and calling logic factored out into `align-and-call.cwl`, shared identically between production patient runs and GIAB benchmark runs. |
| **CWL Standard** | `cwlVersion: v1.1` | `cwlVersion: v1.2` | Bumped across all active CWL tools and workflows to support conditional execution (`when:`) for `dv_mode`. |

---

## 2. Workflow Entry Points and Inputs

The pipeline provides three primary workflow entry points and one legacy baseline:
1. `main-vg.cwl`: Production patient WGS workflow (alignment $\rightarrow$ calling $\rightarrow$ filtering $\rightarrow$ annotation $\rightarrow$ Exomiser $\rightarrow$ ACMG Agent).
2. `align-and-call.cwl`: Shared sub-workflow (KMC $\rightarrow$ vg Giraffe $\rightarrow$ samtools $\rightarrow$ DV/PADV + ExpansionHunter).
3. `benchmark/benchmark-giab.cwl`: Benchmarking workflow (`align-and-call` $\rightarrow$ `hap.py` against GIAB truth).
4. `benchmark/benchmark-giab-legacy.cwl`: Baseline benchmark running the pre-modernisation stack (vg 1.54 + DV 1.6).

```
Entry Points Overview:

[ reads1, reads2, graph, hapl, ref ]
                |
                v
      +--------------------+
      | align-and-call.cwl | <--------------------+
      +--------------------+                      |
         |               |                        |
(CRAM, DV VCF, EH JSON)  (DV VCF)                 |
         |               |                        |
         v               v                        |
  +-------------+  +----------------------+       |
  | main-vg.cwl |  | benchmark-giab.cwl   |       |
  +-------------+  +----------------------+       |
  (Exomiser +      (hap.py vcfeval vs             |
   ACMG Agent)      GIAB HG002/HG003)             |
                                                  |
  +-----------------------------+                 |
  | benchmark-giab-legacy.cwl   | ----------------+ (uses legacy/*.cwl:
  +-----------------------------+                      vg 1.54 + DV 1.6)
```

### 2.1 `main-vg.cwl` (Production Clinical Workflow)

- **Example JSONs in Repository:**
  - `main-vg-pn.json` (NIST small test reads template)
  - `main-vg.TS.e2e.json` (Tay-Sachs clinical patient `id: 148` e2e test)
  - `main-vg.HD.e2e.json` (Huntington clinical patient `id: 149` e2e test)

| Input Name | CWL Type | Default Value | Source Bundle / Recommended Keep Path | Description |
|---|---|---|---|---|
| `reads1` | `File` | *(required)* | `Patient Collection/reads1.fastq.gz` | Forward paired-end FASTQ reads |
| `reads2` | `File` | *(required)* | `Patient Collection/reads2.fastq.gz` | Reverse paired-end FASTQ reads |
| `metadata` | `File` | *(required)* | `Patient Collection/metadata.yaml` | Genoor YAML with sex and HPO phenotype URIs |
| `graph` | `File` | *(required)* | `1d5ad2aefb6e809dddaffea3a0558415+5657/ksa_hg38.gbz` | Saudi pangenome graph (GBZ) |
| `hapl` | `File` | *(required)* | `d704d713f7039bbc4026d55fafec7c33+2743/ksa_hg38.hapl` | vg 1.78 haplotype sampling index |
| `ref` | `File` (sec: `.fai`) | *(required)* | `1d5ad2aefb6e809dddaffea3a0558415+5657/hg38.fa` | UCSC chr-named GRCh38 reference FASTA + index |
| `ref_prefix` | `string` | `"GRCh38"` | N/A | Prefix for reference paths in graph |
| `sample_name` | `string` | `"sample"` | N/A | Sample identifier in BAM/VCF headers |
| `dv_mode` | `string` | `"standard"` | `"standard"` or `"pangenome_aware"` | Small variant caller selection |
| `dv_model_type`| `string` | `"WGS"` | N/A | DeepVariant model type |
| `aligned_bam_output` | `string`| `"aligned.bam"` | N/A | Intermediate BAM output filename |
| `dv_output_vcf`| `string` | `"variants.dv.vcf.gz"`| N/A | DeepVariant output VCF filename |
| `dv_output_gvcf`| `string` | `"variants.dv.gvcf.gz"`| N/A | DeepVariant output gVCF filename |
| `eh_variant_catalog` | `File?` | `null` | Optional override; defaults to image catalog | ExpansionHunter repeat catalog |
| `eh_output_prefix` | `string`| `"eh_repeats"` | N/A | ExpansionHunter output prefix |
| `combined_vcf_name` | `string`| `"combined.vcf.gz"` | N/A | bcftools-concat output filename |
| `normalized_vcf_name`| `string`| `"normalized.vcf.gz"`| N/A | bcftools-norm output filename |
| `slivar_gnomad`| `File` | *(required)* | `148bcb74e5d8268ac486d315e657f70b+3240/gnomad.hg38.genomes.v3.fix.zip` | gnomAD v3 genome frequencies for slivar |
| `slivar_info` | `string` | `"INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt <= 3 && variant.FILTER == 'PASS'"` | Default | Population frequency & quality gate |
| `slivar_sample_expr`| `string` | `"high_quality:sample.GQ >= 20 && sample.DP >= 10 && (sample.het \|\| sample.hom_alt)"` | Default | Sample genotype quality gate |
| `vep_assembly`| `string` | `"GRCh38"` | N/A | Assembly name for VEP |
| `vep_output_file`| `string`| `"vep-output.vcf.gz"`| N/A | VEP Pass 1 output filename |
| `vep_dir` | `Directory` | *(required)* | `8afec8cdb796367144fa2bcec3a3dc59+624229` | Ensembl VEP 116 indexed cache directory |
| `vep_fasta_file`| `File` | *(required)* | `00528ccfa8b07fd08907077b1f4fe3e9+927/fastafile` | Ensembl numeric-chromosome FASTA for VEP |
| `pavs_custom_file`| `File` (sec: `.tbi`)| *(required)* | `148bcb74e5d8268ac486d315e657f70b+3240/sorted.pavs.plugin.38.vcf.gz` | PAVS plugin annotation VCF |
| `pavs_custom_args`| `string`| `"PAVS,vcf,exact"` | Default | PAVS plugin argument string |
| `go_custom_file` | `File` (sec: `.tbi`)| *(required)* | `148bcb74e5d8268ac486d315e657f70b+3240/sorted.plugin.go.bed.gz` | Gene Ontology plugin BED |
| `go_custom_args` | `string`| `"GO_CLASSES,bed,overlap"` | Default | GO plugin argument string |
| `hpo_custom_file`| `File` (sec: `.tbi`)| *(required)* | `148bcb74e5d8268ac486d315e657f70b+3240/sorted.plugin.pheno.bed.gz` | Phenotype plugin BED |
| `hpo_custom_args`| `string`| `"PHENOTYPE,bed,overlap"` | Default | HPO plugin argument string |
| `ppi_custom_file`| `File` (sec: `.tbi`)| *(required)* | `148bcb74e5d8268ac486d315e657f70b+3240/sorted.plugin.ppi.bed.gz` | Protein-protein interaction plugin BED |
| `ppi_custom_args`| `string`| `"PPI,bed,overlap"` | Default | PPI plugin argument string |
| `exomiser_data`| `Directory` | *(required)* | `00cbcbd45419c2f051889639deff3431+39145` | Exomiser 2512 data directory (`data/`) |
| `exomiser_assembly`| `string`| `"hg38"` | N/A | Exomiser reference assembly |
| `exomiser_data_version`| `string`| `"2512"` | N/A | Exomiser data bundle version |
| `hpo_terms` | `string` | `""` | Optional fallback | Curated comma-separated HPO terms |
| `exomiser_top_genes`| `int` | `10` | Default | Number of top Exomiser genes |
| `acmg_clinvar_vcf`| `File` (sec: `.tbi`)| *(required)* | `28a24f20fb9659e2ba4e4b9a17f7bec0+778/clinvar.vcf.gz` | ClinVar GRCh38 VCF + index (2026-11-01) |
| `acmg_gene_table` | `File` | *(required)* | `95d0b6027bd22592ac4ebde6412c0d63+100/gene_table.tsv` | Gene constraint & ClinVar composition |
| `acmg_strchive_json`| `File?` | `null` | `95d0b6027bd22592ac4ebde6412c0d63+100/STRchive-loci.json` | STRchive 82-locus repeat thresholds |
| `acmg_top_candidates`| `int` | `10` | Default | Number of candidate genes reviewed |
| `llm_backend` | `string` | `"openai"` | `"openai"`, `"llama"`, or `"off"` | ACMG agent reviewer backend |
| `llm_endpoint` | `string` | `"http://unimatrix01.kaust.edu.sa:8000/v1"` | N/A | OpenAI-compatible endpoint URL |
| `llm_model` | `string` | `"qwen3.8-27b"` | N/A | Model identifier for completions |
| `llm_api_key` | `string` | `""` | Inject at runtime (`"FROM-ENV"`) | Bearer token for LLM API |
| `llm_max_tokens` | `int` | `4096` | Default | Token budget per completion turn |

- **Workflow Outputs:**
  - `output_aligned_cram`: Final indexed CRAM file (`.cram` + `.crai`).
  - `output_dv_vcf`: Small variant calls from DeepVariant (`.vcf.gz` + `.tbi`).
  - `output_dv_gvcf`: Genomic VCF from DeepVariant (`.gvcf.gz` + `.tbi`).
  - `output_eh_vcf`: ExpansionHunter repeat variant calls (`.vcf.gz` + `.tbi`).
  - `output_eh_json`: Detailed ExpansionHunter repeat analysis (`.json`).
  - `output_combined_vcf`: Combined normalized variant calls (`.vcf.gz`).
  - `output_slivar`: Rare/high-quality filtered VCF (`.vcf.gz`).
  - `output_after_vep`: Pass 1 genome-wide annotated VCF (`.vcf.gz`).
  - `output_exomiser_genes`: Prioritised genes table (`exomiser.genes.tsv`).
  - `output_exomiser_variants`: Prioritised variants table (`exomiser.variants.tsv`).
  - `output_exomiser_html`: Interactive Exomiser HTML report.
  - `output_candidates_vep`: Pass 2 annotated candidate mini-VCF (`candidates.vep.vcf.gz`).
  - `output_acmg_json`: Machine-readable structured ACMG evaluation (`acmg-assessment.json`, schema v1).
  - `output_acmg_report`: Grounded human-readable ACMG clinical report (`acmg-report.md`).

### 2.2 `align-and-call.cwl` (Sub-Workflow)

- **Purpose:** Standalone, reusable sub-workflow containing the compute-intensive alignment and caller steps.
- **Inputs:** `reads1`, `reads2`, `graph`, `hapl`, `ref` (+ `.fai`), `ref_prefix`, `sample_name`, `dv_mode` (`"standard" | "pangenome_aware"`), `aligned_bam_output`, `dv_model_type`, `dv_output_vcf`, `dv_output_gvcf`, `eh_variant_catalog`, `eh_sex`, `eh_output_prefix`.
- **Outputs:** `aligned_cram` (with `.crai`), `dv_vcf` (with `.tbi`), `dv_gvcf` (with `.tbi`), `eh_vcf` (with `.tbi`), `eh_json`.

### 2.3 `benchmark/benchmark-giab.cwl` & `benchmark-giab-legacy.cwl`

- **Purpose:** Validates variant caller accuracy against NIST Genome in a Bottle (GIAB) truth sets using GA4GH `hap.py` with the `vcfeval` engine.
- **Example JSONs in Repository:**
  - `benchmark/benchmark-giab.HG002.standard.json`
  - `benchmark/benchmark-giab.HG002.pangenome_aware.json`
  - `benchmark/benchmark-giab.HG003.standard.json`
  - `benchmark/benchmark-giab.HG003.pangenome_aware.json`
  - `benchmark/benchmark-giab-legacy.HG003.json`
- **Inputs Specific to Benchmarking:**
  - `truth_vcf` (`File`, sec: `.tbi`): GIAB benchmark truth VCF (e.g. `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`).
  - `confident_bed` (`File`): Confident call regions BED (e.g. `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`).
  - `stratification_tsv` (`File?`): Optional stratification TSV (`GRCh38_alldifficultregions.bed.gz` / `notinalldifficultregions`).
  - `stratification_beds` (`File[]?`): Array of stratification BED files staged into the workdir.
  - `happy_prefix` (`string`, default `"happy"`): Prefix for hap.py report tables.
- **Outputs:**
  - `happy_summary_csv`: Summary recall, precision, and F1 scores (`.summary.csv`).
  - `happy_extended_csv`: Per-stratification detailed metrics (`.extended.csv`).
  - `happy_metrics_json`: Compressed metrics object (`.metrics.json.gz`).
  - `happy_roc_all`: Array of ROC curve files (`.roc.*.csv.gz`).

### 2.4 Critical Configuration Semantics

#### `llm_api_key: "FROM-ENV"`
Committed template JSONs set `"llm_api_key": "FROM-ENV"` to prevent credentials from leaking into Git.
- **How runners must handle this:**
  - The submission wrapper or runner must read the API key from a secure local source (such as `~/.config/borg-llm/api-key` or an environment variable `BORG_LLM_API_KEY`) and substitute it into the runtime JSON before execution.
  - In `acmg-agent.cwl`, `llm_api_key` is passed exclusively via `EnvVarRequirement` to `ACMG_LLM_API_KEY`. **The API key never appears in `baseCommand`, `arguments`, or process `argv` lists**, ensuring it is not logged in Arvados Crunch logs or process tables.
  - If `llm_backend` is set to `"llama"` (local CPU GGUF) or `"off"`, `llm_api_key` is ignored.

#### `dv_mode: standard | pangenome_aware`
- `dv_mode: "standard"`: Executes `deepvariant.cwl` (`google/deepvariant:1.10.0`) on the GRCh38-surjected BAM.
- `dv_mode: "pangenome_aware"`: Executes `pangenome-aware-deepvariant.cwl` (`google/deepvariant:pangenome_aware_deepvariant-1.10.0`), passing `--pangenome $(inputs.graph.path)` and running without `/dev/shm`.
- The selection is handled dynamically via CWL v1.2 `when:` directives in `align-and-call.cwl`.

---

## 3. Data Bundles (The Big Table)

Every data bundle referenced in the workflow input JSONs has been inspected, checksummed, and verified directly from Arvados Keep (`cborg.projectnelly.com`).

```
Data Storage Strategy:
Keep Collection (KAUST cborg)
   |  arv-get -r <PDH>/ <local_dir>/
   v
Local Staging Disk
   |  aws s3 cp / oci os object bulk-upload
   v
Cloud Bucket (S3 / OCI Object Storage)
```

| Bundle Name | Keep PDH | Keep UUID | Total Size | Files | Key Contents | Provenance & Version | Licence / Access | Cloud Migration Command |
|---|---|---|---:|---:|---|---|---|---|
| **KSA_GRCh38_pangenome** | `1d5ad2aefb6e809dddaffea3a0558415+5657` | `cborg-4zz18-eh24zxxo8cjx9mj` | 8.00 GiB (8,586,792,607 B) | 6 | `ksa_hg38.gbz`, `hg38.fa`, `hg38.fa.fai`, `ksa_hg38.dist`, `ksa_hg38.snarls`, `ksa_hg38.hapl` | Saudi Pangenome Project (Yang et al.) / UCSC GRCh38 | Academic / Open Research | `arv-get -r 1d5ad2aefb6e809dddaffea3a0558415+5657/ ksa_graph/ && aws s3 cp ksa_graph/ s3://<BUCKET>/inputs/graph/ --recursive` |
| **ksa_hg38-vg1.78-indexes** | `d704d713f7039bbc4026d55fafec7c33+2743` | `cborg-4zz18-g4wjy136eb53bcl` | 3.82 GiB (4,100,880,256 B) | 3 | `ksa_hg38.hapl` (1.24 GB), `ksa_hg38.dist` (1.41 GB), `ksa_hg38.ri` (1.46 GB) | Rebuilt on `ws` with vg 1.78.0 from `ksa_hg38.gbz` (2026-11-09) | Academic / Open Research | `arv-get -r d704d713f7039bbc4026d55fafec7c33+2743/ ksa_vg178/ && aws s3 cp ksa_vg178/ s3://<BUCKET>/inputs/graph/ --recursive` |
| **vep-116-GRCh38-cache** | `8afec8cdb796367144fa2bcec3a3dc59+624229` | `cborg-4zz18-n2p8nd24l66til8` | 25.77 GiB (27,673,136,989 B) | 13,922 | `homo_sapiens/116_GRCh38/` (indexed tabix cache) | Ensembl Release 116 (`homo_sapiens_vep_116_GRCh38.tar.gz`) | Ensembl Terms (Open EMBL-EBI) | `arv-get -r 8afec8cdb796367144fa2bcec3a3dc59+624229/ vep_cache/ && aws s3 cp vep_cache/ s3://<BUCKET>/inputs/vep_cache/ --recursive` |
| **vep102_indexed_fasta** | `00528ccfa8b07fd08907077b1f4fe3e9+927` | `cborg-4zz18-mlr7xwr91tnic76` | 1.12 GiB (1,202,819,800 B) | 3 | `fastafile`, `fastafile.fai`, `fastafile.gzi` | Ensembl numeric-chromosome GRCh38 FASTA with bgzip/gzi index | Ensembl Terms | `arv-get -r 00528ccfa8b07fd08907077b1f4fe3e9+927/ vep_fasta/ && aws s3 cp vep_fasta/ s3://<BUCKET>/inputs/vep_fasta/ --recursive` |
| **exomiser-2512-hg38-phenotype-bundle** | `00cbcbd45419c2f051889639deff3431+39145` | `cborg-4zz18-2msbuc55s0g3ote` | 53.99 GiB (57,969,064,336 B) | 51 | `data/2512_hg38/` (variants, genome, clinvar), `data/2512_phenotype/` | Monarch Initiative Exomiser release 2512 | GPLv3 / Monarch (OMIM & HPO conditions apply) | `arv-get -r 00cbcbd45419c2f051889639deff3431+39145/ exomiser/ && aws s3 cp exomiser/ s3://<BUCKET>/inputs/exomiser/ --recursive` |
| **slivar_vep_fixed_files** | `148bcb74e5d8268ac486d315e657f70b+3240` | `cborg-4zz18-8vf6dsklaefglmt` | 3.81 GiB (4,089,558,370 B) | 9 | `gnomad.hg38.genomes.v3.fix.zip` (4.08 GB), `sorted.pavs.plugin.38.vcf.gz` (+.tbi), `sorted.plugin.{go,pheno,ppi}.bed.gz` (+.tbi) | gnomAD v3.1 + Lab custom annotation plugins | gnomAD Terms / Lab Internal | `arv-get -r 148bcb74e5d8268ac486d315e657f70b+3240/ slivar_vep/ && aws s3 cp slivar_vep/ s3://<BUCKET>/inputs/slivar_vep_custom/ --recursive` |
| **acmg-inputs-2026-11** | `28a24f20fb9659e2ba4e4b9a17f7bec0+778` | `cborg-4zz18-nao9ae036ben5hy` | 0.70 GiB (755,367,936 B) | 5 | `clinvar.vcf.gz` (201 MB), `clinvar.vcf.gz.tbi`, `gnomad.v4.1.constraint_metrics.tsv` (95 MB), `variant_summary.txt.gz` (458 MB) | NCBI ClinVar (release 2026-11-01) + gnomAD v4.1 constraint | Public Domain (NCBI) / gnomAD Terms | `arv-get -r 28a24f20fb9659e2ba4e4b9a17f7bec0+778/ acmg_raw/ && aws s3 cp acmg_raw/ s3://<BUCKET>/inputs/acmg/ --recursive` |
| **acmg-gene-table-strchive-2026-11** | `95d0b6027bd22592ac4ebde6412c0d63+100` | `cborg-4zz18-lurk4nfxtwis1ap` | 1.22 MiB (1,251,754 B) | 2 | `gene_table.tsv` (725 KB, 19,642 genes), `STRchive-loci.json` (526 KB, 82 loci) | Derived from gnomAD v4.1 + ClinVar; STRchive repo data (2026-11-09) | MIT (STRchive) / CC0 (Gene Table) | `arv-get -r 95d0b6027bd22592ac4ebde6412c0d63+100/ acmg_gt/ && aws s3 cp acmg_gt/ s3://<BUCKET>/inputs/acmg/ --recursive` |
| **giab-benchmark-inputs-2026-11** | `0c559a512e87c566c3a3de02d056b4a5+33970` | `cborg-4zz18-bjyz7jqz4o7v4pf` | 49.12 GiB (52,738,403,110 B) | 14 | HG002 30x FASTQs R1/R2 (52.3 GB), GIAB v4.2.1 & v5.0q truth sets, stratifications, `SHA256SUMS` | NIST GIAB / Google brain-genomics-public bucket | Public Domain (NIST) | `arv-get -r 0c559a512e87c566c3a3de02d056b4a5+33970/ giab/ && aws s3 cp giab/ s3://<BUCKET>/inputs/giab/ --recursive` |
| **HG003_fastq** | `4bc83e4703035a3eeb142eb17fcc1ce2+44786` | `cborg-4zz18-mvpgzp27vjy4gx6` | 66.12 GiB (70,998,487,794 B) | 5 | `HG003_R1.fastq.gz` (32.0 GB), `HG003_R2.fastq.gz` (33.5 GB), WES reads, capture BED | NIST HG003 NovaSeq PCR-free 35x | Public Domain (NIST) | `arv-get -r 4bc83e4703035a3eeb142eb17fcc1ce2+44786/ hg003/ && aws s3 cp hg003/ s3://<BUCKET>/inputs/giab_hg003/ --recursive` |
| **148_ID1_148 (Tay-Sachs Patient)** | `5c88241c12f7dc7553ffb3f2ac1ba54b+39418` | `cborg-4zz18-dig5cye2hpmyk3x` | 54.71 GiB (58,743,059,923 B) | 4 | `reads1.fastq.gz` (28.9 GB), `reads2.fastq.gz` (29.8 GB), `metadata.yaml` (`id: 148`), `pedigree.ped` | Patient WGS cohort (Tay-Sachs, HEXA hom) | Protected Health Data (KAUST IRB/IBEC) | `arv-get -r 5c88241c12f7dc7553ffb3f2ac1ba54b+39418/ ts/ && aws s3 cp ts/ s3://<SECURE_BUCKET>/patients/ts/ --recursive` |
| **149_ID2_149 (Huntington Patient)** | `d03401bdd2443e6340903acb552fd0c3+39418` | `cborg-4zz18-2p9btihl28rxn12` | 54.71 GiB (58,740,292,603 B) | 4 | `reads1.fastq.gz` (28.9 GB), `reads2.fastq.gz` (29.8 GB), `metadata.yaml` (`id: 149`), `pedigree.ped` | Patient WGS cohort (Huntington, HTT 18/46) | Protected Health Data (KAUST IRB/IBEC) | `arv-get -r d03401bdd2443e6340903acb552fd0c3+39418/ hd/ && aws s3 cp hd/ s3://<SECURE_BUCKET>/patients/hd/ --recursive` |
| **ksa_hg38.min (Legacy)** | `56e94ed8e0d67d0c8e2616959d23ef2e+23423` | `cborg-4zz18-3jzzem82ma335rf` | 34.78 GiB (37,348,221,952 B) | 1 | `ksa_hg38.min` (minimizer index for vg 1.54) | Built with vg 1.54.0 (Legacy baseline only) | Academic / Open Research | `arv-get -r 56e94ed8e0d67d0c8e2616959d23ef2e+23423/ ksa_min/ && aws s3 cp ksa_min/ s3://<BUCKET>/inputs/graph_legacy/ --recursive` |
| **nist-small (Smoke Reads)** | `b4a8fbf918d6744687ac1b775456ca79+536` | `cborg-4zz18-8a29d7wv4sjujn3` | 0.63 GiB (679,479,388 B) | 2 | `NIST-small_R1.fastq.gz`, `NIST-small_R2.fastq.gz` | NIST test slice for quick dry-runs | Public Domain | `arv-get -r b4a8fbf918d6744687ac1b775456ca79+536/ test_reads/ && aws s3 cp test_reads/ s3://<BUCKET>/inputs/test_reads/ --recursive` |
| **phenofrommetadata_sample_yaml** | `3c9215d6f3acbd11de88305e3187dfc0+57` | `cborg-4zz18-oycqmrllepshtnz` | 887 B | 1 | `sample.yaml` | Sample metadata template | Open Research | `arv-get -r 3c9215d6f3acbd11de88305e3187dfc0+57/ sample_yaml/ && aws s3 cp sample_yaml/ s3://<BUCKET>/inputs/sample_yaml/ --recursive` |

> **Note on Keep PDH Correction:** Early planning notes referenced `b3a8a85e970122c02ea058eb43457b07+177` for `ksa_hg38.hapl`; that collection is actually the Docker image tarball for `coolmaksat/kmc:latest`. The true vg 1.78 index collection is `d704d713f7039bbc4026d55fafec7c33+2743`. All repository JSONs have been updated with this verified PDH.

### 3.1 Verified SHA-256 Checksums for GIAB Benchmark Bundle (`0c559a512e87c566c3a3de02d056b4a5+33970`)

```text
b8d5cadbc01371cd808c455d4edfdb99ecdede117c5c466532cee02762af7c1a  HG002.novaseq.pcr-free.30x.R1.fastq.gz
3c9775cd614a59fffbe20e6d2aa454fa0a70e2ee14574cfe3cdafe79eec9da66  HG002.novaseq.pcr-free.30x.R2.fastq.gz
e6119692ce74779459466c2617563b1da5677cbec7ec4dd5219f08d43a262930  GRCh38_alldifficultregions.bed.gz
ea7c8b0220f0e12bae9feeec8ad62b411bc5a32d34298c7ac9c1b6897f3ec668  GRCh38_notinalldifficultregions.bed.gz
fba9a57c36ec88e5d14ea3e259c8866c7935f4998d3ec0fa2d6c3962da5b5575  HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
23162712553413b10f88916c28bd1238deb703a7237be2910440773beea71894  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
ef7d898dd421ea357bb9db915894ab3e5c6c1689e36ccf619befda4bbf1cef5b  HG002_GRCh38_v5.0q_smvar.benchmark.bed
c7f9d7a4781f1f0a8a42ec54b051d6b95db1935c5e12318bae219740f9c50dc8  HG002_GRCh38_v5.0q_smvar.vcf.gz
4d77f6eb74eebf676572c70d13a0b353e3e8c81468f250739c91ad1343102eb2  HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi
652afd3046705af3200f9c87c255fef11bb212dd76c75a19999c9b2df8a3180c  HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
a4f8fefa826f8d2c9457eacf6f4b2f6bfa29daef429d570ce04254c5e7b1121e  HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
989305947e5af2e81ca8ad1087262a12cc87499b36163c49a7878dc9e18fcaa5  HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
```

---

## 4. Container Images

Every image has been pulled, tested, and registered in Arvados Keep (`arv-keepdocker`). For cloud deployment outside Arvados, images must be mirrored to AWS ECR or Oracle Cloud Infrastructure Registry (OCIR).

### 4.1 Master Container Image Table

| CWL `dockerPull` Tag | Role / Steps | Source | Dockerfile Path | Keepdocker UUID | Keepdocker PDH | Size | Cloud Migration Action |
|---|---|---|---|---|---|---:|---|
| `quay.io/vgteam/vg:v1.78.0` | `vg-paths`, `vg-giraffe-bam` | Public (Quay) | N/A | `cborg-4zz18-b6fx6hspms8ll2b` | `b85cadbf0e8e397cd54bb5fdc809b57e+555` | 667 MB | Re-tag & push to ECR/OCIR |
| `coolmaksat/kmc` | `kmc` (k-mer counting) | Public (Docker Hub) | N/A | `cborg-4zz18-9wm3htobwgzoyge` | `b3a8a85e970122c02ea058eb43457b07+177` | 110 MB | Re-tag & push to ECR/OCIR |
| `coolmaksat/minimap2-samtools` | `samtools-sort`, `samtools-index`, `samtools-markdup` | Public (Docker Hub) | N/A | `cborg-4zz18-upmi5lh8vv7catl` | `85d30a63c9111d644e3c5a9135a05981+512` | 587 MB | Re-tag & push to ECR/OCIR |
| `google/deepvariant:1.10.0` | `deepvariant` (standard mode) | Public (Docker Hub) | N/A | `cborg-4zz18-hntgy76r3e4rn0p` | `37f180d007ba596a54ac44fdac056cb3+3244` | 4.91 GB | Re-tag & push to ECR/OCIR |
| `google/deepvariant:pangenome_aware_deepvariant-1.10.0` | `pangenome-aware-deepvariant` | Public (Docker Hub) | N/A | `cborg-4zz18-io6ryhidk8xhhbs` | `1089704d143bd81268c6772f1bbdade4+4546` | 6.99 GB | Re-tag & push to ECR/OCIR |
| `dunglam2k/expansionhunter-grch38:v5.0.0` | `expansionhunter` | Lab (Docker Hub) | `Dockerfile.expansionhunter` | `cborg-4zz18-qmmxl6516t6hbqn` | `27d7db33349d4971804bbd5d7e88362a+303` | 316 MB | Re-tag & push to ECR/OCIR |
| `staphb/bcftools:1.22` | `bcftools-concat`, `bcftools-norm`, `exclude-eh-regions` | Public (Docker Hub) | N/A | `cborg-4zz18-6c21q30t2xtav3z` | `5519a9ad66b0f99461ca00c46fcf4c75+219` | 146 MB | Re-tag & push to ECR/OCIR |
| `dunglam2k/phenofrommetadata:v1.01` | `phenofrommetadata` | Lab (Docker Hub) | N/A | `cborg-4zz18-552s6289g8hyre5` | `0670dcc33497ad6a43459b1505b98583+766` | 1.01 GB | Re-tag & push to ECR/OCIR |
| `brentp/slivar:v0.3.0` | `slivar` | Public (Docker Hub) | N/A | `cborg-4zz18-yralzaojmwfvskm` | `3bb6ba071c0663b7c2147c629c34826e+934` | 1.26 GB | Re-tag & push to ECR/OCIR |
| `ensemblorg/ensembl-vep:release_116.3` | `vep`, `vep-candidates` | Public (Docker Hub) | N/A | `cborg-4zz18-ss114g8sfa98kjr` | `e491cb49fdb30932e8adb51e3404c2c2+722` | 923 MB | Re-tag & push to ECR/OCIR |
| `nelly-exomiser:v4` | `exomiser` | Lab-built | `exomiser-wrapper/Dockerfile` | `cborg-4zz18-fjxxofhos7m3h0w` | `7712a129b6c89fb29f6f39eb5842ad90+429` | 513 MB | Build locally & push to ECR/OCIR |
| `nelly-acmg-agent:v2` | `acmg-extract-candidates`, `acmg-agent` | Lab-built | `acmg-agent/Dockerfile` | `cborg-4zz18-p2b223ia8ie9968` | `c15a4a032ac6d0bcb2c10e8bd9653a05+1731` | 1.80 GB | Build locally & push to ECR/OCIR |
| `jmcdani20/hap.py:v0.3.12` | `happy` (GIAB benchmark) | Public (Docker Hub) | N/A | `cborg-4zz18-p91vg2azx1aeu8f` | `0ffa70e9fbc7cd1f05a3a58790828d54+1102` | 1.56 GB | Re-tag & push to ECR/OCIR (benchmark only) |
| `quay.io/biocontainers/vg:1.54.0--h9ee0642_0`| `legacy/vg-giraffe-bam` | Public (Quay) | N/A | `cborg-4zz18-rq5fanhkpxfwei4` | `effbc63fa34c2b3f82f960d22063ea70+134` | 140 MB | Legacy baseline only |
| `google/deepvariant:1.6.0` | `legacy/deepvariant` | Public (Docker Hub) | N/A | `cborg-4zz18-zao654gh59tjwkh` | `5e98b419e3fd99e15e6d555da2ddc6f0+3748` | 3.80 GB | Legacy baseline only |

### 4.2 Exact Build Commands for Lab-Built Images

Both lab images strictly follow **CWL Gotcha #4**: `ENTRYPOINT` is unset (`null`), and execution defaults to `CMD ["bash"]` or `CMD ["python"]`.

#### 1. `nelly-exomiser:v4`
Downloads `exomiser-cli-15.1.0-distribution.zip` from GitHub at build time; uses `eclipse-temurin:21-jre`.
```bash
cd pangenome-workflows/GCP/pangenome-GL-wgs-deepvariant-test-pn

docker build -t nelly-exomiser:v4 exomiser-wrapper/

# Verify:
docker run --rm nelly-exomiser:v4 java -jar /exomiser/exomiser-cli-15.1.0.jar --version
# Expected output: 15.1.0
```

#### 2. `nelly-acmg-agent:v2`
Builds `acmg_agent` Python package, compiles `llama-cpp-python 0.3.2` with `CMAKE_ARGS=-DGGML_NATIVE=OFF` (Debian compatibility, Gotcha #18), and bakes in `Qwen2.5-3B-Instruct-Q4_K_M.gguf` at `/model/model.gguf` as the offline fallback.
```bash
cd pangenome-workflows/GCP/pangenome-GL-wgs-deepvariant-test-pn

docker build -t nelly-acmg-agent:v2 acmg-agent/

# Verify:
docker run --rm nelly-acmg-agent:v2 python -m acmg_agent.cli --help
# Expected output: lists 'extract' and 'assess' subcommands
```

### 4.3 Cloud Registry Mirroring (AWS ECR / OCI OCIR)
```bash
# Example for AWS ECR:
ACCOUNT_ID="123456789012"
REGION="us-east-1"
REGISTRY="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/nelly"

aws ecr get-login-password --region ${REGION} | docker login --username AWS --password-stdin ${REGISTRY}

for img in \
  "nelly-acmg-agent:v2" \
  "nelly-exomiser:v4" \
  "quay.io/vgteam/vg:v1.78.0" \
  "google/deepvariant:1.10.0" \
  "google/deepvariant:pangenome_aware_deepvariant-1.10.0" \
  "ensemblorg/ensembl-vep:release_116.3" \
  "staphb/bcftools:1.22" \
  "dunglam2k/expansionhunter-grch38:v5.0.0" \
  "dunglam2k/phenofrommetadata:v1.01" \
  "coolmaksat/kmc:latest" \
  "coolmaksat/minimap2-samtools:latest" \
  "brentp/slivar:v0.3.0" \
  "jmcdani20/hap.py:v0.3.12"
do
  # Extract bare tag
  bare=$(echo ${img} | awk -F'/' '{print $NF}')
  docker tag ${img} ${REGISTRY}:${bare}
  docker push ${REGISTRY}:${bare}
done
```

---

## 5. Compute Sizing & Resource Requirements

### 5.1 Per-Step Resource Allocations and Observed Measurements

Observed metrics combine:
- Historic cborg runs (26–27 Oct 2026, 10-core / 90 GiB nodes).
- Workstation `ws` smoke measurements (56 cores, 125 GiB RAM).
- Phase 1 local tail smoke runs on `cbontsr01` (256 threads, 1 TB RAM, cgroup resource monitoring).
- WP4b live assessment of `nelly-acmg-agent:v2` against `unimatrix01` `qwen3.8-27b`.

| Step Name | Tool CWL | Docker Image | `coresMin` | `ramMin` (MiB) | `tmpdirMin` (MB) | Observed Wall Time | Observed Peak Process RSS | Peak Cgroup Footprint | Sizing Note |
|---|---|---|---:|---:|---:|---|---:|---:|---|
| `phenofrommetadata` | `phenofrommetadata.cwl` | `dunglam2k/phenofrommetadata:v1.01` | 1 | 1,024 | — | 1.1 s | <15 MiB | <50 MiB | Lightweight YAML parser |
| `kmc` | `kmc.cwl` | `coolmaksat/kmc` | 8 | 65,536 | 200,000 | 5.7 s (2M reads) / ~12 min (30x WGS) | 3.8 GiB | 4.2 GiB | Memory scales with k-mer buffer |
| `vg-paths` | `vg-paths.cwl` | `quay.io/vgteam/vg:v1.78.0` | 10 | 36,864 | — | ~45 s | 2.1 GiB | 2.5 GiB | Extracts reference path list from GBZ |
| `vg-giraffe-bam` | `vg-giraffe-bam.cwl` | `quay.io/vgteam/vg:v1.78.0` | 10 | 92,160 | 300,000 | 2,350 s (idx) + 14 h (30x map) | **74.24 GiB** | 78.5 GiB | **Pipeline Peak RAM**. Haplotype sampling building personal graph |
| `samtools-sort` | `samtools-sort.cwl` | `coolmaksat/minimap2-samtools` | 8 | 16,384 | — | ~25 min (30x WGS) | 3.5 GiB | 4.1 GiB | Uses BAM stream |
| `samtools-index-sorted`| `samtools-index.cwl`| `coolmaksat/minimap2-samtools` | 2 | 4,096 | — | ~3 min | 0.8 GiB | 1.2 GiB | Indexing sorted BAM |
| `samtools-markdup` | `samtools-markdup.cwl`| `coolmaksat/minimap2-samtools` | 8 | 16,384 | — | ~35 min (30x WGS) | 4.2 GiB | 4.8 GiB | Duplicate marking |
| `samtools-markdup-index`| `samtools-index.cwl`| `coolmaksat/minimap2-samtools` | 2 | 4,096 | — | ~3 min | 0.8 GiB | 1.2 GiB | Final CRAM index |
| `deepvariant` | `deepvariant.cwl` | `google/deepvariant:1.10.0` | 8 | 92,160 | — | 291 s (chr20:1-10M) / ~17 h (30x WGS) | 5.25 GiB | 6.1 GiB | Standard WGS mode |
| `pangenome-aware-deepvariant` | `pangenome-aware-deepvariant.cwl` | `google/deepvariant:pangenome_aware_deepvariant-1.10.0` | 8 | 92,160 | 100,000 | 376 s (chr20:1-10M) / ~18–20 h (30x WGS) | **40.94 GiB** | 43.9 GiB | 4 shards × ~10.2–12.7 GiB each (no `/dev/shm`) |
| `expansionhunter` | `expansionhunter.cwl` | `dunglam2k/expansionhunter-grch38:v5.0.0`| 8 | 16,384 | — | ~30 min | 3.8 GiB | 4.5 GiB | Repeat expansion profiling |
| `bcftools-concat` | `bcftools-concat.cwl` | `staphb/bcftools:1.22` | 4 | 8,192 | — | ~45 s | 120 MiB | 350 MiB | Concatenates DV and EH calls |
| `bcftools-norm` | `bcftools-norm.cwl` | `staphb/bcftools:1.22` | 4 | 8,192 | — | ~12 min | 250 MiB | 600 MiB | Multiallelic splitting and left-alignment |
| `exclude-eh-regions` | `exclude-eh-regions.cwl`| `staphb/bcftools:1.22` | 2 | 4,096 | — | 8.4–9.4 s | 12.7 MiB | 204.8 MiB | Strips repetitive regions before Exomiser |
| `slivar` | `slivar.cwl` | `brentp/slivar:v0.3.0` | 1 | 8,192 | — | 167.9 s | 2.53 GiB | 4.86 GiB | Evaluates gnomAD popmax frequency |
| `vep` (Pass 1) | `vep.cwl` | `ensemblorg/ensembl-vep:release_116.3` | 2 | 8,192 | — | 605 s solo (719 s conc.) | 860.6 MiB | 12.16 GiB (incl cache) | Annotates ~1,000 candidate variants |
| `exomiser` | `exomiser.cwl` | `nelly-exomiser:v4` | 4 | 32,768 | — | 75.7–76.8 s | 3.58 GiB | 4.49 GiB | Phenotype + variant pathogenicity ranking |
| `acmg-extract-candidates`| `acmg-extract-candidates.cwl`| `nelly-acmg-agent:v2` | 1 | 2,048 | — | ~3 s | 85 MiB | 150 MiB | Generates candidate mini-VCF (~20–50 vars) |
| `vep-candidates` (Pass 2)| `vep.cwl` | `ensemblorg/ensembl-vep:release_116.3` | 2 | 8,192 | — | ~20 s | 550 MiB | 2.5 GiB | Fast consequence annotation on mini-VCF |
| `acmg-agent` | `acmg-agent.cwl` | `nelly-acmg-agent:v2` | 8 | 20,480 | — | 1,437 s (~24 min) via endpoint | 1.2 GiB | 2.1 GiB | Full 28-criteria ClinGen review + LLM |
| `happy` (GIAB benchmark) | `benchmark/happy.cwl` | `jmcdani20/hap.py:v0.3.12` | 8 | 32,768 | — | ~18 min | 8.5 GiB | 12.0 GiB | vcfeval comparison against GIAB truth |
| **Full WGS End-to-End** | `main-vg.cwl` | Multi-container | 10 | 92,160 | 500,000 | *TODO (WP6 phase 2)* | **74.24 GiB** | 80.0 GiB | Currently running on IBEX (jobs 55405248/49) |

### 5.2 Recommended Cloud Instance Types

To satisfy the **74.2 GiB RAM** ceiling imposed by `vg giraffe` haplotype sampling and avoid Out-Of-Memory (OOM) failures:

#### Amazon Web Services (AWS)
- **Single-Node Execution (`cwltool` / Toil) or Worker Nodes (Arvados):**
  - **Recommended Family:** Memory Optimized `r6i` (Intel Ice Lake) or `r7i` (Intel Sapphire Rapids).
  - **Recommended Shape:** `r6i.4xlarge` or `r7i.4xlarge` (16 vCPUs, 128 GiB RAM, up to 12.5 Gbps network).
  - **Storage:** Amazon EBS GP3 volume: minimum **1,000 GB** (1 TB), provisioned at 3,000 IOPS and 250 MB/s throughput to accommodate Giraffe and DeepVariant I/O.
  - **Cost-Optimized Alternative:** `r6a.4xlarge` (16 vCPUs, 128 GiB RAM, AMD EPYC 7003 series).

#### Oracle Cloud Infrastructure (OCI)
- **Flexible Virtual Machines:**
  - **Recommended Shape:** `VM.Standard.E5.Flex` (AMD EPYC Genoa) or `VM.Standard.E4.Flex`.
  - **Allocation:** **16 OCPUs** (equivalent to 32 vCPUs) and **128 GB RAM**.
  - **Storage:** Block Volume: **1,000 GB**, Balanced or Higher Performance (30–75 VPU/GB).

### 5.3 The 64 MB `/dev/shm` Trap on Cloud Container Runtimes
- **Problem:** Container runtimes (Docker, containerd, Kubernetes, Arvados Crunch) allocate a default `/dev/shm` of exactly **64 MB**. PADV 1.10 natively requires loading `ksa_hg38.gbz` (~467 MB) into shared memory (`/dev/shm/GBZ_SHARED_MEMORY`), triggering an immediate `SIGBUS` (`Bus error`).
- **Mitigation Baked into CWL:** `pangenome-aware-deepvariant.cwl` overrides shared memory usage by passing `--make_examples_extra_args=use_loaded_gbz_shared_memory=false` and stubbing `load_gbz_into_shared_memory_command`. Graph data is loaded into the heap of each process shard via `GbzReader`.
- **Warning for Dung:** If you modify `pangenome-aware-deepvariant.cwl` or increase `num_shards` above 4, monitor heap RSS. At 4 shards, heap usage is ~41 GiB. Setting `num_shards: 8` would require `8 × 12.7 GiB ≈ 101.6 GiB` process heap, which would exceed a 128 GiB VM. **Keep `num_shards: 4` on 90–128 GiB instances.**

---

## 6. How to Run on the Cloud

Three distinct deployment patterns are available depending on Dung's infrastructure preferences.

```
Deployment Options Matrix:

Option A: Arvados on AWS / OCI
  + Native environment for Dung
  + Automatic step-level reuse (--enable-reuse)
  + Keepdocker image resolution
  - Requires Arvados cluster deployment on cloud VMs

Option B: Standalone cwltool / Toil on Single Large Cloud VM
  + Zero cluster management overhead; 100% standard Docker
  + Simplest path to an immediate green test run
  - Single VM vertical scaling only; no multi-node distribution

Option C: Managed Cloud Genomics Runner (AWS HealthOmics / OCI Compute)
  + Fully managed serverless execution; pay-per-job
  - Requires translating Keep inputs to S3/OCI URIs and mirroring all images to ECR/OCIR
```

### Option A: Arvados Cluster on Cloud (Recommended for Production)
This mirrors Dung's current operation on `cborg`:
1. **Submit using `arvados-cwl-runner`:**
   ```bash
   submit.sh dev main-vg.cwl main-vg.TS.e2e.json
   ```
2. **Step Reuse (`--enable-reuse`):** Enabled by default in `submit.sh`. When re-running a sample with modified phenotype terms or updated downstream ACMG criteria, Arvados instantly reuses the 14-hour Giraffe BAM and 17-hour DeepVariant VCF, re-running only the 30-minute annotation and interpretation tail.
3. **Workflow Registration:** Register the modernised CWL using the Arvados Python SDK (`workflows().update` or `workflows().create`), supplying `main-vg.cwl` packed with `cwltool --pack`.

### Option B: Standalone `cwltool` on a Single Cloud VM (Best for Staging)
Ideal for bringing up the pipeline quickly on an AWS `r6i.4xlarge` or OCI `VM.Standard.E5.Flex`:
1. Use the automation scripts from `local/` (`local/setup-venv.sh`, `local/stage-inputs.sh`, `local/pull-and-build-images.sh`).
2. Run via `local/run.sh`:
   ```bash
   bash local/run.sh /data/runs/wgs-TS main-vg.cwl local/main-vg.TS.standard.json
   ```
   `local/run.sh` launches `local/sample-resources.py` in the background, recording per-step CPU, RAM, and wall times into `step-resources.md`.

### Option C: Managed Cloud Genomics Service (AWS HealthOmics)
If deploying via AWS HealthOmics Workflows:
1. Push all 12 Docker images to AWS ECR.
2. Mirror all Keep data bundles to an S3 bucket (`s3://my-nelly-data/inputs/...`).
3. Replace `keep:<PDH>/<path>` in the input JSON with `s3://my-nelly-data/inputs/...`.
4. Register the workflow in HealthOmics using CWL v1.2 engine.

### 6.1 LLM Endpoint Requirement for the ACMG Agent

The ACMG agent (`acmg-agent.cwl`) reviews candidate variants through a bounded tool loop.

#### Measured Quality vs. Time Tradeoff

| Backend Option | Setup / Infrastructure | Wall Time (5 Variants) | Success Rate | Assessment Quality & Output Status |
|---|---|---:|---:|---|
| **vLLM Endpoint (`qwen3.8-27b`)** | `unimatrix01.kaust.edu.sa` or cloud GPU VM | **1,437 s (23.95 min)** | **100% (5/5 ok)** | `llm.status: "ok"`. Correctly downgraded RTEL1 BP2 with cited evidence; confirmed HEXA Pathogenic (10 pts). |
| **Local CPU GGUF (`Qwen2.5-3B`)** | Baked into `nelly-acmg-agent:v2` (`/model/model.gguf`) | 3,480 s (58.0 min) | 20% (1/5 ok, rest timeout) | `llm.status: "partial"`. Prompt processing on 8 CPU cores takes ~130s/turn; hits 50 min budget ceiling. |
| **Local CPU GGUF (`Qwen2.5-7B`)** | Built during WP4 prototyping | 2,037 s (timed out) | 0% (timeouts/rejections) | `llm.status: "partial"`. Exceeds 6-min per-variant budget on CPU. |
| **Deterministic Engine Only (`"off"`)** | Set `llm_backend: "off"` | **< 2 s** | **100% deterministic** | `llm.status: "disabled"`. 100% reproducible ClinGen SVI rule-based classification; no narrative LLM summary. |

#### Deployment Recommendations for Dung
1. **Option 1 (Fastest & Best): Cloud GPU Instance for vLLM.**
   - Deploy an AWS `g5.xlarge` (NVIDIA A10G 24GB, ~$1.00/hr) or OCI `VM.GPU.A10.1`.
   - Run vLLM container:
     ```bash
     docker run -d --gpus all -p 8000:8000 --ipc=host vllm/vllm-openai:latest \
       --model Qwen/Qwen2.5-14B-Instruct --max-model-len 16384
     ```
   - Set `llm_backend: "openai"`, `llm_endpoint: "http://<GPU_IP>:8000/v1"`, `llm_model: "Qwen/Qwen2.5-14B-Instruct"`.
2. **Option 2: VPN Tunnel to KAUST `unimatrix01`.**
   - Point to `http://unimatrix01.kaust.edu.sa:8000/v1` with model `qwen3.8-27b` over a site-to-site VPN.
3. **Option 3: Deterministic Fallback (`llm_backend: "off"`).**
   - For high-throughput cloud batch processing without GPU instances, set `llm_backend: "off"`. The deterministic engine provides 100% compliant ClinGen SVI classifications in 2 seconds per patient.

---

## 7. Validation Evidence & Known Gotchas

### 7.1 Summary of Validation Evidence

1. **Unit & Integration Test Suites:**
   - Pipeline directory test suite (`pytest -q`): **235 passed in 0.41s** (offline, no network/LLM/GPU).
   - Package test suite (`acmg-agent/`, `uv run pytest -q`): **134 passed in 0.15s**.
   - Input completeness static analysis (`tests/test_inputs_complete.py`, `tests/test_benchmark_inputs.py`): Passed across all production and benchmark JSON templates.
2. **CWL Validation:**
   - Every active and legacy CWL file passes `cwltool --validate`.
   - `cwltool --pack main-vg.cwl` completes clean with exit code 0.
3. **Phase 1 Local Smoke Runs (`cbontsr01`):**
   - Tay-Sachs tail (`local/test-tail.TS.json`): Exit code 0. `HEXA` ranked #2 in Exomiser, annotated with MANE Select splice donor consequence in VEP 116, classified **Pathogenic** (10 points, PVS1) by ACMG agent.
   - Huntington tail (`local/test-tail.HD.json`): Exit code 0. `exclude-eh-regions` cleanly filtered repeat loci; `HTT` CAG repeat expansion (18/46 repeats) identified as leading finding against STRchive thresholds.
4. **WP4b Container Live Execution:**
   - `nelly-acmg-agent:v2` executed on `ws` against `unimatrix01` in 1,437 s (23.95 min), reviewing all 5 candidates with `agent_status: ok`.
5. **Phase 2 Status:**
   - Alignment and variant calling jobs (`align-TS-standard`, Slurm ID `55405248`; `align-HD-standard`, Slurm ID `55405249`) currently executing on KAUST IBEX. Full WGS e2e wall time will be updated upon job completion.

### 7.2 Known Gotchas for Cloud Deployments (from `notes-cwl-gotchas.md`)

1. **Bash Expression Interpolation in CWL:** Both `$(...)` and `${...}` are evaluated by CWL parsers. In tool shell scripts, reference bash variables as `$VAR` (no braces), use backticks for subshells, and use `$(inputs.param)` only for genuine CWL inputs.
2. **No Backslash Line Continuations:** Single backslashes in multiline command strings are stripped by CWL string evaluators. Use trailing pipes (`|`) or single long lines.
3. **Double-Quoting Env Vars:** Never wrap an input string in quotes within `arguments` if it already contains quotes. Use `EnvVarRequirement` (`ACMG_LLM_API_KEY`) to safely pass tokens without shell parsing.
4. **No Container `ENTRYPOINT`:** Wrapper images (`nelly-exomiser:v4`, `nelly-acmg-agent:v2`) MUST have `ENTRYPOINT: null` and `CMD: ["bash"]` or `CMD: ["python"]`. Arvados prepends `ENTRYPOINT` to `baseCommand`, duplicating arguments.
5. **Optional File Inputs:** When passing optional files (e.g. `strchive_json`), never access `$(inputs.strchive_json.path)` directly. Use `valueFrom: "$(inputs.strchive_json ? inputs.strchive_json.path : null)"` to avoid null pointer exceptions.
6. **Exomiser Output Glob:** Exomiser 15.1.0 writes JSON output as `results/exomiser.jsonl` (JSON Lines), not `.json`. `exomiser.cwl` specifically globs `*.jsonl`.
7. **Read Group Format for Giraffe:** `vg giraffe --read-group` expects a bare ID string (`sample`), NOT a tab-separated `@RG\tID:...` line. Tab-containing strings crash downstream `samtools sort`.
8. **Reference Genome Consistency:** All steps downstream of Giraffe (`samtools`, `deepvariant`, `expansionhunter`, `bcftools-norm`, `exomiser`) MUST use the UCSC chr-named `hg38.fa` (`1d5ad2aefb6e809dddaffea3a0558415+5657/hg38.fa`). The numeric-chromosome `fastafile` is strictly reserved for VEP (`vep_fasta_file`).
9. **Phenopacket Proband ID:** Patient Genoor metadata must set proband ID equal to the VCF sample name (`sample`), or Exomiser silently skips variant scoring.
10. **VEP Secondary Files on Read-Only Mounts:** Under `cwltool`, `vep_fasta_file` must declare `secondaryFiles: [.fai, .gzi]`. If omitted, VEP attempts to create a `.lock` file in the read-only mount and fails.

---

## 8. Checklist for Dung (Step-by-Step to a First Green Run)

Follow this sequence to deploy and verify the workflow on AWS or OCI:

- [ ] **Step 1: Stage Data Bundles to Cloud Bucket**
  - Allocate an S3 bucket (e.g. `s3://nelly-cloud-data/`) or OCI Object Storage bucket.
  - From a KAUST host with Arvados access (`ws` or laptop with `.arvados-venv`), download the 8 core reference and database collections:
    - Pangenome graph & UCSC ref (`1d5ad2aefb6e809dddaffea3a0558415+5657`)
    - vg 1.78 indexes (`d704d713f7039bbc4026d55fafec7c33+2743`)
    - VEP 116 cache (`8afec8cdb796367144fa2bcec3a3dc59+624229`)
    - VEP numeric FASTA (`00528ccfa8b07fd08907077b1f4fe3e9+927`)
    - Exomiser 2512 bundle (`00cbcbd45419c2f051889639deff3431+39145`)
    - slivar & VEP plugins (`148bcb74e5d8268ac486d315e657f70b+3240`)
    - ACMG databases (`28a24f20fb9659e2ba4e4b9a17f7bec0+778`)
    - ACMG gene table & STRchive (`95d0b6027bd22592ac4ebde6412c0d63+100`)
  - Upload directories to the bucket matching the paths in `handover/data-bundles.tsv`.
  - For Arvados cloud deployments, run `arv-put` to upload into your cloud Keep instance.

- [ ] **Step 2: Mirror Docker Images**
  - If using AWS/OCI without Arvados Keep, create ECR/OCIR repositories for the 12 active images listed in `handover/container-images.tsv`.
  - Build `nelly-exomiser:v4` from `exomiser-wrapper/Dockerfile`.
  - Build `nelly-acmg-agent:v2` from `acmg-agent/Dockerfile`.
  - Re-tag and push all images to your cloud registry.

- [ ] **Step 3: Provision Compute Instance(s)**
  - Ensure compute nodes have at least **16 vCPUs (or 8–16 OCPUs)** and **128 GiB RAM** (e.g. AWS `r6i.4xlarge` or OCI `VM.Standard.E5.Flex` with 128 GB RAM).
  - Attach at least **1 TB EBS/Block Volume** scratch disk mounted at `/data/tmp` or similar (never `/tmp`).

- [ ] **Step 4: Configure LLM Backend**
  - If using the Unimatrix endpoint, verify network connectivity to `http://unimatrix01.kaust.edu.sa:8000/v1`.
  - If deploying in cloud VPC, launch a vLLM GPU container (`Qwen/Qwen2.5-14B-Instruct` or `32B`) or use an OpenAI-compatible API endpoint.
  - Test the endpoint:
    ```bash
    curl -H "Authorization: Bearer <API_KEY>" http://<LLM_HOST>:8000/v1/models
    ```
  - For zero-GPU batch operations, set `llm_backend: "off"`.

- [ ] **Step 5: Run Phase 1 Tail Smoke Test (Validation Gate 1)**
  - Before running a 14-hour alignment, run the lightweight tail workflow (`local/test-tail.cwl` with `local/test-tail.TS.json` or `local/test-tail.HD.json`).
  - Expected execution time: ~15–20 minutes.
  - **Success Criteria:**
    - Process exits 0.
    - `exomiser.genes.tsv` places `HEXA` at Rank #1 or #2.
    - `acmg-assessment.json` classifies `HEXA` as `PATHOGENIC` (10 points, PVS1 met).
    - Huntington run identifies `HTT` repeat expansion (46 repeats) as leading finding.

- [ ] **Step 6: Launch Full End-to-End WGS Run (Validation Gate 2)**
  - Submit `main-vg.cwl` with `main-vg.TS.e2e.json` (or `main-vg.HD.e2e.json`).
  - Inject the API key into the runtime JSON or environment (`ACMG_LLM_API_KEY`).
  - **Success Criteria:**
    - Full pipeline completes exit 0.
    - `aligned.cram` is generated and indexed.
    - `output_acmg_json` and `output_acmg_report` are generated.
    - In Tay-Sachs: `HEXA` leads the report as Pathogenic.
    - In Huntington: `HTT` repeat expansion leads the report.

- [ ] **Step 7: Execute GIAB Benchmark (Validation Gate 3)**
  - Run `benchmark/benchmark-giab.cwl` across both caller modes:
    - `benchmark/benchmark-giab.HG003.standard.json`
    - `benchmark/benchmark-giab.HG003.pangenome_aware.json`
  - Collect results using `benchmark/collect.py`:
    ```bash
    python3 benchmark/collect.py --runs runs.json --out benchmark/results-2026-11.md
    ```
  - Confirm SNP F1 $\ge$ 0.9970 on HG003 and select the winning caller for default `dv_mode`.
