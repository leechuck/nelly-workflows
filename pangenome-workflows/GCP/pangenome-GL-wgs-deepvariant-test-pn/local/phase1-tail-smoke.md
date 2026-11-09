# Phase 1 — Tail-Only Smoke Run Log (`local/test-tail.cwl`)

**Host:** `cbontsr01.kaust.edu.sa` (`10.67.24.207`, AMD EPYC, 256 threads, 1 TB RAM; Ubuntu 20.04 -> upgraded in-place to Ubuntu 24.04.5 LTS, kernel `5.15.0-203-generic`)  
**Engine:** `cwltool 3.3.20261106131234` (`--parallel --timestamps`)  
**Date:** 2026-11-09

## 1. Phase 1a Runs (Pre-WP4 Merge: `phenofrommetadata` + `exclude-eh-regions` + `slivar` -> `vep` -> `exomiser`)

1. **Tay-Sachs (`tail-TS`)**:
   - Command:
     ```bash
     bash local/run.sh /data/nelly-local/runs/tail-TS local/test-tail.cwl local/test-tail.TS.json
     ```
   - Inputs:
     - `normalized.vcf.gz`(+.tbi), `eh_repeats.vcf.gz`(+.tbi), `eh_repeats.json`: `keep:b2c8c2b940de5bf4d7881b652e12c24f+18026`
     - `metadata.yaml`: `keep:5c88241c12f7dc7553ffb3f2ac1ba54b+39418/metadata.yaml` (`HP:0002376`, `HP:0001252`, `HP:0010780`; sex `MALE`)
   - Exit status: **`0` (`Final process status is success`)**

2. **Huntington (`tail-HD`)**:
   - Command:
     ```bash
     bash local/run.sh /data/nelly-local/runs/tail-HD local/test-tail.cwl local/test-tail.HD.json
     ```
   - Inputs:
     - `normalized.vcf.gz`(+.tbi), `eh_repeats.vcf.gz`(+.tbi), `eh_repeats.json`: `keep:6d04e84a0ab500c5036ec6713756783b+17943`
     - `metadata.yaml`: `keep:d03401bdd2443e6340903acb552fd0c3+39418/metadata.yaml` (`HP:0000726`, `HP:0002072`; sex `MALE`)
   - Exit status: **`0` (`Final process status is success`)**

## 2. Phase 1b Runs (Full Tail with Merged ACMG Agent Step + Unimatrix Qwen3.8-27B)

Following merge `87fbe22` (WP4/WP4b), `local/test-tail.cwl` was extended with `acmg-extract-candidates -> vep-candidates -> acmg-agent` (`nelly-acmg-agent:v2`, image ID `8ed50b5e81a5`, loaded from `keep:c15a4a032ac6d0bcb2c10e8bd9653a05+1731`), wired to the Unimatrix vLLM endpoint (`http://unimatrix01.kaust.edu.sa:8000/v1`, model `qwen3.8-27b`, `llm_max_tokens: 4096`), `clinvar.vcf.gz`(+.tbi) from `keep:28a24f20fb9659e2ba4e4b9a17f7bec0+778`, and `gene_table.tsv` + `STRchive-loci.json` from `keep:95d0b6027bd22592ac4ebde6412c0d63+100`.

### 2.1 Fail-Closed Network-Isolated Run (`tail-TS-acmg-netnone` & `tail-HD-acmg-netnone`)
- On the first execution, because `acmg-agent.cwl` did not declare `NetworkAccess: { networkAccess: true }`, `cwltool` passed `--net=none` to `docker run`.
- Both runs exited **`0`** (`Final process status is success`), empirically proving the WP4 fail-closed fallback contract: after 3 connection retries (~15.8 s wall time, 71.3 MiB RSS), `acmg-agent` fell back cleanly to the deterministic engine (`llm.status: "partial"`), producing complete `acmg-assessment.json` and `acmg-report.md` outputs with **`HEXA` `Likely pathogenic` (8 pts)** leading for TS and **`HTT` `Pathogenic` repeat expansion (`18/46`)** leading for HD.
- Fix applied (`commit 5a9e992`): Added `NetworkAccess: { networkAccess: true }` to `local/test-tail.cwl` and defaulted `JobBase.networkaccess = True` in `local/run.sh` (which also prevents `cwltool --singularity` on IBEX from passing root-only `--net --network none`).

### 2.2 Live Unimatrix Qwen3.8-27B Run (`tail-TS-acmg` & `tail-HD-acmg`)
1. **Tay-Sachs (`tail-TS-acmg`)**:
   - Output dir: `/data/nelly-local/runs/tail-TS-acmg/` (log: `/data/nelly-local/logs/tail-TS-acmg.log`)
   - Command:
     ```bash
     bash local/run.sh /data/nelly-local/runs/tail-TS-acmg local/test-tail.cwl local/test-tail.TS.json
     ```
   - Exit status: **`0` (`[2026-11-09 21:14:44] INFO Final process status is success`)**
2. **Huntington (`tail-HD-acmg`)**:
   - Output dir: `/data/nelly-local/runs/tail-HD-acmg/` (log: `/data/nelly-local/logs/tail-HD-acmg.log`)
   - Command:
     ```bash
     bash local/run.sh /data/nelly-local/runs/tail-HD-acmg local/test-tail.cwl local/test-tail.HD.json
     ```
   - Exit status: **`0` (`[2026-11-09 21:16:26] INFO Final process status is success`)**

## 3. Per-Step Resource Usage Across Full Tail (`step-resources.md`)

| Step | Sample | Container ID | Image | Wall Time (s) | Peak Process RSS (MiB) | Peak Cgroup Mem incl. Page Cache (MiB) |
|---|---|---|---|---:|---:|---:|
| `phenofrommetadata` | TS | `ffu6194y` | `dunglam2k/phenofrommetadata:v1.01` | 1.0 | <15.0 | <50.0 |
| `phenofrommetadata` | HD | `zwpa53td` | `dunglam2k/phenofrommetadata:v1.01` | 1.0 | <15.0 | <50.0 |
| `exclude-eh-regions` | TS | `3a71613052ab` | `staphb/bcftools:1.22` | 8.7 | 15.7 | 105.7 |
| `exclude-eh-regions` | HD | `44c5c5392817` | `staphb/bcftools:1.22` | 8.7 | 12.7 | 102.8 |
| `exomiser` | TS | `459ae8d0335c` | `nelly-exomiser:v4` (15.1.0 / 2512) | 69.4 | 4509.2 | 4548.4 |
| `exomiser` | HD | `228d470e1d9d` | `nelly-exomiser:v4` (15.1.0 / 2512) | 68.1 | 4040.5 | 4108.1 |
| `slivar` | TS | `a8b0347c44cc` | `brentp/slivar:v0.3.0` | 166.6 | 2534.4 | 2542.6 |
| `slivar` | HD | `f244d6c17dc4` | `brentp/slivar:v0.3.0` | 165.4 | 2534.4 | 2542.7 |
| `vep` (full VCF) | TS | `74180e74bc17` | `ensemblorg/ensembl-vep:release_116.3` | 703.3 | 837.4 | 852.3 |
| `vep` (full VCF) | HD | `192d5acf7b59` | `ensemblorg/ensembl-vep:release_116.3` | 704.2 | 860.3 | 875.4 |
| `acmg-extract-candidates` | TS | `c2aabc347cb3` | `nelly-acmg-agent:v2` | 1.0 | 23.8 | 23.8 |
| `acmg-extract-candidates` | HD | `074e626f386c` | `nelly-acmg-agent:v2` | 1.0 | 23.8 | 23.8 |
| `vep-candidates` (mini-VCF) | TS | `e76e770f93c8` | `ensemblorg/ensembl-vep:release_116.3` | 1.3 (4.0 cwltool) | 254.3 | 257.6 |
| `vep-candidates` (mini-VCF) | HD | `32acb81713a8` | `ensemblorg/ensembl-vep:release_116.3` | 1.3 (4.0 cwltool) | 215.5 | 219.2 |
| `acmg-agent` (live Qwen3.8-27B) | TS | `f18e8c1970d7` | `nelly-acmg-agent:v2` | **1603.1** | **71.2** | **73.2** |
| `acmg-agent` (live Qwen3.8-27B) | HD | `9192a0344c38` | `nelly-acmg-agent:v2` | **1703.2** | **71.1** | **73.2** |
| `acmg-agent` (`--net=none` fallback) | TS / HD | `e8ad4f8aec7c` / `4124d1b694c8` | `nelly-acmg-agent:v2` | 15.8 / 15.9 | 71.3 / 71.1 | 73.3 / 73.1 |

## 4. Output & Grounding Verification

### 4.1 Tay-Sachs (`tail-TS-acmg`)
- **Exomiser 15.1.0 + 2512 ranking (`exomiser.genes.tsv`)**:
  - Rank #1: `SEMA6B` (`AD`, combined `0.7898`, pheno `0.6717`, variant `1.0000`)
  - **Rank #2: `HEXA` (`AR`, combined `0.5359`, pheno `0.5450`, variant `1.0000`, `WHITELIST_VARIANT=1`)** — meets top-2 gate.
- **`acmg-assessment.json` & `acmg-report.md`**:
  - **Leading candidate**:
    ```json
    {
      "kind": "snv",
      "gene": "HEXA",
      "variant": "chr15:72346234 C>G",
      "hgvs": "HEXA:ENST00000268097.10:c.1421+1G>C:p.?",
      "zygosity": "hom",
      "moi": "AR",
      "disease": "GM2-gangliosidosis, several forms (OMIM:272800)",
      "classification": "Likely pathogenic",
      "points": 8,
      "rule_based_class": "VUS",
      "conflict": true
    }
    ```
  - **LLM status**: `{"backend": "openai", "model": "qwen3.8-27b", "status": "partial"}` — 4 of the 5 selected review candidates (`#1 HEXA`, `#2 SEMA6B`, `#3 PAK2`, `#5 RTEL1` — **including the leading candidate `HEXA`**) completed with `agent_status == "ok"` and full grounded LLM `agent_rationale` and `agent_narrative`; `#4 PKD1` hit the per-variant tool-call/validation budget and fell back cleanly to its deterministic classification (`Likely pathogenic`, 8 pts), resulting in run-level `status: "partial"`.
  - **Grounding check**: Programmatic token intersection against all HGNC symbols in `gene_table.tsv` (19,212 genes) confirmed **0 out-of-pack genes** in `acmg-report.md` and all `agent_rationale` / `agent_narrative` strings.

### 4.2 Huntington (`tail-HD-acmg`)
- **`acmg-assessment.json` & `acmg-report.md`**:
  - **Leading candidate**:
    ```json
    {
      "kind": "expansion",
      "gene": "HTT",
      "variant": "HTT repeat (18/46)",
      "hgvs": null,
      "zygosity": "het",
      "moi": "AD",
      "disease": "Huntington disease",
      "classification": "Pathogenic",
      "points": null,
      "rule_based_class": "Pathogenic",
      "conflict": false,
      "allele1": 18,
      "allele2": 46,
      "genotypes": "18/46",
      "ci": "18-18/46-58",
      "normal_max": 26,
      "pathogenic_min": 36,
      "status": "pathogenic_range"
    }
    ```
  - **LLM status**: **`{"backend": "openai", "model": "qwen3.8-27b", "status": "ok"}`** — **all 5 selected review candidates (`#1 SEMA6B`, `#2 ABCA7`, `#3 PKD1`, `#4 PAK2`, `#5 RTEL1`) succeeded with `agent_status == "ok"`!**
  - **Validated LLM adjustments**: On `PKD1` and `RTEL1`, the LLM agent legitimately adjusted `BP2` from `met (supporting)` to `not_met` (`adjusted_by_agent: ["BP2"]`, changing points from `8` to `9`, `Likely pathogenic`) with a contract-validated reason citing the `moi` (`AD`) and `other_var_in_gene` evidence fields (`single-sample data cannot establish in-trans configuration`).
  - **Grounding check**: Programmatic token intersection against `gene_table.tsv` confirmed **0 out-of-pack genes** in `acmg-report.md` and all `agent_rationale` / `agent_narrative` strings.
