# Phase 1 — Tail-Only Smoke Run Log (`local/test-tail.cwl`)

**Host:** `cbontsr01.kaust.edu.sa` (AMD EPYC, 256 threads, 1 TB RAM, Ubuntu 20.04, kernel `5.15.0-203-generic`)  
**Engine:** `cwltool 3.3.20261106131234` (`--parallel --timestamps`)  
**Date:** 2026-11-09

## 1. Runs Executed

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

## 2. Per-Step Resource Usage (`step-resources.md`)

| Step | Sample | Container ID | Image | Wall Time (s) | Peak Process RSS (MiB) | Peak Cgroup Mem incl. Page Cache (MiB) |
|---|---|---|---|---:|---:|---:|
| `phenofrommetadata` | TS | `ffu6194y` | `dunglam2k/phenofrommetadata:v1.01` | 1.1 | <15.0 | <50.0 |
| `phenofrommetadata` | HD | `zwpa53td` | `dunglam2k/phenofrommetadata:v1.01` | 1.1 | <15.0 | <50.0 |
| `exclude-eh-regions` | TS | `24074de83c6f` | `staphb/bcftools:1.22` | 9.4 | 12.4 | 199.7 |
| `exclude-eh-regions` | HD | `dde85b121548` | `staphb/bcftools:1.22` | 8.4 | 12.7 | 204.8 |
| `exomiser` | TS | `4b657efea3fd` | `nelly-exomiser:v4` (15.1.0 / 2512) | 76.8 | 3545.4 | 4204.4 |
| `exomiser` | HD | `aa8e0da7a9fd` | `nelly-exomiser:v4` (15.1.0 / 2512) | 75.7 | 3584.9 | 4491.0 |
| `slivar` | TS | `80eec8ed061f` | `brentp/slivar:v0.3.0` | 167.9 | 2534.4 | 4865.0 |
| `slivar` | HD | `6b1d7f88bb64` | `brentp/slivar:v0.3.0` | 167.9 | 2534.4 | 4131.8 |
| `vep` | TS | `0adfb7322055` | `ensemblorg/ensembl-vep:release_116.3` | 719.1 | 837.5 | 4944.7 |
| `vep` | HD | `622dd7cf8a9f` | `ensemblorg/ensembl-vep:release_116.3` | 719.1 | 860.6 | 12158.6 |

*(Note: When `vep.cwl` ran solo without a concurrent second `vep` instance, wall time was 605 s and peak process RSS was ~850 MiB.)*

## 3. Output Verification

### 3.1 Tay-Sachs (`tail-TS`)
- **Exomiser 15.1.0 + 2512 ranking (`exomiser.genes.tsv`)**:
  - Rank #1: `SEMA6B` (`AD`, combined `0.7898`, pheno `0.6717`, variant `1.0000`)
  - **Rank #2: `HEXA` (`AR`, combined `0.5359`, pheno `0.5450`, variant `1.0000`, `WHITELIST_VARIANT=1`)** — meets the top-2 gate.
  - Contributing variant in `exomiser.variants.tsv`: `15-72346234-C-G_AR` (`rs147324677`, homozygous `1/1`).
- **VEP 116.3 output (`vep-output.vcf.gz`)**:
  - `CSQ` format verified: `Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|ENSP|SOURCE|SIFT|PolyPhen|HGVS_OFFSET|AF|gnomADe_AF|...|gnomADg_AF|...|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PAVS|GO_CLASSES|PHENOTYPE|PPI`
  - `chr15:72346234 C>G` (`HEXA`): annotated as `splice_donor_variant|HIGH|HEXA|ENSG00000213614|Transcript|ENST00000268097.10|protein_coding||12/13|ENST00000268097.10:c.1421+1G>C|...|MANE_Select|NM_000520.6`.
- **Deterministic `acmg-agent` check (`acmg-report.det.md`)**:
  - **Leading candidate**: `HEXA` `chr15:72346234 C>G` (`HEXA:ENST00000268097.10:c.1421+1G>C:p.?`, `hom`, `AR`, `GM2-gangliosidosis, several forms (OMIM:272800)`), **Likely pathogenic** (8 pts, `PVS1_VeryStrong`).

### 3.2 Huntington (`tail-HD`)
- **`exclude-eh-regions`**:
  - Excluded ExpansionHunter repeat regions (pad=20 bp) cleanly into `exomiser_input.vcf.gz` (92,014,078 bytes).
- **ExpansionHunter repeat call (`eh_repeats.vcf.gz` / `eh_repeats.json`)**:
  - `chr4:3074876` `REPID=HTT`, `RU=CAG`, `REF=19`, `GT=1/2`, `REPCN=18/46`, `REPCI=18-18/46-58` (pathogenic threshold >= 36 repeats).
- **Deterministic `acmg-agent` check (`acmg-report.det.md`)**:
  - **Leading finding**: pathogenic-range repeat expansion in **`HTT`** (alleles `18/46`, 95% CI `18-18/46-58`; pathogenic threshold >= 36 repeats; MOI `AD`) associated with **Huntington disease**, with STRchive loci thresholds active (`strchive_used: true`).
