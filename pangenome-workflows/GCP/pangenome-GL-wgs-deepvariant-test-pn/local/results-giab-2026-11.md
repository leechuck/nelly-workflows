# GIAB Benchmark Evaluation Results

## Performance Summary (Confident Regions, Filter = PASS)

| sample | graph | caller | SNP recall | SNP precision | SNP F1 | INDEL recall | INDEL precision | INDEL F1 |
|---|---|---|---|---|---|---|---|---|
| HG003 | Saudi (KSA) | vg 1.78 + DV 1.10 (standard) | - | - | - | - | - | - |
| HG003 | Saudi (KSA) | vg 1.78 + PADV 1.10 (pangenome_aware) | - | - | - | - | - | - |
| HG003 | Saudi (KSA) | legacy (vg 1.54 + DV 1.6) | - | - | - | - | - | - |
| HG002 | Saudi (KSA) | vg 1.78 + DV 1.10 (standard) | - | - | - | - | - | - |
| HG002 | Saudi (KSA) | vg 1.78 + PADV 1.10 (pangenome_aware) | - | - | - | - | - | - |

## Published Reference Baselines (Yang et al. — vg 1.54 + DV 1.6, HG003 vs GIAB v4.2.1)

| Reference Study | Graph | Caller / Aligner | SNP F1 | SNP Errors | INDEL Errors | Difficult Regions SNP Recall |
|---|---|---|---|---|---|---|
| HG003 (Yang et al.) | GRCh38 | BWA-MEM2 + DV 1.6 | 0.9956 | 29,052 | 7,064 | 0.964 |
| HG003 (Yang et al.) | HPRC v1.0 | vg 1.54 + DV 1.6 | 0.9973 | 18,099 | 6,609 | 0.984 |
| HG003 (Yang et al.) | Saudi (KSA) | vg 1.54 + DV 1.6 | 0.9971 | 18,996 | 6,582 | 0.983 |
