# WP4 → supervisor: intended `main-vg.cwl` / input-JSON wiring

WP1 was unfinished at kickoff, so WP4 did **not** touch `main-vg.cwl` or the
`main-vg*.json` files. Everything in this file is ready to paste; the two tool
CWLs (`acmg-extract-candidates.cwl`, `acmg-agent.cwl`) are committed and pass
`cwltool --validate` — apply the blocks below in one pass and re-validate.

## 1. `main-vg.cwl` — workflow inputs

```yaml
  # ACMG agent (will replace genome-linter): independently refreshable Keep PDHs
  acmg_clinvar_vcf:
    type: File
    secondaryFiles:
      - .tbi
  acmg_gene_table:
    type: File
  acmg_top_candidates:
    type: int
    default: 10
  # WK ... optional future input (WP3 Task 3.4 repeats): acmg_strchive_json: File?
```

The `genomelinter_output_name` / `genomelinter_model` / `openrouter_api_key` /
`openrouter_*` inputs (LLM-ranker era) are removed together with the step.

## 2. `main-vg.cwl` — step replacements

`exomiser` keeps its position behind `exclude-eh-regions`. Then:

```yaml
  # NEW (replaces step `genomelinter`): candidate extraction -> mini-VEP -> agent
  acmg-extract-candidates:
    in:
      exomiser_genes_tsv: exomiser/exomiser_genes_tsv
      exomiser_variants_tsv: exomiser/exomiser_variants_tsv
      top_candidates: acmg_top_candidates
      # optional: reference_index_fai (fai from the workflow reference) gives
      # ##contig lines + chr naming; leave null to let extract use chr-prefix
      reference_fai: reference_fai
    out: [candidates_vcf]
    run: acmg-extract-candidates.cwl

  # Second use of the SAME vep.cwl on the mini-VCF (design 8.1). All inputs are
  # the same as the existing `vep` step except vep_input and vep_output_file.
  # When WP5 lands, vep.cwl gains --af_gnomade/--af_gnomadg/--max_af, which the
  # ACMG builder reads from this output.
  vep-candidates:
    in:
      vep_input: acmg-extract-candidates/candidates_vcf
      vep_assembly: vep_assembly
      vep_output_file: vep_candidates_output_file   # see workflow inputs below
      pavs_custom_file: pavs_custom_file
      pavs_custom_args: pavs_custom_args
      go_custom_file: go_custom_file
      go_custom_args: go_custom_args
      hpo_custom_file: hpo_custom_file
      hpo_custom_args: hpo_custom_args
      ppi_custom_file: ppi_custom_file
      ppi_custom_args: ppi_custom_args
      vep_dir: vep_dir
      vep_fasta_file: vep_fasta_file
    out: [vep_console_out, vep_output]
    run: vep.cwl

  # ACMG evaluation agent. The original `vep` step and workflow output
  # `output_after_vep` stay unchanged.
  acmg-agent:
    in:
      vep_vcf: vep-candidates/vep_output
      clinvar_vcf: acmg_clinvar_vcf
      gene_table: acmg_gene_table
      exomiser_genes_tsv: exomiser/exomiser_genes_tsv
      exomiser_variants_tsv: exomiser/exomiser_variants_tsv
      eh_json: expansionhunter/json
      slivar_vcf: slivar/slivar_output
      # phenofrommetadata outputs (gotcha 11: sample must match the VCF column)
      sample_name: sample
      hpo_terms: phenofrommetadata/hpo_terms            # File -> see note
      phenotype: phenofrommetadata/pheno_output
      sample_sex: phenofrommetadata/exomiser_sex
      top_candidates: acmg_top_candidates
    out: [acmg_json, acmg_report]
    run: acmg-agent.cwl
```

Notes for the supervisor:

* **`hpo_terms` type mismatch.** `phenofrommetadata/hpo_terms` is a `File`
  (`hpo_terms.txt`) but `acmg-agent.cwl` declares `hpo_terms: string`. Add a
  tiny sink in the step with valueFrom reading the file contents, or extend
  `phenofrommetadata.cwl` with a term-list output; either is one line. (The
  agent currently uses HPO only for the JSON header.)
* **`vep_candidates_output_file`**: declare as workflow input
  `type: string, default: "candidates.vep.vcf.gz"` (mirrors the existing
  `vep_output_file` input).
* **`expansionhunter/json`** is `eh_repeats.json` from the `expansionhunter`
  step — already referenced by `genomelinter`, same binding.
* **`sample`** is the DeepVariant sample-name input already used by Exomiser
  (gotcha 11).
* **`reference_fai`** does not exist yet in `main-vg.cwl`; if absent, pass
  nothing (the extract tool defaults to chr-prefixed contig names, which the
  UCSC hg38 ref lineage expects).

## 3. `main-vg.cwl` — workflow outputs

```yaml
  # REPLACES output_after_genomelinter (remove that output)
  output_acmg_json:
    type: File
    outputSource: acmg-agent/acmg_json
  output_acmg_report:
    type: File
    outputSource: acmg-agent/acmg_report
```

## 4. Input JSONs (`main-vg-pn.json`, `main-vg.HD.e2e.json`, `main-vg.TS.e2e.json`)

Remove `genomelinter_output_name` / `genomelinter_model` / `openrouter_*` and
add:

```json
  "acmg_clinvar_vcf": {"class": "File", "location": "keep:TODO-ACMG-PDH/clinvar.vcf.gz",
                       "secondaryFiles": [{"class": "File", "location": "keep:TODO-ACMG-PDH/clinvar.vcf.gz.tbi"}]},
  "acmg_gene_table": {"class": "File", "location": "keep:TODO-ACMG-PDH/gene_table.tsv"},
  "acmg_top_candidates": 10
```

Keep PDHs are unknown to WP4 — the supervisor stages them with
`benchmark/stage_clinvar.sh` (ClinVar GRCh38 VCF+tbi) and
`benchmark/build_gene_table.py --out gene_table.tsv`, then `arv-put`.
`tests/test_inputs_complete.py` stays green once `main-vg.cwl` declares the
three new inputs and the JSONs supply them (extend it if it asserts removals).

## 5. Decisions recorded at wiring time

* `acmg-extract-candidates` defaults to `--top 15` while `acmg_top_candidates`
  (assess) defaults to 10: VEP on a handful of extra records is cheap and
  rank ties outside the top-10 (same gene, second contributing variant) keep
  their chance to be evaluated.
* The merged `acmg_agent.cli extract` writes an **uncompressed** candidates
  VCF, so the tool WIP outputs/candidates to `candidates.vcf` (`vep.cwl` takes
  plain VCF via `-input_file` and its output stays bgzipped). If VEP 116
  (WP5) requires a `.gz` literal, `bcftools view -Oz` can be inserted later.
* WP3's `cli assess` has no `--exomiser-jsonl` flag yet (design §4.3 lists
  it); `acmg-agent.cwl` therefore omits it. The `SerialExomiserJson`
  Exomiser output stays available for WP3 follow-ups.
* Workflow-level `hpo_terms`/`phenotype`/`sex` travel: HPO + sex as plain
  argv, phenotype via `EnvVarRequirement ACMG_PHENOTYPE` (the CLI reads it as
  the `--phenotype` default), so labels with semicolons/quotes can't break
  the shell (gotchas #1/#3).
* Genomelinter removal: delete step `genomelinter`, workflow outputs
  `output_after_genomelinter`, inputs `genomelinter_*`, `openrouter_*`; the
  grounded-narrative role moves into the agent (status: narrative generated
  by the LLM from the final criteria list, classification decided in code).
