cwlVersion: v1.1
class: CommandLineTool
# ACMG evaluation agent (replaces genome-linter-llm.cwl). Deterministic
# ACMG/AMP 2015 + ClinGen SVI engine plus a bounded LLM review pass over the
# evidence pack (max 12 tool calls and 6 min per variant, 50 min global,
# tier 0-2 + top-5 candidates, at most 8 variants). All contract rules are
# enforced in code inside the container (acmg_agent.validate_submission);
# a rejected/failed/timed-out LLM assessment falls back to the deterministic
# result, so the step always produces a complete report.
#
# The LLM is local (Qwen2.5-7B-Instruct Q4_K_M baked into the image, CPU only,
# no network): ACMG_THREADS must receive $(runtime.cores) so llama.cpp sizes
# its thread pool from the assigned container cores.
#
# NOTE (WP4 supervisor integration): while WP1 is unfinished this tool is
# standalone; the intended main-vg.cwl step/inputs/outputs blocks live in
# acmg-agent/WIRING.md (including the removal of the genomelinter step).
baseCommand: [python, -m, acmg_agent.cli, assess]
requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 8
    ramMin: 20480
  EnvVarRequirement:
    envDef:
      PYTHONPATH: /app
      # llama.cpp thread count from the runtime cores (design 4.3 LLM budget)
      ACMG_THREADS: $(runtime.cores)
      # Phenotype labels travel via env, never through argv shell-quoting
      # (notes-cwl-gotchas.md #1/#3); read by acmg_agent.cli --phenotype default.
      ACMG_PHENOTYPE: $(inputs.phenotype)
  DockerRequirement:
    dockerPull: 'nelly-acmg-agent:v1'
arguments:
  # VEP-annotated mini-VCF from the vep-candidates step: the ACMG evidence
  # builder's primary annotation source (consequence, MANE/HGVS, AFs).
  - prefix: --vep-vcf
    valueFrom: $(inputs.vep_vcf.path)
  - prefix: --clinvar-vcf
    valueFrom: $(inputs.clinvar_vcf.path)
  - prefix: --gene-table
    valueFrom: $(inputs.gene_table.path)
  - prefix: --sample
    valueFrom: $(inputs.sample_name)
  - prefix: --top
    valueFrom: $(inputs.top_candidates)
  # HPO term IDs and sex as plain arguments; free-text labels via env only.
  - prefix: --hpo-terms
    valueFrom: $(inputs.hpo_terms)
  - prefix: --sex
    valueFrom: $(inputs.sample_sex)
  # Exomiser outputs (guarded per gotcha 5 -- optional inputs must not be
  # dereferenced when null). Always wired in main-vg, but tolerable as null.
  - prefix: --genes-tsv
    valueFrom: "$(inputs.exomiser_genes_tsv ? inputs.exomiser_genes_tsv.path : null)"
  - prefix: --variants-tsv
    valueFrom: "$(inputs.exomiser_variants_tsv ? inputs.exomiser_variants_tsv.path : null)"
  # ExpansionHunter repeats JSON (optional): pathogenic-range expansions
  # surface as tier-0 candidates; omitted when not wired.
  - prefix: --eh-json
    valueFrom: "$(inputs.eh_json ? inputs.eh_json.path : null)"
  # slivar-annotated normalized VCF (optional secondary AF/nhomalt source).
  - prefix: --slivar-vcf
    valueFrom: "$(inputs.slivar_vcf ? inputs.slivar_vcf.path : null)"
  - prefix: --output-json
    valueFrom: acmg-assessment.json
  - prefix: --output-md
    valueFrom: acmg-report.md
inputs:
  vep_vcf:
    type: File
  clinvar_vcf:
    type: File
    secondaryFiles:
      - .tbi
  gene_table:
    type: File
  exomiser_genes_tsv:
    type: File?
  exomiser_variants_tsv:
    type: File?
  eh_json:
    type: File?
  slivar_vcf:
    type: File?
  sample_name:
    type: string
    default: "proband"
  hpo_terms:
    type: string
    default: ""
  # Free-text phenotype labels from phenofrommetadata (read via ACMG_PHENOTYPE)
  phenotype:
    type: string
    default: ""
  sample_sex:
    type: string
    default: ""
  top_candidates:
    type: int
    default: 10
outputs:
  acmg_json:
    type: File
    format: "http://edamontology.org/format_3464"   # JSON
    outputBinding:
      glob: acmg-assessment.json
  acmg_report:
    type: File
    format: "http://edamontology.org/format_1960"   # Markdown
    outputBinding:
      glob: acmg-report.md
