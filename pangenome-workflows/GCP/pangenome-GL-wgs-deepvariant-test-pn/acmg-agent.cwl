cwlVersion: v1.2
class: CommandLineTool
# ACMG evaluation agent (replaces genome-linter-llm.cwl). Deterministic
# ACMG/AMP 2015 + ClinGen SVI engine plus a bounded LLM review pass over the
# evidence pack (max 12 tool calls and 6 min per variant, 50 min global,
# tier 0-2 + top-5 candidates, at most 8 variants). All contract rules are
# enforced in code inside the container (acmg_agent.validate_submission);
# a rejected/failed/timed-out LLM assessment falls back to the deterministic
# result, so the step always produces a complete report.
#
# LLM backend is configured via workflow inputs mapped to environment variables
# (openai endpoint by default, with local 3B GGUF via llama.cpp as offline
# fallback inside the image). The API key travels only via EnvVarRequirement
# and never appears in argv.
baseCommand: [python, -m, acmg_agent.cli, assess]
$namespaces:
  cwltool: "http://commonwl.org/cwltool#"
hints:
  # Keeps the LLM key out of stored inputs and logs (cwltool masks it; Arvados
  # stores it in secret_mounts instead of the container request's visible inputs).
  cwltool:Secrets:
    secrets: [llm_api_key]
requirements:
  # The LLM reviewer calls an OpenAI-compatible endpoint; without this, Docker-based
  # runners (Arvados, cwltool --docker) start the container with --net=none and the
  # step silently falls back to the deterministic result (llm.status failed).
  NetworkAccess:
    networkAccess: true
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 8
    ramMin: 20480
  EnvVarRequirement:
    envDef:
      PYTHONPATH: /app
      # llama.cpp thread count from the runtime cores (design 4.3 LLM budget)
      ACMG_THREADS: $(String(runtime.cores))
      # Phenotype labels travel via env, never through argv shell-quoting
      # (notes-cwl-gotchas.md #1/#3); read by acmg_agent.cli --phenotype default.
      ACMG_PHENOTYPE: $(inputs.phenotype)
      # LLM backend configuration via env vars (never shell-reparsed; API key
      # never appears on the command line).
      ACMG_LLM_BACKEND: $(inputs.llm_backend)
      ACMG_LLM_ENDPOINT: $(inputs.llm_endpoint)
      ACMG_LLM_MODEL: $(inputs.llm_model)
      ACMG_LLM_API_KEY: $(inputs.llm_api_key)
      ACMG_LLM_MAX_TOKENS: $(String(inputs.llm_max_tokens))
      ACMG_PIPELINE_VERSION: $(inputs.pipeline_version)
  DockerRequirement:
    dockerPull: 'nelly-acmg-agent:v8'
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
  - prefix: --patient-id
    valueFrom: $(inputs.patient_id)
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
  # STRchive repeat locus definitions JSON (optional, guarded per gotcha 5).
  - prefix: --strchive-json
    valueFrom: "$(inputs.strchive_json ? inputs.strchive_json.path : null)"
  # slivar-annotated normalized VCF (optional secondary AF/nhomalt source).
  - prefix: --slivar-vcf
    valueFrom: "$(inputs.slivar_vcf ? inputs.slivar_vcf.path : null)"
  - prefix: --output-json
    valueFrom: acmg-assessment.json
  - prefix: --output-md
    valueFrom: acmg-report.md
  - prefix: --output-pdf
    valueFrom: acmg-report.pdf
inputs:
  vep_vcf:
    type: File
  clinvar_vcf:
    type: File
    secondaryFiles:
      - .tbi
  gene_table:
    type: File
  strchive_json:
    type: File?
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
    default: "sample"
  patient_id:
    type: string
    default: ""
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
  llm_backend:
    type: string
    default: "openai"
  llm_endpoint:
    type: string
    default: "http://unimatrix01.kaust.edu.sa:8000/v1"
  llm_model:
    type: string
    default: "qwen3.8-27b"
  llm_api_key:
    type: string
    default: ""
  llm_max_tokens:
    type: int
    default: 4096
  pipeline_version:
    type: string
    default: "nelly-workflows modernise-2026-11"
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
  acmg_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"   # PDF
    outputBinding:
      glob: acmg-report.pdf
  llm_trace:
    type: File?
    outputBinding:
      glob: llm-trace.jsonl
