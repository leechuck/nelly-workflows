cwlVersion: v1.1
class: CommandLineTool
# Extract the ACMG candidate set from Exomiser's ranked output into a tiny VCF
# (the contributing variants of the top-N re-ranked genes) so a second, cheap
# VEP invocation can annotate exactly those records before the ACMG agent runs.
# Only a few dozen records reach VEP, so the mini-VEP stays well inside the
# 60-min step budget (design section 8.1 'acmg-extract-candidates').
#
# NOTE (WP4 supervisor integration): while WP1 is not finished this tool is
# standalone; the intended main-vg.cwl step block lives in acmg-agent/WIRING.md.
baseCommand: [python, -m, acmg_agent.cli, extract]
requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 1
    ramMin: 2048
  EnvVarRequirement:
    envDef:
      PYTHONPATH: /app
  DockerRequirement:
    dockerPull: 'nelly-acmg-agent:v1'
arguments:
  - prefix: --genes-tsv
    valueFrom: $(inputs.exomiser_genes_tsv.path)
  - prefix: --variants-tsv
    valueFrom: $(inputs.exomiser_variants_tsv.path)
  - prefix: --top
    valueFrom: $(inputs.top_candidates)
  - prefix: --output-vcf
    valueFrom: candidates.vcf
  # Reference FAI supplies ##contig lines and the chr-prefixed contig naming
  # VEP expects; omitted from the command line when not wired (gotcha 5).
  - prefix: --fai
    valueFrom: "$(inputs.reference_fai ? inputs.reference_fai.path : null)"
inputs:
  exomiser_genes_tsv:
    type: File
  exomiser_variants_tsv:
    type: File
  top_candidates:
    type: int
    default: 15   # a few more than acmg_top_candidates=10 so rank ties survive
  reference_fai:
    type: File?
outputs:
  candidates_vcf:
    type: File
    outputBinding:
      glob: candidates.vcf
