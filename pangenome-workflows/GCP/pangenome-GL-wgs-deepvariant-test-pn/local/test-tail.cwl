#!/usr/bin/env cwl-runner
# Standalone validation of the modernised annotation and interpretation tail on an
# already-produced normalized.vcf.gz + ExpansionHunter outputs (eh_repeats.vcf.gz,
# eh_repeats.json) + patient metadata.
#
# Mirrors main-vg.cwl's tail wiring:
#   phenofrommetadata + exclude-eh-regions + slivar -> vep -> exomiser
#   (and, once WP4 is merged: -> acmg-extract-candidates -> vep-candidates -> acmg-agent)
cwlVersion: v1.2
class: Workflow
requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  NetworkAccess:
    networkAccess: true
  ResourceRequirement:
    ramMin: $(1 * 1024)
    coresMin: 1

inputs:
  sample_name:
    type: string
    default: sample
  normalized_vcf:
    type: File?
    secondaryFiles:
      - .tbi
  dv_vcf:
    type: File?
    secondaryFiles:
      - .tbi
  ref:
    type: File?
    secondaryFiles:
      - .fai
  combined_vcf_name:
    type: string
    default: combined.vcf.gz
  normalized_vcf_name:
    type: string
    default: normalized.vcf.gz
  eh_vcf:
    type: File
    secondaryFiles:
      - .tbi
  eh_json:
    type: File
  metadata:
    type: File
  hpo_terms:
    type: string
    default: ""
  slivar_gnomad:
    type: File
  slivar_info:
    type: string
    default: "INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt <= 3 && variant.FILTER == 'PASS'"
  slivar_sample_expr:
    type: string
    default: "high_quality:sample.GQ >= 20 && sample.DP >= 10 && (sample.het || sample.hom_alt)"
  vep_assembly:
    type: string
    default: GRCh38
  vep_output_file:
    type: string
    default: vep-output.vcf.gz
  pavs_custom_file:
    type: File
    secondaryFiles:
      - .tbi
  pavs_custom_args:
    type: string
    default: "PAVS,vcf,exact"
  go_custom_file:
    type: File
    secondaryFiles:
      - .tbi
  go_custom_args:
    type: string
    default: "GO_CLASSES,bed,overlap"
  hpo_custom_file:
    type: File
    secondaryFiles:
      - .tbi
  hpo_custom_args:
    type: string
    default: "PHENOTYPE,bed,overlap"
  ppi_custom_file:
    type: File
    secondaryFiles:
      - .tbi
  ppi_custom_args:
    type: string
    default: "PPI,bed,overlap"
  vep_dir:
    type: Directory
  vep_fasta_file:
    type: File
    secondaryFiles:
      - .fai
      - .gzi
  exomiser_data:
    type: Directory
  exomiser_assembly:
    type: string
    default: "hg38"
  exomiser_data_version:
    type: string
    default: "2512"
  exomiser_top_genes:
    type: int
    default: 10
  # ACMG evaluation agent + unimatrix LLM endpoint:
  reference_fai:
    type: File?
  acmg_clinvar_vcf:
    type: File
    secondaryFiles:
      - .tbi
  acmg_gene_table:
    type: File
  acmg_strchive_json:
    type: File?
  acmg_top_candidates:
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

outputs:
  output_filtered_vcf:
    type: File
    outputSource: exclude-eh-regions/filtered_vcf
  output_slivar:
    type: File
    outputSource: slivar/slivar_output
  output_after_vep:
    type: File
    outputSource: vep/vep_output
  output_exomiser_genes:
    type: File
    outputSource: exomiser/exomiser_genes_tsv
  output_exomiser_variants:
    type: File
    outputSource: exomiser/exomiser_variants_tsv
  output_exomiser_json:
    type: File
    outputSource: exomiser/exomiser_json
  output_exomiser_html:
    type: File
    outputSource: exomiser/exomiser_html
  output_candidates_vep:
    type: File
    outputSource: vep-candidates/vep_output
  output_acmg_json:
    type: File
    outputSource: acmg-agent/acmg_json
  output_acmg_report:
    type: File
    outputSource: acmg-agent/acmg_report

steps:
  bcftools-concat:
    when: $(inputs.vcf_a != null)
    in:
      vcf_a: dv_vcf
      vcf_b: eh_vcf
      output_name: combined_vcf_name
    out: [combined_vcf]
    run: ../bcftools-concat.cwl

  bcftools-norm:
    when: $(inputs.input_vcf != null)
    in:
      input_vcf: bcftools-concat/combined_vcf
      ref: ref
      output_name: normalized_vcf_name
    out: [normalized_vcf]
    run: ../bcftools-norm.cwl

  phenofrommetadata:
    in:
      metadata: metadata
      fallback_hpo_terms: hpo_terms
    out: [pheno_output, hpo_terms, expansionhunter_sex, exomiser_sex]
    run: ../phenofrommetadata.cwl

  exclude-eh-regions:
    in:
      input_vcf:
        source:
          - normalized_vcf
          - bcftools-norm/normalized_vcf
        pickValue: first_non_null
      eh_vcf: eh_vcf
    out: [filtered_vcf]
    run: ../exclude-eh-regions.cwl

  slivar:
    in:
      slivar_gnomad: slivar_gnomad
      slivar_input:
        source:
          - normalized_vcf
          - bcftools-norm/normalized_vcf
        pickValue: first_non_null
      slivar_info: slivar_info
      slivar_sample_expr: slivar_sample_expr
    out: [slivar_output]
    run: ../slivar.cwl

  vep:
    in:
      vep_input: slivar/slivar_output
      vep_assembly: vep_assembly
      vep_output_file: vep_output_file
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
    run: ../vep.cwl

  exomiser:
    in:
      exomiser_vcf: exclude-eh-regions/filtered_vcf
      exomiser_data: exomiser_data
      exomiser_assembly: exomiser_assembly
      exomiser_data_version: exomiser_data_version
      hpo_terms: phenofrommetadata/hpo_terms
      sample_sex: phenofrommetadata/exomiser_sex
    out: [exomiser_genes_tsv, exomiser_variants_tsv, exomiser_json, exomiser_html]
    run: ../exomiser.cwl

  acmg-extract-candidates:
    in:
      exomiser_genes_tsv: exomiser/exomiser_genes_tsv
      exomiser_variants_tsv: exomiser/exomiser_variants_tsv
      top_candidates: acmg_top_candidates
      reference_fai: reference_fai
    out: [candidates_vcf]
    run: ../acmg-extract-candidates.cwl

  vep-candidates:
    in:
      vep_input: acmg-extract-candidates/candidates_vcf
      vep_assembly: vep_assembly
      vep_output_file:
        default: "candidates.vep.vcf.gz"
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
    run: ../vep.cwl

  acmg-agent:
    in:
      vep_vcf: vep-candidates/vep_output
      clinvar_vcf: acmg_clinvar_vcf
      gene_table: acmg_gene_table
      strchive_json: acmg_strchive_json
      exomiser_genes_tsv: exomiser/exomiser_genes_tsv
      exomiser_variants_tsv: exomiser/exomiser_variants_tsv
      eh_json: eh_json
      slivar_vcf: slivar/slivar_output
      sample_name: sample_name
      hpo_terms: phenofrommetadata/hpo_terms
      phenotype: phenofrommetadata/pheno_output
      sample_sex: phenofrommetadata/exomiser_sex
      top_candidates: acmg_top_candidates
      llm_backend: llm_backend
      llm_endpoint: llm_endpoint
      llm_model: llm_model
      llm_api_key: llm_api_key
      llm_max_tokens: llm_max_tokens
    out: [acmg_json, acmg_report]
    run: ../acmg-agent.cwl
