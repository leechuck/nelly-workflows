cwlVersion: v1.2
class: Workflow
requirements:
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}
  ResourceRequirement:
    ramMin: $(1 * 1024)
    coresMin: 1
  MultipleInputFeatureRequirement: {}

inputs:
  ref_prefix:
    type: string
    default: GRCh38
  sample_name:
    type: string
    default: sample
  graph:
    type: File
  hapl:
    type: File
  dv_mode:
    type: string
    default: standard
  ref:
    type: File
    secondaryFiles:
      - .fai
  reads1: File
  reads2: File
  aligned_bam_output:
    type: string
    default: aligned.bam
  dv_model_type:
    type: string
    default: WGS
  dv_output_vcf:
    type: string
    default: variants.dv.vcf.gz
  dv_output_gvcf:
    type: string
    default: variants.dv.gvcf.gz
  eh_variant_catalog: File?
  eh_output_prefix:
    type: string
    default: eh_repeats
  combined_vcf_name:
    type: string
    default: combined.vcf.gz
  normalized_vcf_name:
    type: string
    default: normalized.vcf.gz
  slivar_gnomad: File
  slivar_info:
    type: string
    # Rare-variant gate: gnomAD popmax AF < 0.1% (was 1%) AND <=3 homozygotes in
    # gnomAD AND a PASS call. This is a rarity/quality gate, not an impact ranker
    # (slivar runs before VEP here, so no consequence annotation is available);
    # impact + phenotype ranking happens downstream (genome-linter / Exomiser).
    default: "INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt <= 3 && variant.FILTER == 'PASS'"
  slivar_sample_expr:
    type: string
    # Genotype-quality gate: keep only confidently-genotyped het/hom-alt calls and
    # tag them in INFO (high_quality=<sample>) for downstream inheritance reasoning.
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
  pavs_custom_args: string
  go_custom_file:
    type: File
    secondaryFiles:
      - .tbi
  go_custom_args: string
  hpo_custom_file:
    type: File
    secondaryFiles:
      - .tbi
  hpo_custom_args: string
  ppi_custom_file:
    type: File
    secondaryFiles:
      - .tbi
  ppi_custom_args: string
  vep_dir: Directory
  # --- Exomiser (phenotype-driven gene/variant prioritisation) ---
  exomiser_data:
    type: Directory          # Keep collection: data/2512_hg38 + data/2512_phenotype
  exomiser_assembly:
    type: string
    default: "hg38"
  exomiser_data_version:
    type: string
    default: "2512"
  hpo_terms:
    # Legacy fallback only. Genoor metadata with coded HPO phenotype URIs takes
    # precedence; this keeps older records with null URIs runnable.
    type: string
    default: ""
  exomiser_top_genes:
    type: int
    default: 10
  vep_fasta_file: File
  metadata: File
  genomelinter_output_name:
    type: string
    default: "genomelinter-output.txt"
  genomelinter_model:
    type: string
    default: "openai/gpt-oss-20b:free"
  openrouter_api_key:
    type: string
    default: ""

outputs:
  output_aligned_cram:
    type: File
    outputSource: align-and-call/aligned_cram
  output_dv_vcf:
    type: File
    outputSource: align-and-call/dv_vcf
  output_dv_gvcf:
    type: File
    outputSource: align-and-call/dv_gvcf
  output_eh_vcf:
    type: File
    outputSource: align-and-call/eh_vcf
  output_eh_json:
    type: File
    outputSource: align-and-call/eh_json
  output_combined_vcf:
    type: File
    outputSource: bcftools-norm/normalized_vcf
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
  output_exomiser_html:
    type: File
    outputSource: exomiser/exomiser_html
  output_after_genomelinter:
    type: File
    outputSource: genomelinter/genomelinter_output

steps:
  align-and-call:
    run: align-and-call.cwl
    in:
      reads1: reads1
      reads2: reads2
      graph: graph
      hapl: hapl
      ref_prefix: ref_prefix
      ref: ref
      sample_name: sample_name
      eh_variant_catalog: eh_variant_catalog
      eh_sex: phenofrommetadata/expansionhunter_sex
      dv_mode: dv_mode
      aligned_bam_output: aligned_bam_output
      dv_model_type: dv_model_type
      dv_output_vcf: dv_output_vcf
      dv_output_gvcf: dv_output_gvcf
      eh_output_prefix: eh_output_prefix
    out: [aligned_cram, dv_vcf, dv_gvcf, eh_vcf, eh_json]

  bcftools-concat:
    in:
      vcf_a: align-and-call/dv_vcf
      vcf_b: align-and-call/eh_vcf
      output_name: combined_vcf_name
    out: [combined_vcf]
    run: bcftools-concat.cwl

  bcftools-norm:
    in:
      input_vcf: bcftools-concat/combined_vcf
      ref: ref
      output_name: normalized_vcf_name
    out: [normalized_vcf]
    run: bcftools-norm.cwl

  slivar:
    in:
      slivar_gnomad: slivar_gnomad
      slivar_input: bcftools-norm/normalized_vcf
      slivar_info: slivar_info
      slivar_sample_expr: slivar_sample_expr
    out: [slivar_output]
    run: slivar.cwl

  vep:
    in:
      vep_input: [slivar/slivar_output]
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
    run: vep.cwl

  phenofrommetadata:
    in:
      metadata: metadata
      fallback_hpo_terms: hpo_terms
    out: [pheno_output, hpo_terms, expansionhunter_sex, exomiser_sex]
    run: phenofrommetadata.cwl

  # Strip records inside ExpansionHunter ReferenceRegions before Exomiser: DeepVariant
  # mis-calls these STR loci as coding frameshift indels (a spurious ATXN3 frameshift was
  # ranking #1 in every sample) and the genuine expansion is an unscorable symbolic <STR>
  # allele -- both pollute the ranking. The real repeat genotypes go to the genome-linter
  # via eh_repeats.json instead.
  exclude-eh-regions:
    in:
      input_vcf: bcftools-norm/normalized_vcf
      eh_vcf: align-and-call/eh_vcf
    out: [filtered_vcf]
    run: exclude-eh-regions.cwl

  # Phenotype-driven prioritisation: Exomiser ranks genes/variants from rarity +
  # pathogenicity + HPO similarity (grounded), replacing the phenotype-label-only
  # LLM ranker that buried the true causative gene.
  exomiser:
    in:
      exomiser_vcf: exclude-eh-regions/filtered_vcf
      exomiser_data: exomiser_data
      exomiser_assembly: exomiser_assembly
      exomiser_data_version: exomiser_data_version
      hpo_terms: phenofrommetadata/hpo_terms
      # Exomiser uses upper-case MALE/FEMALE for the same metadata value.
      sample_sex: phenofrommetadata/exomiser_sex
    out: [exomiser_genes_tsv, exomiser_variants_tsv, exomiser_json, exomiser_html]
    run: exomiser.cwl

  # The LLM writes a grounded narrative of Exomiser's ranked output PLUS any pathogenic-
  # range repeat expansion from ExpansionHunter (which Exomiser cannot score). Local
  # model, no OpenRouter; it cannot recall or bury genes.
  genomelinter:
    in:
      exomiser_genes_tsv: exomiser/exomiser_genes_tsv
      exomiser_variants_tsv: exomiser/exomiser_variants_tsv
      eh_json: align-and-call/eh_json
      phenotype: phenofrommetadata/pheno_output
      top_genes: exomiser_top_genes
    out: [genomelinter_output]
    run: genome-linter-llm.cwl
