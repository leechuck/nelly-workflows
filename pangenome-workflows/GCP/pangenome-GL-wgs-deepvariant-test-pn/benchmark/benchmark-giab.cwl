cwlVersion: v1.2
class: Workflow
requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}

doc: "GIAB benchmark workflow comparing align-and-call variant predictions against GIAB truth using hap.py"

inputs:
  graph:
    type: File
    doc: "Pangenome graph in GBZ format"
  hapl:
    type: File
    doc: "Haplotype index (.hapl) for giraffe haplotype sampling"
  ref:
    type: File
    secondaryFiles:
      - .fai
    doc: "Reference genome FASTA with .fai index"
  reads1:
    type: File
    doc: "Forward FASTQ reads (.fastq.gz)"
  reads2:
    type: File
    doc: "Reverse FASTQ reads (.fastq.gz)"
  ref_prefix:
    type: string
    default: "GRCh38"
    doc: "Reference prefix for graph paths"
  sample_name:
    type: string
    default: "sample"
    doc: "Sample identifier"
  dv_mode:
    type: string
    doc: "DeepVariant calling mode: standard or pangenome_aware"
  dv_model_type:
    type: string
    default: "WGS"
    doc: "DeepVariant model type"
  dv_output_vcf:
    type: string
    default: "variants.dv.vcf.gz"
    doc: "DeepVariant output VCF filename"
  dv_output_gvcf:
    type: string
    default: "variants.dv.gvcf.gz"
    doc: "DeepVariant output gVCF filename"
  eh_variant_catalog:
    type: File?
    doc: "Optional custom ExpansionHunter catalog"
  eh_output_prefix:
    type: string
    default: "eh_repeats"
    doc: "ExpansionHunter output prefix"
  eh_sex:
    type: string
    default: "male"
    doc: "Sample sex for ExpansionHunter (male|female)"
  truth_vcf:
    type: File
    secondaryFiles:
      - .tbi
    doc: "GIAB truth benchmark VCF with .tbi index"
  confident_bed:
    type: File
    doc: "GIAB confident call regions BED"
  stratification_tsv:
    type: File?
    secondaryFiles:
      - pattern: "^*"
        required: false
    doc: "Optional GIAB stratification TSV file (name <tab> bed_path)"
  stratification_beds:
    type: File[]?
    doc: "Optional stratification BED files referenced in stratification_tsv to stage into workdir"
  happy_prefix:
    type: string
    default: "happy"
    doc: "hap.py report output prefix"

outputs:
  aligned_cram:
    type: File
    outputSource: align-and-call/aligned_cram
  dv_vcf:
    type: File
    outputSource: align-and-call/dv_vcf
  dv_gvcf:
    type: File
    outputSource: align-and-call/dv_gvcf
  eh_vcf:
    type: File
    outputSource: align-and-call/eh_vcf
  eh_json:
    type: File
    outputSource: align-and-call/eh_json
  happy_summary_csv:
    type: File
    outputSource: happy/summary_csv
  happy_extended_csv:
    type: File
    outputSource: happy/extended_csv
  happy_metrics_json:
    type: File
    outputSource: happy/metrics_json
  happy_roc_all:
    type: File[]
    outputSource: happy/roc_all

steps:
  align-and-call:
    run: ../align-and-call.cwl
    in:
      graph: graph
      hapl: hapl
      ref: ref
      reads1: reads1
      reads2: reads2
      ref_prefix: ref_prefix
      sample_name: sample_name
      dv_mode: dv_mode
      dv_model_type: dv_model_type
      dv_output_vcf: dv_output_vcf
      dv_output_gvcf: dv_output_gvcf
      eh_variant_catalog: eh_variant_catalog
      eh_output_prefix: eh_output_prefix
      eh_sex: eh_sex
    out:
      - aligned_cram
      - dv_vcf
      - dv_gvcf
      - eh_vcf
      - eh_json

  happy:
    run: happy.cwl
    in:
      truth_vcf: truth_vcf
      query_vcf: align-and-call/dv_vcf
      confident_bed: confident_bed
      ref: ref
      prefix: happy_prefix
      stratification_tsv: stratification_tsv
      stratification_beds: stratification_beds
    out:
      - summary_csv
      - extended_csv
      - metrics_json
      - roc_all
