cwlVersion: v1.2
class: Workflow

doc: "GIAB benchmark workflow comparing legacy pipeline variant predictions (vg 1.54 Giraffe with dist/min + DV 1.6) against GIAB truth using hap.py"

inputs:
  graph:
    type: File
    doc: "Pangenome graph in GBZ format"
  dist:
    type: File
    doc: "Giraffe distance index (.dist)"
  min:
    type: File
    doc: "Giraffe minimizer index (.min)"
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
  aligned_bam_output:
    type: string
    default: "aligned.bam"
    doc: "Intermediate BAM output filename"
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
    outputSource: samtools-markdup-index/aligned_reads_indexed
  dv_vcf:
    type: File
    outputSource: deepvariant/vcf
  dv_gvcf:
    type: File
    outputSource: deepvariant/gvcf
  eh_vcf:
    type: File
    outputSource: expansionhunter/vcf
  eh_json:
    type: File
    outputSource: expansionhunter/json
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
  vg-paths:
    run: ../legacy/vg-paths.cwl
    in:
      graph: graph
      ref_prefix: ref_prefix
    out: [path_list]

  vg-giraffe-bam:
    run: ../legacy/vg-giraffe-bam.cwl
    in:
      graph: graph
      reads1: reads1
      reads2: reads2
      ref_paths: vg-paths/path_list
      dist: dist
      min: min
      sample_name: sample_name
      bam_output: aligned_bam_output
    out: [aligned_reads]

  samtools-sort:
    run: ../legacy/samtools-sort.cwl
    in:
      aligned_reads: vg-giraffe-bam/aligned_reads
      ref: ref
      ref_prefix: ref_prefix
      sample_name: sample_name
    out: [aligned_reads_sorted]

  samtools-index-sorted:
    run: ../legacy/samtools-index.cwl
    in:
      aligned_reads: samtools-sort/aligned_reads_sorted
    out: [aligned_reads_indexed]

  samtools-markdup:
    run: ../legacy/samtools-markdup.cwl
    in:
      aligned_reads: samtools-index-sorted/aligned_reads_indexed
      ref: ref
    out: [marked_reads]

  samtools-markdup-index:
    run: ../legacy/samtools-index.cwl
    in:
      aligned_reads: samtools-markdup/marked_reads
    out: [aligned_reads_indexed]

  deepvariant:
    run: ../legacy/deepvariant.cwl
    in:
      model_type: dv_model_type
      ref: ref
      aligned_reads: samtools-markdup-index/aligned_reads_indexed
      output_vcf: dv_output_vcf
      output_gvcf: dv_output_gvcf
    out: [vcf, gvcf]

  expansionhunter:
    run: ../legacy/expansionhunter.cwl
    in:
      aligned_reads: samtools-markdup-index/aligned_reads_indexed
      ref: ref
      variant_catalog: eh_variant_catalog
      sex: eh_sex
      output_prefix: eh_output_prefix
    out: [vcf, json]

  happy:
    run: happy.cwl
    in:
      truth_vcf: truth_vcf
      query_vcf: deepvariant/vcf
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
