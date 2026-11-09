cwlVersion: v1.2
class: Workflow
requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  reads1: File
  reads2: File
  graph: File
  hapl: File
  ref_prefix:
    type: string
    default: GRCh38
  ref:
    type: File
    secondaryFiles:
      - .fai
  sample_name:
    type: string
    default: sample
  eh_variant_catalog: File?
  eh_sex:
    type: string
    default: female
  dv_mode:
    type: string
    default: standard
    doc: "'standard' runs google/deepvariant:1.10.0; 'pangenome_aware' runs google/deepvariant:pangenome_aware_deepvariant-1.10.0 with --pangenome=<graph>"
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
  eh_output_prefix:
    type: string
    default: eh_repeats

outputs:
  aligned_cram:
    type: File
    secondaryFiles:
      - .crai
    outputSource: samtools-markdup-index/aligned_reads_indexed
  dv_vcf:
    type: File
    secondaryFiles:
      - .tbi
    outputSource:
      - deepvariant/vcf
      - pangenome-aware-deepvariant/vcf
    pickValue: first_non_null
  dv_gvcf:
    type: File
    secondaryFiles:
      - .tbi
    outputSource:
      - deepvariant/gvcf
      - pangenome-aware-deepvariant/gvcf
    pickValue: first_non_null
  eh_vcf:
    type: File
    secondaryFiles:
      - .tbi
    outputSource: expansionhunter/vcf
  eh_json:
    type: File
    outputSource: expansionhunter/json

steps:
  vg-paths:
    run: vg-paths.cwl
    in:
      graph: graph
      ref_prefix: ref_prefix
    out: [path_list]

  kmc:
    run: kmc.cwl
    in:
      reads1: reads1
      reads2: reads2
      sample_name: sample_name
    out: [kff]

  vg-giraffe-bam:
    run: vg-giraffe-bam.cwl
    in:
      graph: graph
      hapl: hapl
      kff: kmc/kff
      ref_paths: vg-paths/path_list
      reads1: reads1
      reads2: reads2
      sample_name: sample_name
      bam_output: aligned_bam_output
    out: [aligned_reads]

  samtools-sort:
    run: samtools-sort.cwl
    in:
      aligned_reads: vg-giraffe-bam/aligned_reads
      ref: ref
      ref_prefix: ref_prefix
      sample_name: sample_name
    out: [aligned_reads_sorted]

  samtools-index-sorted:
    run: samtools-index.cwl
    in:
      aligned_reads: samtools-sort/aligned_reads_sorted
    out: [aligned_reads_indexed]

  samtools-markdup:
    run: samtools-markdup.cwl
    in:
      aligned_reads: samtools-index-sorted/aligned_reads_indexed
      ref: ref
    out: [marked_reads]

  samtools-markdup-index:
    run: samtools-index.cwl
    in:
      aligned_reads: samtools-markdup/marked_reads
    out: [aligned_reads_indexed]

  deepvariant:
    run: deepvariant.cwl
    when: $(inputs.dv_mode == "standard")
    in:
      dv_mode: dv_mode
      model_type: dv_model_type
      ref: ref
      aligned_reads: samtools-markdup-index/aligned_reads_indexed
      output_vcf: dv_output_vcf
      output_gvcf: dv_output_gvcf
    out: [vcf, gvcf]

  pangenome-aware-deepvariant:
    run: pangenome-aware-deepvariant.cwl
    when: $(inputs.dv_mode == "pangenome_aware")
    in:
      dv_mode: dv_mode
      model_type: dv_model_type
      ref: ref
      aligned_reads: samtools-markdup-index/aligned_reads_indexed
      graph: graph
      ref_name_pangenome: ref_prefix
      output_vcf: dv_output_vcf
      output_gvcf: dv_output_gvcf
    out: [vcf, gvcf]

  expansionhunter:
    run: expansionhunter.cwl
    in:
      aligned_reads: samtools-markdup-index/aligned_reads_indexed
      ref: ref
      variant_catalog: eh_variant_catalog
      sex: eh_sex
      output_prefix: eh_output_prefix
    out: [vcf, json, realigned_bam]
