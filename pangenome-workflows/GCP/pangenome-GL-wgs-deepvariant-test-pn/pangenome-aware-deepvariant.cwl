cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "google/deepvariant:pangenome_aware_deepvariant-1.10.0"
  ResourceRequirement:
    coresMin: 8
    ramMin: $(90 * 1024)
    tmpdirMin: 100000

# On Arvados, crunch-run launches containers with Docker's default 64 MB /dev/shm
# (--shm-size cannot be set). By default, run_pangenome_aware_deepvariant.py
# unconditionally runs /opt/deepvariant/bin/load_gbz_into_shared_memory whenever
# --pangenome ends in .gbz, which mmaps /dev/shm/GBZ_SHARED_MEMORY (~467 MB for
# ksa_hg38.gbz) and crashes with SIGBUS (Bus error) on a 64 MB /dev/shm.
#
# To run without /dev/shm on the cluster:
# 1. Pass --make_examples_extra_args=use_loaded_gbz_shared_memory=false so each
#    make_examples_pangenome_aware_dv shard passes --nouse_loaded_gbz_shared_memory
#    and loads the GBZ directly into normal process heap RAM instead of /dev/shm.
# 2. Stub out load_gbz_into_shared_memory_command in run_pangenome_aware_deepvariant
#    (which otherwise ignores use_loaded_gbz_shared_memory=false and still invokes
#    the shared-memory loader).
baseCommand:
  - python3
  - -u
  - -c
  - |
    import sys
    sys.path.insert(0, '/opt/deepvariant/bin')
    import run_pangenome_aware_deepvariant as r
    r.load_gbz_into_shared_memory_command = lambda *a, **kw: ('true', None)
    r.app.run(r.main)

inputs:
  model_type:
    type: string
    default: WGS
    inputBinding:
      prefix: --model_type
  ref:
    type: File
    secondaryFiles:
      - .fai
    inputBinding:
      prefix: --ref
  aligned_reads:
    type: File
    secondaryFiles:
      - .crai
    inputBinding:
      prefix: --reads
  graph:
    type: File
    doc: "GBZ pangenome graph used for alignment (e.g. ksa_hg38.gbz)"
    inputBinding:
      prefix: --pangenome
  ref_name_pangenome:
    type: string
    default: GRCh38
    doc: "Reference path prefix in the GBZ (our GBZ uses GRCh38#0#chr* with sample name GRCh38)"
    inputBinding:
      prefix: --ref_name_pangenome
  output_vcf:
    type: string
    default: variants.vcf.gz
    inputBinding:
      prefix: --output_vcf
  output_gvcf:
    type: string
    default: variants.gvcf.gz
    inputBinding:
      prefix: --output_gvcf
  num_shards:
    type: int
    default: 4
    doc: "Number of parallel make_examples_pangenome_aware_dv shards. Without /dev/shm, each shard loads ksa_hg38.gbz into process RAM (~12.7 GiB/shard measured on chr20); 4 shards (~51 GiB peak) stays safely under the 90 GiB container limit (8 shards would reach ~102 GiB and OOM)."
    inputBinding:
      prefix: --num_shards
  make_examples_extra_args:
    type: string
    default: "use_loaded_gbz_shared_memory=false"
    doc: "Extra arguments for make_examples_pangenome_aware_dv; disables /dev/shm GBZ loading by default for Arvados compatibility"
    inputBinding:
      prefix: --make_examples_extra_args
  intermediate_results_dir:
    type: string
    default: dv_intermediate
    inputBinding:
      prefix: --intermediate_results_dir
  regions:
    type: File?
    inputBinding:
      prefix: --regions

outputs:
  vcf:
    type: File
    secondaryFiles:
      - .tbi
    outputBinding:
      glob: $(inputs.output_vcf)
  gvcf:
    type: File
    secondaryFiles:
      - .tbi
    outputBinding:
      glob: $(inputs.output_gvcf)
