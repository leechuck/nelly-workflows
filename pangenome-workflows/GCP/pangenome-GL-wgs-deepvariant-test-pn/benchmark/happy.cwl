cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "jmcdani20/hap.py:v0.3.12"
  ResourceRequirement:
    coresMin: 8
    ramMin: 32768
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: |
      ${
        var entries = [];
        if (inputs.stratification_tsv) {
          entries.push(inputs.stratification_tsv);
        }
        if (inputs.stratification_beds) {
          for (var i = 0; i < inputs.stratification_beds.length; i++) {
            entries.push(inputs.stratification_beds[i]);
          }
        }
        return entries;
      }

inputs:
  truth_vcf:
    type: File
    secondaryFiles:
      - .tbi
    doc: "GIAB truth VCF (bgzipped)"
  query_vcf:
    type: File
    secondaryFiles:
      - .tbi
    doc: "Query VCF from caller (bgzipped)"
  confident_bed:
    type: File
    doc: "GIAB confident call regions BED"
  ref:
    type: File
    secondaryFiles:
      - .fai
    doc: "Reference genome FASTA"
  prefix:
    type: string
    default: "happy"
    doc: "Output report prefix"
  stratification_tsv:
    type: File?
    secondaryFiles:
      - pattern: "^*"
        required: false
    inputBinding:
      position: 9
      prefix: --stratification
    doc: "Optional GIAB stratification TSV file (name <tab> bed_path)"
  stratification_beds:
    type: File[]?
    doc: "Optional stratification BED files referenced in stratification_tsv to stage into workdir"

baseCommand:
  - /opt/hap.py/bin/hap.py

arguments:
  - position: 1
    valueFrom: $(inputs.truth_vcf.path)
  - position: 2
    valueFrom: $(inputs.query_vcf.path)
  - position: 3
    prefix: -f
    valueFrom: $(inputs.confident_bed.path)
  - position: 4
    prefix: -r
    valueFrom: $(inputs.ref.path)
  - position: 5
    prefix: -o
    valueFrom: $(inputs.prefix)
  - position: 6
    valueFrom: --engine=vcfeval
  - position: 7
    valueFrom: --pass-only
  - position: 8
    prefix: --threads
    valueFrom: $(runtime.cores)

outputs:
  summary_csv:
    type: File
    outputBinding:
      glob: $(inputs.prefix).summary.csv
  extended_csv:
    type: File
    outputBinding:
      glob: $(inputs.prefix).extended.csv
  metrics_json:
    type: File
    outputBinding:
      glob: $(inputs.prefix).metrics.json.gz
  roc_all:
    type: File[]
    outputBinding:
      glob: $(inputs.prefix).roc.*
