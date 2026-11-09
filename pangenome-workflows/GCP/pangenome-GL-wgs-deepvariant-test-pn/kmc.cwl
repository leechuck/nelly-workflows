cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "coolmaksat/kmc"
  ResourceRequirement:
    coresMin: 8
    ramMin: 65536
    tmpdirMin: 200000
  InitialWorkDirRequirement:
    listing:
      - entryname: reads.txt
        entry: |-
          $(inputs.reads1.path)
          $(inputs.reads2.path)

baseCommand: [kmc]
arguments:
  - -k29
  - -m64
  - -okff
  - -t$(runtime.cores)
  - "@reads.txt"
  - $(inputs.sample_name)
  - $(runtime.tmpdir)

inputs:
  reads1:
    type: File
  reads2:
    type: File
  sample_name:
    type: string
    default: sample

outputs:
  kff:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).kff
