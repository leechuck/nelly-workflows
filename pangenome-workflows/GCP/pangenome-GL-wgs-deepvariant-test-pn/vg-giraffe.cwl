cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: "quay.io/vgteam/vg:v1.78.0"
  ResourceRequirement:
    coresMin: 10
    ramMin: $(64 * 1024)
baseCommand: vg
arguments: [giraffe, -p, -t, $(runtime.cores)]

inputs:
  graph:
    type: File
    inputBinding:
      prefix: -Z
  ref_paths:
    type: File
    inputBinding:
      prefix: --ref-paths
  reads1:
    type: File
    inputBinding:
      prefix: -f
  reads2:
    type: File?
    inputBinding:
      prefix: -f
  dist:
    type: File
    inputBinding:
      prefix: -d
  min:
    type: File
    inputBinding:
      prefix: -m
  output:
    type: string
    default: output.gam
outputs:
  aligned_reads:
    type: stdout

stdout: $(inputs.output)
