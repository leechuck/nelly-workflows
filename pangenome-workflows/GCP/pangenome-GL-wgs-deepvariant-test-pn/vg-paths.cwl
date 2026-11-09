cwlVersion: v1.2
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: "quay.io/vgteam/vg:v1.78.0"
  ResourceRequirement:
    coresMin: 10
    ramMin: $(36 * 1024)
baseCommand: vg
arguments: [paths, -x, $(inputs.graph.path), -L, -Q, $(inputs.ref_prefix)]
inputs:
  graph:
    type: File
  ref_prefix:
    type: string
    default: GRCh38

outputs:
  path_list:
    type: stdout

stdout: $(inputs.ref_prefix).path_list.txt