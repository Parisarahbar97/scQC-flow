process CELLBENDER {
  label "process_cellbender"
  tag { sampleName }
  container "us.gcr.io/broad-dsde-methods/cellbender:latest"
  publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true

  input:
  tuple val(sampleName), path(mappingDir)

  output:
  tuple val(sampleName), path("${sampleName}_cellbender_output"), emit: cellbender_output
  path "${sampleName}_cellbender_output/cellbender_out.h5", emit: h5_file
  path "${sampleName}_cellbender_output/summary.txt", emit: summary

  script:
  """
  set -euo pipefail
  echo "Running CellBender (CPU) for sample: ${sampleName}"
  echo "Mapping directory: ${mappingDir}"

  mkdir -p "${sampleName}_cellbender_output"

  cellbender remove-background \\
    --input  "${mappingDir}/outs/raw_feature_bc_matrix.h5" \\
    --output "${sampleName}_cellbender_output/cellbender_out.h5"

  {
    echo "CellBender (CPU) completed"
    echo "sample: ${sampleName}"
    echo "note: ran with CellBender defaults (no expected-cells / TDI / fpr / epochs specified)"
  } > "${sampleName}_cellbender_output/summary.txt"

  echo "CellBender (CPU) completed for ${sampleName}"
  """
}

process CELLBENDER_GPU {
  label "process_gpu"
  tag { sampleName }
  container "us.gcr.io/broad-dsde-methods/cellbender:latest"
  publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true

  input:
  tuple val(sampleName), path(mappingDir)

  output:
  tuple val(sampleName), path("${sampleName}_cellbender_output"), emit: cellbender_output
  path "${sampleName}_cellbender_output/cellbender_out.h5", emit: h5_file
  path "${sampleName}_cellbender_output/summary.txt", emit: summary

  script:
  """
  set -euo pipefail
  echo "Running CellBender (GPU) for sample: ${sampleName}"
  echo "Mapping directory: ${mappingDir}"

  mkdir -p "${sampleName}_cellbender_output"

  cellbender remove-background \\
    --input  "${mappingDir}/outs/raw_feature_bc_matrix.h5" \\
    --output "${sampleName}_cellbender_output/cellbender_out.h5" \\
    --cuda

  {
    echo "CellBender (GPU) completed"
    echo "sample: ${sampleName}"
    echo "note: ran with CellBender defaults (no expected-cells / TDI / fpr / epochs specified)"
  } > "${sampleName}_cellbender_output/summary.txt"

  echo "CellBender (GPU) completed for ${sampleName}"
  """
}

process CELLBENDER_H5_CONVERT {
    label "process_low"
    tag { sampleName }
    container "ah3918/pilot-analyses:latest"
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sampleName), path(cellbender_output)

    output:
    tuple val(sampleName), path("${sampleName}_cellbender_output_seurat.h5"), emit: seurat_h5

    script:
    """
    set -euo pipefail
    echo "Converting CellBender H5 to Seurat-compatible format for: ${sampleName}"
    
    # Use ptrepack to create Seurat-compatible H5 file, overwriting nodes if needed
     ptrepack --complevel 5 ${cellbender_output}/cellbender_out_filtered.h5:/matrix ${sampleName}_cellbender_output_seurat.h5:/matrix

    echo "H5 conversion completed for ${sampleName}"
    """
}