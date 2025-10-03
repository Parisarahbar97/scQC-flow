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

    mkdir -p ${sampleName}_cellbender_output

    METRICS="${mappingDir}/outs/metrics_summary.csv"

    # 1) Read 'Estimated Number of Cells' (first data row, first column), strip quotes/commas
    if [ -f "\${METRICS}" ]; then
      EXPECTED=\$(awk -F',' 'NR==2{gsub(/[^0-9]/,"",\$1); print \$1}' "\${METRICS}")
    fi
    # Fallback if parsing failed or value empty
    [ -z "\${EXPECTED:-}" ] && EXPECTED=5000

    # 2) Choose total-droplets-included ~ 8x expected, with clamps for stability
    TDI=\$(( EXPECTED * 8 ))
    [ "\$TDI" -lt 20000 ] && TDI=20000
    [ "\$TDI" -gt 120000 ] && TDI=120000

    echo "Derived: expected-cells=\$EXPECTED  total-droplets-included=\$TDI"

    cellbender remove-background \\
        --input  ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
        --output ${sampleName}_cellbender_output/cellbender_out.h5 \\
        --expected-cells \$EXPECTED \\
        --total-droplets-included \$TDI \\
        --fpr 0.01 \\
        --epochs 150

    {
      echo "CellBender (CPU) completed"
      echo "sample: ${sampleName}"
      echo "expected_cells: \$EXPECTED"
      echo "total_droplets_included: \$TDI"
    } > ${sampleName}_cellbender_output/summary.txt

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

    mkdir -p ${sampleName}_cellbender_output

    METRICS="${mappingDir}/outs/metrics_summary.csv"

    # Read Estimated Number of Cells (2nd line, 1st column) and strip quotes/commas
    if [ -f "\${METRICS}" ]; then
      EXPECTED=\$(awk -F',' 'NR==2{gsub(/[^0-9]/,"",\$1); print \$1}' "\${METRICS}")
    fi
    # Fallback if missing
    [ -z "\${EXPECTED:-}" ] && EXPECTED=5000

    # total-droplets-included ≈ 8 × expected (clamped 20k–120k)
    TDI=\$(( EXPECTED * 8 ))
    [ "\$TDI" -lt 20000 ] && TDI=20000
    [ "\$TDI" -gt 120000 ] && TDI=120000

    echo "Derived: expected-cells=\$EXPECTED  total-droplets-included=\$TDI"

    cellbender remove-background \\
      --input  ${mappingDir}/outs/raw_feature_bc_matrix.h5 \\
      --output ${sampleName}_cellbender_output/cellbender_out.h5 \\
      --expected-cells \$EXPECTED \\
      --total-droplets-included \$TDI \\
      --fpr 0.01 \\
      --epochs 150 \\
      --cuda

    {
      echo "CellBender (GPU) completed"
      echo "sample: ${sampleName}"
      echo "expected_cells: \$EXPECTED"
      echo "total_droplets_included: \$TDI"
    } > ${sampleName}_cellbender_output/summary.txt

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
    echo "Converting CellBender H5 to Seurat-compatible format for: ${sampleName}"
    
    # Use ptrepack to create Seurat-compatible H5 file, overwriting nodes if needed
     ptrepack --complevel 5 ${cellbender_output}/cellbender_out_filtered.h5:/matrix ${sampleName}_cellbender_output_seurat.h5:/matrix

    echo "H5 conversion completed for ${sampleName}"
    """
}