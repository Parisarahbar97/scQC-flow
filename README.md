
# 🧬 scQC-flow

<img src="https://img.shields.io/badge/Nextflow-v22.10.0+-green.svg" alt="Nextflow Version">
<img src="https://img.shields.io/badge/Containers-Docker%2FSingularity-orange.svg" alt="Container Support">

**A Nextflow pipeline for single-cell RNA-seq quality control and reporting**

scQC-flow automates nuclear fraction analysis (DropletQC), doublet detection (scDblFinder), Seurat object creation (pre/post QC), and generates per-sample and combined HTML reports using Quarto.

---

## 🚀 Quick Start

### Installation

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
```

### Running the Pipeline (from the website)

```bash
nextflow run johnsonlab-ic/scQC-flow \
	--mapping_dirs <path/to/mapping_dirs.csv> \
	--outputDir <output_directory> \
	-profile <profile> \
	-resume
```

Replace `<path/to/mapping_dirs.csv>`, `<output_directory>`, and `<profile>` as needed. For example:

```bash
nextflow run johnsonlab-ic/scQC-flow \
	--mapping_dirs personal/test_mapping_dirs.csv \
	--outputDir personal/outs \
	-profile offline \
	-resume
```

---

## 📁 Input Files & Parameters

### Input and Output Paths

| Parameter         | Description                                 | Example/Default |
|-------------------|---------------------------------------------|-----------------|
| `--mapping_dirs`  | CSV with `samplename` and `path` columns    | `personal/mapping_dirs.csv` |
| `--outputDir`     | Output directory for results                | `results`       |

### Required Input Format

The mapping CSV should look like:

```csv
samplename,path
LM0094_SN1,/path/to/LM0094_SN1_snRNA_mapped
LM0094_SN2,/path/to/LM0094_SN2_snRNA_mapped
```

---

## 📋 Overview

The pipeline performs these key steps:

1. **Input Mapping**: Reads a CSV listing sample names and mapping directories (Cell Ranger outputs)
2. **DropletQC**: Computes nuclear fraction metrics for each sample
3. **scDblFinder**: Detects doublets and annotates cells
4. **Seurat Object Creation**: Builds Seurat objects (pre/post QC) with DropletQC and doublet metadata
5. **Reporting**: Generates per-sample and combined HTML reports using Quarto

---

## 📦 Outputs

| Output File/Folder                        | Description                                 |
|-------------------------------------------|---------------------------------------------|
| `*_dropletqc_metrics.csv`                 | Per-sample DropletQC metrics                |
| `*_scdbl_metrics.csv`                     | Per-sample scDblFinder metrics              |
| `*_qc_report.html`                        | Per-sample HTML QC reports                  |
| `combined_qc_book/_book/index.html`       | Combined multi-sample HTML report/book      |
| `*_preQC.rds`, `*_postQC.rds`             | Seurat RDS files (pre/post QC)              |

---

## ⚙️ Requirements

- Nextflow
- Docker or Singularity (for containers)

---

## 📊 Aggregate QC Reporting (new)

The script `scripts/qc_summary_report.R` summarises the pre/post QC Seurat objects produced by the pipeline and creates CSV tables, markdown summaries, and bar charts. A companion Docker image ensures the required R packages (Seurat, tidyverse, optparse, etc.) are available out-of-the-box.

### Build the container

```bash
docker build -t scqc-flow-report .
```

### Generate reports for a cohort

Bind-mount the project (for the script) and the directories that contain the per-sample outputs:

```bash
# Healthy cohort example
docker run --rm \
  -v /home/pr422/RDS/live/Users/Parisa/scQC-flow:/workspace \
  -v /home/pr422/RDS/live/Users/Parisa/EPILEP:/data \
  scqc-flow-report \
  --input /data/healthy/qc/output_latest \
  --output /data/QC_reports \
  --label healthy

# Diseased cohort example
docker run --rm \
  -v /home/pr422/RDS/live/Users/Parisa/scQC-flow:/workspace \
  -v /home/pr422/RDS/live/Users/Parisa/EPILEP:/data \
  scqc-flow-report \
  --input /data/diseased/qc/output_latest \
  --output /data/QC_reports \
  --label diseased
```

Each run writes:
- `*_qc_summary_per_sample.csv`, `*_qc_breakdown_per_sample.csv`, `*_qc_totals_all_samples.csv`, `*_qc_breakdown_totals.csv`
- `*_cells_after_qc_per_sample.png`
- `*_cells_pre_vs_post_per_sample.png`
- `*_qc_summary_report.md`

The `--max-mito`, `--min-nuclear`, `--min-features`, and `--min-counts` CLI flags mirror the thresholds used during QC and default to the pipeline's recommended values.

### Harmony Integration & Leiden clustering

Use `scripts/harmony_integration.R` to merge the post-QC Seurat objects, perform Harmony integration, run Leiden clustering, generate UMAPs, and export an AnnData copy. The script is container-friendly; the example below uses the Johnson Lab Landmark image you referenced, mounting the project plus healthy/diseased QC outputs:

```bash
docker run --rm \
  -v /home/pr422/RDS/live/Users/Parisa/scQC-flow:/workspace \
  -v /home/pr422/RDS/live/Users/Parisa/EPILEP:/data \
  -w /workspace \
  ghcr.io/johnsonlab-ic/landmark-sc_image:with-scanpy \
  Rscript scripts/harmony_integration.R \
    --healthy /data/healthy/qc/output_latest \
    --diseased /data/diseased/qc/output_latest \
    --output /data/harmony_integration \
    --dims 50 \
    --resolution 0.2
```

Outputs written to `--output` (e.g. `/home/pr422/RDS/live/Users/Parisa/EPILEP/harmony_integration`) include:
- `integrated_harmony_seurat.rds` (integrated Seurat object)
- `integrated_harmony_seurat.h5seurat` and `integrated_harmony_seurat.h5ad` (if SeuratDisk available)
- `integrated_metadata_with_umap.csv` (metadata + UMAP coordinates)
- UMAP PNGs coloured by clusters, cohort, and sample
- `integration_summary.txt` logging key parameters

### SingleR cell-type annotation

`scripts/singler_annotation.R` augments the integrated Seurat object with SingleR predictions using the `humanHippocampus2024` reference (all five annotation granularities). All outputs land in `/home/pr422/RDS/live/Users/Parisa/EPILEP/singleR_annotation` with a single container call:

```bash
docker run --rm \
  -v /home/pr422/RDS/live/Users/Parisa/scQC-flow:/workspace \
  -v /home/pr422/RDS/live/Users/Parisa/EPILEP:/data \
  -w /workspace \
  ghcr.io/johnsonlab-ic/landmark-sc_image:with-scanpy \
  bash -lc 'Rscript -e "if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org\"); \
                    BiocManager::install(c(\"SingleR\",\"SummarizedExperiment\",\"SingleCellExperiment\",\"scuttle\",\"ExperimentHub\",\"humanHippocampus2024\"), ask=FALSE, update=FALSE); \
                    install.packages(c(\"optparse\",\"magick\"), repos=\"https://cloud.r-project.org\")" && \
            Rscript scripts/singler_annotation.R \
              --seurat /data/harmony_integration/integrated_harmony_seurat.rds \
              --output /data/singleR_annotation \
              --individual-col Sample_ID \
              --default-level mid.cell.class \
              --cluster-col harmony_leiden \
              --levels "cell.type2,mid.cell.class,broad.cell.class" \
              --ncores 16 \
              --bp-type multicore'
```

Artifacts include:
- `integrated_with_singleR.rds` (Seurat object with `singler_*` annotations, `cell_type`, `individual`)
- Per-level label tables (`singleR_cells_*`, `singleR_counts_*`)
- UMAPs coloured by each annotation level (`umap_singler_*`)
- `singleR_summary.txt` documenting parameters

Flags of interest:
- `--cluster-col` (default `harmony_leiden`) lets SingleR annotate one cluster at a time (faster). Set to `''` for per-cell mode.
- `--levels` accepts a comma-separated list so you can run coarse levels first.
- `--ncores` + `--bp-type` control BiocParallel (use `snow` if fork isn't allowed).
- `--fine-tune` / `--prune` mirror SingleR options when you need more refinement.

---

## ⚠️ Notes & Warnings

> **System Requirements**: This pipeline is computationally intensive and best run on HPC systems. It is optimized for the Imperial College HPC system but can be adapted to other systems with sufficient resources.

---

## 📚 Repository Information

<img src="https://img.shields.io/badge/GitHub-scQC--flow-lightgrey?logo=github" alt="GitHub Repo">

This pipeline is maintained in a public repository:
- [johnsonlab-ic/scQC-flow](https://github.com/johnsonlab-ic/scQC-flow)

---

## 📞 Support

For questions or issues, please:
- Open an issue on the [GitHub repository](https://github.com/johnsonlab-ic/scQC-flow/issues)
- Contact the Johnson Lab at Imperial College London
