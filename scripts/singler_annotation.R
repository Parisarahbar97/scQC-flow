#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(ExperimentHub)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--seurat", type = "character", help = "Path to integrated Seurat .rds"),
  make_option("--output", type = "character", help = "Output directory for SingleR annotations"),
  make_option("--default-level", type = "character", default = "mid.cell.class",
              help = "Annotation level stored in `cell_type` [default %default]"),
  make_option("--individual-col", type = "character", default = "Sample_ID",
              help = "Meta-data column copied into `individual` [default %default]"),
  make_option("--ncores", type = "integer", default = 1L,
              help = "Cores for SingleR (use 1 inside most containers)"),
  make_option("--assay", type = "character", default = "RNA",
              help = "Seurat assay to operate on [default %default]")
)

parser <- OptionParser(
  option_list = option_list,
  description = "Annotate an integrated Seurat object with SingleR using the humanHippocampus2024 reference."
)
opts <- parse_args(parser)

if (is.null(opts$seurat) || is.null(opts$output)) {
  print_help(parser)
  stop("Missing required arguments: --seurat and --output", call. = FALSE)
}

dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# Load integrated data and convert to SCE
# ----------------------------------------------------------------------
message("Loading Seurat object: ", opts$seurat)
seu <- readRDS(opts$seurat)
DefaultAssay(seu) <- opts$assay

# Join layered assays (Seurat v5 multi-sample) into a single counts layer for downstream conversion
layer_names <- tryCatch(Layers(seu[[opts$assay]]), error = function(e) character(0))
if (length(layer_names) > 0 && any(grepl("^counts[.]", layer_names))) {
  message("Joining layered assay '", opts$assay, "' into a single counts matrix for SingleR ...")
  seu <- Seurat::JoinLayers(seu, assay = opts$assay)
  DefaultAssay(seu) <- opts$assay
}

to_sce <- function(seu, assay = "RNA") {
  sce <- tryCatch({
    as.SingleCellExperiment(seu)
  }, error = function(e) NULL)

  if (is.null(sce)) {
    counts <- tryCatch({
      Seurat::GetAssayData(seu[[assay]], layer = "counts")
    }, error = function(e) {
      Seurat::GetAssayData(seu[[assay]], slot = "counts")
    })
    sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)))
    colData(sce) <- DataFrame(seu@meta.data)
  }

  if (!"logcounts" %in% assayNames(sce)) {
    sce <- scuttle::logNormCounts(sce)
  }
  sce
}

sce_test <- to_sce(seu, assay = opts$assay)
message("Test SCE dims: ", paste(dim(sce_test), collapse = " x "))

# ----------------------------------------------------------------------
# Load reference (humanHippocampus2024 EH9606)
# ----------------------------------------------------------------------
message("Fetching humanHippocampus2024 reference (EH9606) ...")
eh <- ExperimentHub::ExperimentHub()
files <- ExperimentHub::query(eh, "humanHippocampus2024")
sce_ref <- files[["EH9606"]]
message("Reference SCE dims: ", paste(dim(sce_ref), collapse = " x "))

if (!"logcounts" %in% assayNames(sce_ref)) {
  sce_ref <- scuttle::logNormCounts(sce_ref)
}

# ----------------------------------------------------------------------
# Align genes between query and reference
# ----------------------------------------------------------------------
align_genes <- function(test_sce, ref_sce) {
  common <- intersect(rownames(test_sce), rownames(ref_sce))
  if (length(common) >= 200) {
    return(list(test = test_sce[common, ], ref = ref_sce[common, ]))
  }

  test_uc <- toupper(rownames(test_sce))
  ref_uc  <- toupper(rownames(ref_sce))
  common_uc <- intersect(test_uc, ref_uc)
  common_uc <- common_uc[!duplicated(common_uc)]

  if (length(common_uc) < 200) {
    stop("Too few overlapping genes between the dataset and reference. Ensure gene identifiers are compatible.")
  }

  test_idx <- match(common_uc, test_uc)
  ref_idx  <- match(common_uc, ref_uc)
  test_sce <- test_sce[test_idx, ]
  ref_sce  <- ref_sce[ref_idx, ]
  rownames(test_sce) <- rownames(ref_sce)
  list(test = test_sce, ref = ref_sce)
}

aligned <- align_genes(sce_test, sce_ref)
sce_test_al <- aligned$test
sce_ref_al  <- aligned$ref

# ----------------------------------------------------------------------
# Run SingleR at multiple annotation resolutions
# ----------------------------------------------------------------------
label_cols <- c(
  "cell.type2",
  "superfine.cell.class",
  "fine.cell.class",
  "mid.cell.class",
  "broad.cell.class"
)

missing_cols <- setdiff(label_cols, colnames(colData(sce_ref_al)))
if (length(missing_cols)) {
  stop("Reference is missing label columns: ", paste(missing_cols, collapse = ", "))
}

run_singler <- function(test_sce, ref_sce, labels, ncores = 1L) {
  if (!requireNamespace("SingleR", quietly = TRUE))
    stop("Package 'SingleR' is required.")
  if (!requireNamespace("BiocParallel", quietly = TRUE))
    stop("Package 'BiocParallel' is required.")

  bp <- if (ncores > 1) {
    BiocParallel::MulticoreParam(workers = ncores)
  } else {
    BiocParallel::SerialParam()
  }

  SingleR::SingleR(
    test = test_sce,
    ref = ref_sce,
    labels = labels,
    assay.type.test = "logcounts",
    assay.type.ref = "logcounts",
    BPPARAM = bp
  )
}

preds <- list()
for (col in label_cols) {
  message("Running SingleR for level: ", col)
  pred <- run_singler(
    test_sce = sce_test_al,
    ref_sce = sce_ref_al,
    labels = colData(sce_ref_al)[[col]],
    ncores = opts$ncores
  )
  preds[[col]] <- pred
}

add_pred_to_meta <- function(seu, pred, prefix) {
  md <- seu@meta.data
  md[[paste0("singler_", prefix)]] <- pred$labels
  md[[paste0("singler_", prefix, "_score")]] <- apply(pred$scores, 1, max, na.rm = TRUE)
  seu@meta.data <- md
  seu
}

for (col in names(preds)) {
  key <- gsub("[.]", "_", col)
  seu <- add_pred_to_meta(seu, preds[[col]], key)
}

# ----------------------------------------------------------------------
# Ensure required meta columns
# ----------------------------------------------------------------------
default_key <- gsub("[.]", "_", opts$`default-level`)
default_col <- paste0("singler_", default_key)
if (!default_col %in% colnames(seu@meta.data)) {
  stop("Default level column not found: ", default_col)
}
seu$cell_type <- seu[[default_col]][, 1, drop = TRUE]

ind_source <- opts$`individual-col`
if (!ind_source %in% colnames(seu@meta.data)) {
  warning("Requested individual column '", ind_source, "' not found. Using Sample_ID if available.")
  ind_source <- if ("Sample_ID" %in% colnames(seu@meta.data)) "Sample_ID" else NULL
}
seu$individual <- if (!is.null(ind_source)) seu[[ind_source]][, 1, drop = TRUE] else "unknown"

# ----------------------------------------------------------------------
# Save annotated Seurat object
# ----------------------------------------------------------------------
annot_path <- file.path(opts$output, "integrated_with_singleR.rds")
saveRDS(seu, annot_path)
message("Saved annotated Seurat object: ", annot_path)

# ----------------------------------------------------------------------
# Export per-level predictions and counts
# ----------------------------------------------------------------------
for (nm in names(preds)) {
  key <- gsub("[.]", "_", nm)
  pred <- preds[[nm]]
  per_cell <- data.frame(
    cell_barcode = rownames(pred),
    label = pred$labels,
    score = apply(pred$scores, 1, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    per_cell,
    file = file.path(opts$output, paste0("singleR_cells_", key, ".csv")),
    row.names = FALSE
  )

  counts <- as.data.frame(sort(table(pred$labels), decreasing = TRUE))
  colnames(counts) <- c("label", "n_cells")
  utils::write.csv(
    counts,
    file = file.path(opts$output, paste0("singleR_counts_", key, ".csv")),
    row.names = FALSE
  )
}

# ----------------------------------------------------------------------
# UMAP plots by annotation level
# ----------------------------------------------------------------------
if ("umap" %in% Reductions(seu)) {
  for (nm in names(preds)) {
    key <- gsub("[.]", "_", nm)
    coln <- paste0("singler_", key)
    plot_path <- file.path(opts$output, paste0("umap_singler_", key, ".png"))
    p <- DimPlot(
      seu,
      reduction = "umap",
      group.by = coln,
      label = FALSE,
      raster = TRUE
    ) +
      ggtitle(paste0("UMAP: SingleR ", nm)) +
      theme_minimal(base_size = 12)
    ggsave(plot_path, plot = p, width = 10, height = 8, dpi = 300)
    message("Saved UMAP: ", plot_path)
  }
} else {
  message("UMAP reduction not found; skipping UMAP plots.")
}

# ----------------------------------------------------------------------
# Summary log
# ----------------------------------------------------------------------
summary_file <- file.path(opts$output, "singleR_summary.txt")
summary_lines <- c(
  paste0("SingleR annotation completed: ", Sys.time()),
  paste0("Seurat input: ", opts$seurat),
  paste0("Output dir:  ", normalizePath(opts$output)),
  paste0("Default level assigned to `cell_type`: ", opts$`default-level`),
  paste0("`individual` column sourced from: ", opts$`individual-col`)
)
writeLines(summary_lines, summary_file)
message("Summary written to: ", summary_file)
