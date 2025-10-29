#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  ok_sw <- requireNamespace("SeuratWrappers", quietly = TRUE)
  ok_harmony <- requireNamespace("harmony", quietly = TRUE)
  ok_seuratdisk <- requireNamespace("SeuratDisk", quietly = TRUE)
  library(dplyr)
  library(ggplot2)
  library(glue)
})

option_list <- list(
  make_option("--healthy", type = "character", help = "Directory containing healthy *_seurat_object_postqc.rds files"),
  make_option("--diseased", type = "character", help = "Directory containing diseased *_seurat_object_postqc.rds files"),
  make_option("--output", type = "character", help = "Destination directory for integration outputs"),
  make_option("--variable-features", type = "integer", default = 3000L,
              help = "Number of HVGs to select prior to PCA [default %default]"),
  make_option("--dims", type = "integer", default = 40L,
              help = "Number of PCs / Harmony dimensions to use [default %default]"),
  make_option("--resolution", type = "double", default = 0.2,
              help = "Leiden clustering resolution [default %default]"),
  make_option("--umap-min-dist", type = "double", default = 0.3,
              help = "UMAP min.dist parameter [default %default]"),
  make_option("--assay", type = "character", default = "RNA",
              help = "Seurat assay to operate on (must contain counts layer) [default %default]"),
  make_option("--seed", type = "integer", default = 12345L,
              help = "Random seed for PCA/Harmony/UMAP [default %default]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Integrate post-QC Seurat objects with Harmony and cluster using Leiden.")
opts <- parse_args(parser)

if (is.null(opts$healthy) || is.null(opts$diseased) || is.null(opts$output)) {
  print_help(parser)
  stop("Missing required arguments --healthy, --diseased or --output", call. = FALSE)
}

stopifnot(dir.exists(opts$healthy), dir.exists(opts$diseased))
dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)

message("Healthy directory:  ", normalizePath(opts$healthy))
message("Diseased directory: ", normalizePath(opts$diseased))
message("Output directory:   ", normalizePath(opts$output))

list_postqc_files <- function(root) {
  list.files(root, pattern = "_seurat_object_postqc[.]rds$", recursive = TRUE, full.names = TRUE)
}

files_healthy <- list_postqc_files(opts$healthy)
files_diseased <- list_postqc_files(opts$diseased)

if (!length(files_healthy)) {
  stop("No *_seurat_object_postqc.rds files found under: ", opts$healthy)
}
if (!length(files_diseased)) {
  stop("No *_seurat_object_postqc.rds files found under: ", opts$diseased)
}

message("Found ", length(files_healthy), " healthy and ", length(files_diseased), " diseased post-QC files.")

read_and_label <- function(path, cohort) {
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) {
    stop("File does not contain a Seurat object: ", path)
  }

  DefaultAssay(obj) <- opts$assay

  samp <- basename(dirname(path))
  obj$Sample_ID <- samp
  obj$Cohort <- cohort

  assay_layers <- Layers(obj[[opts$assay]])
  if (is.null(assay_layers) || !"counts" %in% assay_layers) {
    stop("Assay '", opts$assay, "' is missing a 'counts' layer in ", path)
  }

  obj
}

objs <- c(
  lapply(files_healthy, read_and_label, cohort = "healthy"),
  lapply(files_diseased, read_and_label, cohort = "diseased")
)

message("Merging ", length(objs), " Seurat objects ...")
merged <- Reduce(function(a, b) merge(a, b), objs)
message("Merged object dimensions (features x cells): ", paste(dim(merged), collapse = " x "))
message("Cells per cohort:")
print(table(merged$Cohort))

counts_split <- any(grepl("^counts\\.", Layers(merged[[opts$assay]])))
if (!counts_split) {
  message("Splitting assay layers by Sample_ID ...")
  merged[[opts$assay]] <- split(merged[[opts$assay]], f = merged$Sample_ID)
} else {
  message("Assay already contains per-sample layers; continuing.")
}

set.seed(opts$seed)
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = opts$`variable-features`)
merged <- ScaleData(merged)
merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = opts$dims)

reduction_to_use <- "harmony"

if (!ok_harmony) {
  stop("The 'harmony' package is required but not installed in this environment.")
}

if (ok_sw && exists("HarmonyIntegration")) {
  message("Running IntegrateLayers(method = HarmonyIntegration) ...")
  merged <- IntegrateLayers(
    object = merged,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = reduction_to_use,
    dims = 1:opts$dims,
    verbose = TRUE
  )
} else {
  message("SeuratWrappers not found; falling back to harmony::RunHarmony ...")
  hm_fun <- getS3method("RunHarmony", "Seurat")
  hm_args <- names(formals(hm_fun))

  hm_call <- list(
    object = merged,
    `group.by.vars` = c("Sample_ID", "Cohort"),
    verbose = TRUE
  )

  if ("reduction" %in% hm_args) {
    hm_call$reduction <- "pca"
  } else {
    hm_call$`reduction.use` <- "pca"
    if ("reduction.save" %in% hm_args) {
      hm_call$`reduction.save` <- "harmony"
    }
  }

  if ("dims.use" %in% hm_args) {
    hm_call$`dims.use` <- 1:opts$dims
  } else {
    hm_call$dims <- 1:opts$dims
  }

  if ("assay.use" %in% hm_args) {
    hm_call$`assay.use` <- opts$assay
  } else if ("assay" %in% hm_args) {
    hm_call$assay <- opts$assay
  }

  hm_call$`project.dim` <- FALSE

  merged <- do.call(harmony::RunHarmony, hm_call)
}

if (!reduction_to_use %in% Reductions(merged)) {
  stop("Harmony reduction '", reduction_to_use, "' not present after integration.")
}

merged <- FindNeighbors(merged, reduction = reduction_to_use, dims = 1:opts$dims)
merged <- FindClusters(merged,
                       resolution = opts$resolution,
                       algorithm = 4,
                       verbose = TRUE,
                       cluster.name = "harmony_leiden")

merged <- RunUMAP(
  merged,
  reduction = reduction_to_use,
  dims = 1:opts$dims,
  return.model = FALSE,
  seed.use = opts$seed,
  min.dist = opts$`umap-min-dist`
)

obj_path <- file.path(opts$output, "integrated_harmony_seurat.rds")
saveRDS(merged, obj_path)
message("Saved integrated Seurat object: ", obj_path)

emb_umap <- Embeddings(merged, "umap")
meta_out <- merged@meta.data %>%
  mutate(cell_barcode = rownames(.)) %>%
  bind_cols(as.data.frame(emb_umap))

meta_path <- file.path(opts$output, "integrated_metadata_with_umap.csv")
write.csv(meta_out, meta_path, row.names = FALSE)
message("Saved metadata + UMAP coordinates: ", meta_path)

plots <- list(
  DimPlot(merged, reduction = "umap", group.by = "harmony_leiden", label = TRUE) + ggtitle("UMAP: Harmony/Leiden clusters"),
  DimPlot(merged, reduction = "umap", group.by = "Cohort") + ggtitle("UMAP: Cohort"),
  DimPlot(merged, reduction = "umap", group.by = "Sample_ID") + ggtitle("UMAP: Sample ID")
)

plot_names <- c("umap_clusters.png", "umap_cohort.png", "umap_sample.png")

for (i in seq_along(plots)) {
  ggplot2::ggsave(
    filename = file.path(opts$output, plot_names[i]),
    plot = plots[[i]],
    width = 8,
    height = 6,
    dpi = 300
  )
  message("Saved plot: ", plot_names[i])
}

if (ok_seuratdisk) {
  message("Exporting h5Seurat and h5ad ...")
  h5seurat_path <- file.path(opts$output, "integrated_harmony_seurat.h5seurat")
  h5ad_path <- file.path(opts$output, "integrated_harmony_seurat.h5ad")
  SeuratDisk::SaveH5Seurat(merged, filename = h5seurat_path, overwrite = TRUE)
  SeuratDisk::Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
  if (file.exists(h5ad_path)) {
    message("Saved AnnData file: ", h5ad_path)
  } else {
    message("Conversion completed; h5ad file should be located alongside the h5Seurat file.")
  }
} else {
  message("Package 'SeuratDisk' not installed; skipping h5Seurat/h5ad export.")
}

summary_csv <- file.path(opts$output, "integration_summary.txt")
summary_lines <- c(
  glue("Harmony integration completed on {Sys.time()}"),
  glue("Healthy files: {length(files_healthy)}"),
  glue("Diseased files: {length(files_diseased)}"),
  glue("Total cells: {ncol(merged)}"),
  glue("Number of Harmony clusters: {length(unique(merged$harmony_leiden))}"),
  glue("Parameters: HVGs={opts$`variable-features`}, dims={opts$dims}, resolution={opts$resolution}")
)
writeLines(summary_lines, summary_csv)
message("Integration summary written to: ", summary_csv)

message("Done.")
