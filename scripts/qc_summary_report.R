#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(knitr)
  library(purrr)
  library(stringr)
  library(scales)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Root directory with sample subdirectories containing Seurat RDS outputs"),
  make_option(c("-o", "--output"), type = "character", default = ".", help = "Destination directory for summary artefacts [default: %default]"),
  make_option(c("-l", "--label"), type = "character", default = NULL, help = "Label prefix for output files (e.g. healthy, diseased)"),
  make_option("--max-mito", type = "double", default = 10.0, help = "Mitochondrial percentage threshold used during QC"),
  make_option("--min-nuclear", type = "double", default = 0.4, help = "Minimum nuclear fraction threshold used during QC"),
  make_option("--min-features", type = "integer", default = 200L, help = "Minimum feature threshold used during QC"),
  make_option("--min-counts", type = "integer", default = 500L, help = "Minimum count threshold used during QC"),
  make_option("--plot-width", type = "double", default = 14.0, help = "Width of the bar plot in inches"),
  make_option("--plot-height", type = "double", default = 7.0, help = "Height of the bar plot in inches"),
  make_option("--dpi", type = "integer", default = 200L, help = "Resolution (DPI) for saved figures")
)

parser <- OptionParser(option_list = option_list, description = "Summarise pre- and post-QC Seurat objects.")
opts <- parse_args(parser)

if (is.null(opts$input)) {
  print_help(parser)
  stop("Missing required argument: --input", call. = FALSE)
}

if (!dir.exists(opts$input)) {
  stop("Input directory does not exist: ", opts$input)
}

output_dir <- opts$output
if (is.null(output_dir) || !nzchar(output_dir)) {
  output_dir <- "."
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

label <- opts$label
prefix <- if (!is.null(label) && nzchar(label)) {
  paste0(label, "_")
} else {
  ""
}

message("Discovering samples under: ", opts$input)
cand_dirs <- list.dirs(opts$input, recursive = FALSE, full.names = TRUE)
samples <- basename(cand_dirs)

pre_rds <- file.path(cand_dirs, paste0(samples, "_seurat_object.rds"))
post_rds <- file.path(cand_dirs, paste0(samples, "_seurat_object_postqc.rds"))

keep <- file.exists(pre_rds) & file.exists(post_rds)

if (!any(keep)) {
  stop("No sample directories with paired Seurat objects were found under: ", opts$input)
}

cand_dirs <- cand_dirs[keep]
samples <- samples[keep]
pre_rds <- pre_rds[keep]
post_rds <- post_rds[keep]

rows <- vector("list", length(samples))
breakdown_rows <- vector("list", length(samples))

gcol <- function(df, candidates) {
  nm <- candidates[candidates %in% colnames(df)][1]
  if (is.na(nm)) {
    return(NULL)
  }
  df[[nm]]
}

for (i in seq_along(samples)) {
  sample <- samples[i]
  message("Processing sample: ", sample)

  pre <- readRDS(pre_rds[i])
  post <- readRDS(post_rds[i])

  pre_n <- ncol(pre)
  if (length(pre_n) == 0 || is.na(pre_n)) pre_n <- 0L
  post_n <- ncol(post)
  if (length(post_n) == 0 || is.na(post_n)) post_n <- 0L
  lost <- pre_n - post_n

  md <- pre@meta.data

  doublet_flag <- if (!is.null(md$doublet_class)) md$doublet_class == "doublet" else rep(FALSE, pre_n)
  nuclear_frac <- gcol(md, c("nuclear_fraction", "NuclearFraction", "nucleic_fraction"))
  nuclear_low <- if (!is.null(nuclear_frac)) nuclear_frac < opts$`min-nuclear` else rep(FALSE, pre_n)
  mito_pct <- gcol(md, c("percent.mt", "percent_mt"))
  mito_high <- if (!is.null(mito_pct)) mito_pct > opts$`max-mito` else rep(FALSE, pre_n)
  low_features <- if (!is.null(md$nFeature_RNA)) md$nFeature_RNA < opts$`min-features` else rep(FALSE, pre_n)
  low_counts <- if (!is.null(md$nCount_RNA)) md$nCount_RNA < opts$`min-counts` else rep(FALSE, pre_n)

  flagged_counts <- c(
    Doublet = sum(doublet_flag, na.rm = TRUE),
    Nuclear.Low = sum(nuclear_low, na.rm = TRUE),
    Mito.High = sum(mito_high, na.rm = TRUE),
    Low.Features = sum(low_features, na.rm = TRUE),
    Low.Counts = sum(low_counts, na.rm = TRUE)
  )

  lost_barcodes <- setdiff(colnames(pre), colnames(post))
  if (length(lost_barcodes) > 0) {
    lost_df <- data.frame(
      bc = lost_barcodes,
      Low.Features = if (!is.null(md$nFeature_RNA)) md[lost_barcodes, "nFeature_RNA"] < opts$`min-features` else FALSE,
      Low.Counts = if (!is.null(md$nCount_RNA)) md[lost_barcodes, "nCount_RNA"] < opts$`min-counts` else FALSE,
      Mito.High = if (!is.null(mito_pct)) mito_pct[lost_barcodes] > opts$`max-mito` else FALSE,
      Nuclear.Low = if (!is.null(nuclear_frac)) nuclear_frac[lost_barcodes] < opts$`min-nuclear` else FALSE,
      Doublet = if (!is.null(md$doublet_class)) md[lost_barcodes, "doublet_class"] == "doublet" else FALSE,
      stringsAsFactors = FALSE
    )

    lost_df$reason <- NA_character_
    order_names <- c("Low.Features", "Low.Counts", "Mito.High", "Nuclear.Low", "Doublet")

    for (nm in order_names) {
      sel <- is.na(lost_df$reason) & lost_df[[nm]]
      lost_df$reason[sel] <- nm
    }

    lost_df$reason[is.na(lost_df$reason)] <- "Other/NA"
    exclusive_counts <- table(factor(lost_df$reason, levels = c(order_names, "Other/NA")))
  } else {
    exclusive_counts <- setNames(rep(0, 6), c("Low.Features", "Low.Counts", "Mito.High", "Nuclear.Low", "Doublet", "Other/NA"))
  }

  rows[[i]] <- tibble(
    sample = sample,
    pre_cells = pre_n,
    post_cells = post_n,
    cells_lost = lost,
    pct_retained = round(100 * post_n / pre_n, 2),
    max_mito = opts$`max-mito`,
    min_nuclear = opts$`min-nuclear`,
    min_features = opts$`min-features`,
    min_counts = opts$`min-counts`
  )

  metrics_order <- c("Doublet", "Nuclear.Low", "Mito.High", "Low.Features", "Low.Counts", "Other/NA")
  flagged_vec <- c(
    as.integer(flagged_counts[c("Doublet", "Nuclear.Low", "Mito.High", "Low.Features", "Low.Counts")]),
    NA_integer_
  )
  flagged_pct_vec <- if (pre_n > 0) {
    c(round(100 * flagged_vec[seq_len(5)] / pre_n, 2), NA_real_)
  } else {
    rep(NA_real_, length(flagged_vec))
  }

  exclusive_vec <- as.integer(exclusive_counts[metrics_order])
  exclusive_pct_vec <- if (pre_n > 0) {
    round(100 * exclusive_vec / pre_n, 2)
  } else {
    rep(NA_real_, length(exclusive_vec))
  }

  breakdown_rows[[i]] <- tibble(
    sample = sample,
    metric = metrics_order,
    flagged_cells = flagged_vec,
    flagged_pct_of_pre = flagged_pct_vec,
    exclusive_lost_cells = exclusive_vec,
    exclusive_pct_of_pre = exclusive_pct_vec
  )

  rm(pre, post)
  gc()
}

summary_tbl <- bind_rows(rows) %>% arrange(sample)
breakdown_tbl <- bind_rows(breakdown_rows) %>%
  arrange(sample, metric)

total_pre <- sum(summary_tbl$pre_cells)
total_post <- sum(summary_tbl$post_cells)
total_lost <- total_pre - total_post

totals_tbl <- tibble(
  metric = c("Pre-QC cells", "Post-QC cells", "Cells retained", "Cells lost"),
  value = c(total_pre, total_post, total_post, total_lost),
  pct = c(NA_real_, NA_real_, round(100 * total_post / total_pre, 2), round(100 * total_lost / total_pre, 2))
)

breakdown_totals <- breakdown_tbl %>%
  group_by(metric) %>%
  summarise(
    flagged_cells = sum(flagged_cells, na.rm = TRUE),
    exclusive_lost_cells = sum(exclusive_lost_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    flagged_pct_of_all_pre = round(100 * flagged_cells / total_pre, 2),
    exclusive_pct_of_all_pre = round(100 * exclusive_lost_cells / total_pre, 2)
  ) %>%
  arrange(desc(flagged_cells))

# Write outputs -------------------------------------------------------------
write_csv(summary_tbl, file.path(output_dir, paste0(prefix, "qc_summary_per_sample.csv")))
write_csv(breakdown_tbl, file.path(output_dir, paste0(prefix, "qc_breakdown_per_sample.csv")))
write_csv(totals_tbl, file.path(output_dir, paste0(prefix, "qc_totals_all_samples.csv")))
write_csv(breakdown_totals, file.path(output_dir, paste0(prefix, "qc_breakdown_totals.csv")))

plot_title <- if (!is.null(label) && nzchar(label)) {
  paste0("Cells per Sample After QC (", label, ")")
} else {
  "Cells per Sample After QC"
}

p <- ggplot(summary_tbl, aes(x = sample, y = post_cells)) +
  geom_col(width = 0.7, fill = "#4472C4") +
  geom_text(aes(label = comma(post_cells)), vjust = -0.3, size = 3, color = "#1F1F1F") +
  labs(
    title = plot_title,
    x = "Sample",
    y = "Cells after QC"
  ) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_dir, paste0(prefix, "cells_after_qc_per_sample.png")),
  plot = p,
  width = opts$`plot-width`,
  height = opts$`plot-height`,
  dpi = opts$dpi
)

md_path <- file.path(output_dir, paste0(prefix, "qc_summary_report.md"))

summary_tbl_pretty <- summary_tbl %>%
  mutate(
    pre_cells = comma(pre_cells),
    post_cells = comma(post_cells),
    cells_lost = comma(cells_lost)
  )

breakdown_exclusive_wide <- breakdown_tbl %>%
  select(sample, metric, exclusive_lost_cells) %>%
  pivot_wider(names_from = metric, values_from = exclusive_lost_cells)

sink(md_path)
cat("# Multi-sample QC Summary\n\n")
cat("*Directory:* ", normalizePath(opts$input), "\n\n", sep = "")
if (!is.null(label) && nzchar(label)) {
  cat("*Group label:* ", label, "\n\n", sep = "")
}
cat("## Filtering thresholds\n\n")
cat(sprintf("- Max mitochondrial percentage: %.2f\n", opts$`max-mito`))
cat(sprintf("- Minimum nuclear fraction: %.2f\n", opts$`min-nuclear`))
cat(sprintf("- Minimum genes detected: %d\n", opts$`min-features`))
cat(sprintf("- Minimum counts: %d\n\n", opts$`min-counts`))

cat("## Totals across all samples\n\n")
cat(sprintf("- Pre-QC cells: %s\n", comma(total_pre)))
cat(sprintf("- Post-QC cells: %s\n", comma(total_post)))
cat(sprintf("- Cells retained: %s (%.2f%%)\n", comma(total_post), 100 * total_post / total_pre))
cat(sprintf("- Cells lost: %s (%.2f%%)\n\n", comma(total_lost), 100 * total_lost / total_pre))

cat("## Per-sample summary\n\n")
print(kable(summary_tbl_pretty, format = "pipe"))
cat("\n\n")

cat("## Reasons for filtering (exclusive counts)\n\n")
print(kable(breakdown_exclusive_wide, format = "pipe"))
cat("\n\n![Cells after QC](", paste0(prefix, "cells_after_qc_per_sample.png"), ")\n", sep = "")
sink()

message("All artefacts written to: ", normalizePath(output_dir))
