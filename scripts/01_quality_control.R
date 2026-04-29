# ============================================================
# Script: 01_quality_control.R
# Project: Adipose scRNA-seq analysis
# Purpose: Perform sample-level quality control for scRNA-seq data
# Author: Yi Hu
# ============================================================

# This script provides a public skeleton of the QC workflow used for
# single-cell RNA-seq analysis. Raw and processed count matrices are
# available through GEO. 

# -----------------------------
# 1. Load packages
# -----------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(cowplot)
})

# -----------------------------
# 2. Define QC parameters
# -----------------------------

qc_params <- list(
  min_cells = 3,
  min_features_initial = 200,
  min_genes = 500,
  min_umis = 800,
  max_mito_percent = 18,
  min_novelty = 0.8
)

# -----------------------------
# 3. Define sample metadata
# -----------------------------

# Replace the example sample IDs and paths with the corresponding
# GEO-downloaded or locally stored count matrices.

sample_info <- tibble::tribble(
  ~sample_id,   ~condition, ~matrix_dir,
  "Sample_01",  "Chow",     "data/raw/Sample_01/",
  "Sample_02",  "HFSC",     "data/raw/Sample_02/"
)

# -----------------------------
# 4. Helper function for QC plots
# -----------------------------

plot_qc_metrics <- function(seurat_obj, sample_id, output_dir, params) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  p1 <- ggplot(seurat_obj[[]], aes(x = nFeature_RNA)) +
    geom_histogram(binwidth = 50) +
    geom_vline(
      xintercept = params$min_genes,
      color = "red",
      linetype = "dashed"
    ) +
    xlab("Number of detected genes")
  
  ggsave(
    filename = file.path(output_dir, paste0(sample_id, "_nFeature_histogram.png")),
    plot = p1,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  p2 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) +
    geom_hline(
      yintercept = params$min_genes,
      color = "red",
      linetype = "dashed"
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(sample_id, "_nFeature_violin.png")),
    plot = p2,
    width = 5,
    height = 4,
    dpi = 300
  )
  
  p3 <- ggplot(seurat_obj[[]], aes(x = nCount_RNA)) +
    geom_histogram(binwidth = 50) +
    geom_vline(
      xintercept = params$min_umis,
      color = "red",
      linetype = "dashed"
    ) +
    xlab("Number of UMIs") +
    xlim(0, 5000)
  
  ggsave(
    filename = file.path(output_dir, paste0(sample_id, "_nCount_histogram.png")),
    plot = p3,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  p4 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) +
    geom_hline(
      yintercept = params$min_umis,
      color = "red",
      linetype = "dashed"
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(sample_id, "_nCount_violin.png")),
    plot = p4,
    width = 5,
    height = 4,
    dpi = 300
  )
  
  p5 <- ggplot(seurat_obj[[]], aes(x = percent.mt)) +
    geom_histogram(binwidth = 1) +
    geom_vline(
      xintercept = params$max_mito_percent,
      color = "red",
      linetype = "dashed"
    ) +
    xlab("Mitochondrial reads (%)") +
    xlim(0, 100)
  
  ggsave(
    filename = file.path(output_dir, paste0(sample_id, "_percent_mt_histogram.png")),
    plot = p5,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  p6 <- ggplot(seurat_obj[[]], aes(x = log10(nCount_RNA), y = log10(nFeature_RNA), color = percent.mt)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_gradient(limits = c(0, 100)) +
    stat_smooth(method = "lm", se = FALSE) +
    xlab("log10 UMIs") +
    ylab("log10 detected genes")
  
  ggsave(
    filename = file.path(output_dir, paste0(sample_id, "_nCount_vs_nFeature.png")),
    plot = p6,
    width = 5,
    height = 4,
    dpi = 300
  )
  
  p7 <- ggplot(seurat_obj[[]], aes(x = novelty)) +
    geom_histogram(bins = 50) +
    xlab("Library novelty")
  
  ggsave(
    filename = file.path(output_dir, paste0(sample_id, "_novelty_histogram.png")),
    plot = p7,
    width = 6,
    height = 4,
    dpi = 300
  )
}

# -----------------------------
# 5. Helper function for one sample
# -----------------------------

run_sample_qc <- function(sample_id, condition, matrix_dir, params) {
  
  message("Processing sample: ", sample_id)
  
  # Load 10x Genomics count matrix
  counts <- Read10X(data.dir = matrix_dir)
  
  # Remove duplicated gene symbols if present
  counts <- counts[!duplicated(rownames(counts)), ]
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = params$min_cells,
    min.features = params$min_features_initial
  )
  
  # Add sample-level metadata
  seurat_obj$sample_id <- sample_id
  seurat_obj$condition <- condition
  
  # Calculate QC metrics
  seurat_obj[["total_reads"]] <- Matrix::colSums(GetAssayData(seurat_obj, slot = "counts"))
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj[["novelty"]] <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  
  initial_cells <- ncol(seurat_obj)
  
  # Save QC plots before filtering
  plot_qc_metrics(
    seurat_obj = seurat_obj,
    sample_id = paste0(sample_id, "_before_filtering"),
    output_dir = "results/qc/before_filtering",
    params = params
  )
  
  # Apply QC filtering
  seurat_obj_filtered <- subset(
    seurat_obj,
    subset =
      nFeature_RNA > params$min_genes &
      nCount_RNA > params$min_umis &
      percent.mt < params$max_mito_percent &
      novelty > params$min_novelty
  )
  
  final_cells <- ncol(seurat_obj_filtered)
  
  # Save QC plots after filtering
  plot_qc_metrics(
    seurat_obj = seurat_obj_filtered,
    sample_id = paste0(sample_id, "_after_filtering"),
    output_dir = "results/qc/after_filtering",
    params = params
  )
  
  # Summarize QC result
  qc_summary <- tibble::tibble(
    sample_id = sample_id,
    condition = condition,
    initial_cells = initial_cells,
    cells_after_filtering = final_cells,
    percent_retained = round(final_cells / initial_cells * 100, 2),
    median_nCount_RNA = median(seurat_obj_filtered$nCount_RNA),
    median_nFeature_RNA = median(seurat_obj_filtered$nFeature_RNA),
    median_percent_mt = median(seurat_obj_filtered$percent.mt),
    median_novelty = median(seurat_obj_filtered$novelty)
  )
  
  # Save filtered object
  # Large Seurat objects should not be committed to GitHub.
  # They should be saved locally or deposited in GEO.
  dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    seurat_obj_filtered,
    file = file.path("data/processed", paste0(sample_id, "_filtered_seurat.rds"))
  )
  
  return(qc_summary)
}

# -----------------------------
# 6. Run QC for all samples
# -----------------------------

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

qc_summary_list <- purrr::pmap(
  sample_info,
  function(sample_id, condition, matrix_dir) {
    run_sample_qc(
      sample_id = sample_id,
      condition = condition,
      matrix_dir = matrix_dir,
      params = qc_params
    )
  }
)

qc_summary <- dplyr::bind_rows(qc_summary_list)

write.csv(
  qc_summary,
  file = "results/qc/qc_summary.csv",
  row.names = FALSE
)

print(qc_summary)

# -----------------------------
# 7. Session information
# -----------------------------

sink("results/qc/session_info_QC.txt")
sessionInfo()
sink()