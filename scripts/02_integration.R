# ============================================================
# Script: 02_integration.R
# Project: Adipose scRNA-seq analysis
# Purpose: Integrate QC-filtered scRNA-seq samples using Seurat
# Author: Yi Hu
# ============================================================

# This script provides the core integration workflow used for the
# adipose scRNA-seq analysis. Large Seurat objects and raw count matrices
# are not included in this public repository.

# -----------------------------
# 1. Load packages
# -----------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
})

# -----------------------------
# 2. Define input and output paths
# -----------------------------

input_dir <- "data/processed"
output_dir <- "results/integration"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Example sample table.
# Replace these example sample IDs and file paths with locally generated
# QC-filtered Seurat objects or GEO-downloaded processed objects.

sample_info <- tibble::tribble(
  ~sample_id,   ~condition, ~file_path,
  "Sample_01",  "Chow",     file.path(input_dir, "Sample_01_filtered_seurat.rds"),
  "Sample_02",  "HFSC",     file.path(input_dir, "Sample_02_filtered_seurat.rds")
)

# -----------------------------
# 3. Load QC-filtered Seurat objects
# -----------------------------

seurat_list <- purrr::pmap(
  sample_info,
  function(sample_id, condition, file_path) {
    
    message("Loading sample: ", sample_id)
    
    obj <- readRDS(file_path)
    
    # Add sample-level metadata.
    obj$sample_id <- sample_id
    obj$condition <- condition
    
    return(obj)
  }
)

names(seurat_list) <- sample_info$sample_id

# -----------------------------
# 4. Normalize each sample with SCTransform
# -----------------------------

# SCTransform was used for normalization before integration.
# Each sample is normalized separately before anchor identification.

seurat_list <- lapply(
  X = seurat_list,
  FUN = function(obj) {
    SCTransform(
      object = obj,
      assay = "RNA",
      verbose = FALSE
    )
  }
)

# -----------------------------
# 5. Select integration features
# -----------------------------

integration_features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = 3000
)

# -----------------------------
# 6. Prepare SCT integration
# -----------------------------

seurat_list <- PrepSCTIntegration(
  object.list = seurat_list,
  anchor.features = integration_features
)

# -----------------------------
# 7. Find integration anchors
# -----------------------------

integration_anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = integration_features
)

# -----------------------------
# 8. Integrate samples
# -----------------------------

seurat_integrated <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT"
)

DefaultAssay(seurat_integrated) <- "integrated"

# -----------------------------
# 9. Dimensionality reduction
# -----------------------------

seurat_integrated <- RunPCA(
  object = seurat_integrated,
  verbose = FALSE
)

# A simple elbow plot is generated only to document PC selection.
# Final figure aesthetics were handled separately and are not included here.

pdf(file.path(output_dir, "elbow_plot.pdf"), width = 5, height = 4)
ElbowPlot(seurat_integrated, ndim = 30)
dev.off()

# Number of PCs used for downstream clustering.
# Adjust this value based on elbow plot and dataset complexity.

n_pcs <- 20

seurat_integrated <- RunUMAP(
  object = seurat_integrated,
  reduction = "pca",
  dims = 1:n_pcs
)

# -----------------------------
# 10. Graph-based clustering
# -----------------------------

seurat_integrated <- FindNeighbors(
  object = seurat_integrated,
  reduction = "pca",
  dims = 1:n_pcs,
  verbose = FALSE
)

# Multiple resolutions can be tested during exploratory analysis.
# The public workflow keeps the resolution used for downstream annotation.

selected_resolution <- 0.5

seurat_integrated <- FindClusters(
  object = seurat_integrated,
  resolution = selected_resolution,
  verbose = FALSE
)

# Store final cluster identities explicitly.

seurat_integrated$seurat_cluster <- Idents(seurat_integrated)

# -----------------------------
# 11. Minimal UMAP visualization
# -----------------------------

# These plots are intended only as quality-control outputs.
# Publication-quality figures were generated separately.

p_cluster <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  label = TRUE
) +
  ggtitle("Integrated clusters")

ggsave(
  filename = file.path(output_dir, "integrated_umap_by_cluster.pdf"),
  plot = p_cluster,
  width = 6,
  height = 5
)

p_condition <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "condition"
) +
  ggtitle("Integrated samples by condition")

ggsave(
  filename = file.path(output_dir, "integrated_umap_by_condition.pdf"),
  plot = p_condition,
  width = 6,
  height = 5
)

# -----------------------------
# 12. Identify cluster markers
# -----------------------------

# Cluster markers are used for downstream cell type annotation.

DefaultAssay(seurat_integrated) <- "SCT"

seurat_integrated <- PrepSCTFindMarkers(seurat_integrated)

Idents(seurat_integrated) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(
  object = seurat_integrated,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  cluster_markers,
  file = file.path(output_dir, "cluster_markers_all.csv"),
  row.names = FALSE
)

top20_cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE)

write.csv(
  top20_cluster_markers,
  file = file.path(output_dir, "cluster_markers_top20.csv"),
  row.names = FALSE
)

# -----------------------------
# 13. Save integrated object
# -----------------------------

# Large RDS files should be generated locally and should not be committed
# to the public GitHub repository.

saveRDS(
  seurat_integrated,
  file = file.path(output_dir, "seurat_integrated_sct.rds")
)

# -----------------------------
# 14. Session information
# -----------------------------

sink(file.path(output_dir, "session_info_integration.txt"))
sessionInfo()
sink()