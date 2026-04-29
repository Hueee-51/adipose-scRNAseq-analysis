# ============================================================
# Script: 03_clustering_annotation.R
# Project: Adipose snRNA-seq analysis
# Purpose: Annotate major cell types after Seurat integration
# Author: Yi Hu
# ============================================================

# This script provides the core workflow for annotating major cell types
# from integrated snRNA-seq clusters. Cluster annotation was based on
# canonical marker gene expression and cluster-level marker analysis.
#
# Input:
#   results/integration/seurat_integrated_sct.rds
#
# Output:
#   results/annotation/seurat_integrated_annotated.rds
#   results/annotation/cluster_to_celltype_annotation.csv
#   results/annotation/celltype_counts_by_condition.csv
#   results/annotation/celltype_fraction_by_condition.csv
#   results/annotation/major_celltype_markers_all.csv
#   results/annotation/major_celltype_markers_top20.csv

# -----------------------------
# 1. Load packages
# -----------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(forcats)
})

# -----------------------------
# 2. Define input and output paths
# -----------------------------

input_file <- "results/integration/seurat_integrated_sct.rds"
output_dir <- "results/annotation"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3. Load integrated Seurat object
# -----------------------------

seurat_obj <- readRDS(input_file)

# Use the cluster identity generated in the integration step.
Idents(seurat_obj) <- "seurat_clusters"

# -----------------------------
# 4. Review canonical marker genes
# -----------------------------

# Marker genes used to support major cell type annotation.
# This list is intentionally concise and does not include all exploratory
# marker checks performed during the full analysis.

major_marker_genes <- c(
  # Adipose stem and progenitor cells
  "Pdgfra", "Dpp4", "Ly6a",
  
  # Endothelial cells
  "Pecam1", "Cdh5", "Kdr",
  
  # Macrophages / myeloid cells
  "Ptprc", "Lyz2", "Adgre1", "Mafb",
  
  # T cells
  "Cd3d", "Cd3e",
  
  # B cells
  "Ms4a1", "Cd79a",
  
  # Mesothelial cells
  "Msln", "Wt1",
  
  # Dendritic cells
  "Flt3", "Itgax",
  
  # Neutrophils
  "S100a8", "S100a9"
)

major_marker_genes <- major_marker_genes[
  major_marker_genes %in% rownames(seurat_obj)
]

# Minimal marker dot plot for annotation record.
# Publication-quality figure styling was performed separately.

p_marker_dotplot <- DotPlot(
  object = seurat_obj,
  features = major_marker_genes
) +
  RotatedAxis() +
  ggtitle("Canonical markers for major cell type annotation")

ggsave(
  filename = file.path(output_dir, "major_marker_dotplot.pdf"),
  plot = p_marker_dotplot,
  width = 10,
  height = 5
)

# -----------------------------
# 5. Define cluster-to-cell type annotation
# -----------------------------

# Cluster annotation was assigned based on canonical marker genes and
# cluster-level marker analysis.
# Clusters sharing the same canonical marker profile were merged into
# biologically interpretable major cell types.

cluster_annotation <- tibble::tribble(
  ~cluster, ~cell_type,
  "0",      "ASPC",
  "1",      "Endothelial Cell",
  "2",      "Endothelial Cell",
  "3",      "T Cell",
  "4",      "B Cell",
  "5",      "Macrophage",
  "6",      "Macrophage",
  "7",      "Macrophage",
  "8",      "Macrophage",
  "9",      "Endothelial Cell",
  "10",     "ASPC",
  "11",     "Endothelial Cell",
  "12",     "Macrophage",
  "13",     "Mesothelial Cell",
  "14",     "Macrophage",
  "15",     "ASPC",
  "16",     "Dendritic Cell",
  "17",     "Macrophage",
  "18",     "B Cell",
  "19",     "Neutrophil",
  "20",     "B Cell"
)

write.csv(
  cluster_annotation,
  file = file.path(output_dir, "cluster_to_celltype_annotation.csv"),
  row.names = FALSE
)

# -----------------------------
# 6. Add cell type annotation to metadata
# -----------------------------

metadata <- seurat_obj@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  left_join(
    cluster_annotation,
    by = c("seurat_clusters" = "cluster")
  ) %>%
  tibble::column_to_rownames("barcode")

metadata$cell_type <- forcats::fct_explicit_na(
  metadata$cell_type,
  na_level = "Unknown"
)

seurat_obj$cell_type <- metadata[
  rownames(seurat_obj@meta.data),
  "cell_type"
]

seurat_obj$cell_type <- factor(
  seurat_obj$cell_type,
  levels = c(
    "ASPC",
    "Endothelial Cell",
    "Macrophage",
    "B Cell",
    "T Cell",
    "Mesothelial Cell",
    "Dendritic Cell",
    "Neutrophil",
    "Unknown"
  )
)

# -----------------------------
# 7. Create condition-specific DE groups
# -----------------------------

# The condition variable should have been added during QC/integration.
# Example DE_group values:
#   ASPC-Chow
#   ASPC-HFSC
#   Endothelial Cell-Chow
#   Endothelial Cell-HFSC

seurat_obj$DE_group <- paste0(
  seurat_obj$cell_type,
  "-",
  seurat_obj$condition
)

# -----------------------------
# 8. Summarize cell type composition
# -----------------------------

celltype_counts <- seurat_obj@meta.data %>%
  count(condition, cell_type, name = "n_cells")

write.csv(
  celltype_counts,
  file = file.path(output_dir, "celltype_counts_by_condition.csv"),
  row.names = FALSE
)

celltype_fraction <- celltype_counts %>%
  group_by(condition) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  celltype_fraction,
  file = file.path(output_dir, "celltype_fraction_by_condition.csv"),
  row.names = FALSE
)

# -----------------------------
# 9. Minimal UMAP visualization after annotation
# -----------------------------

# These plots are intended only as workflow/QC outputs.
# Final manuscript figures were customized separately.

Idents(seurat_obj) <- "cell_type"

p_celltype <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE
) +
  ggtitle("Annotated major cell types")

ggsave(
  filename = file.path(output_dir, "umap_annotated_celltypes.pdf"),
  plot = p_celltype,
  width = 6,
  height = 5
)

p_celltype_condition <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  split.by = "condition",
  label = TRUE
) +
  ggtitle("Annotated major cell types by condition")

ggsave(
  filename = file.path(output_dir, "umap_annotated_celltypes_by_condition.pdf"),
  plot = p_celltype_condition,
  width = 10,
  height = 5
)

# -----------------------------
# 10. Identify major cell type markers
# -----------------------------

# This marker table documents genes enriched in each annotated major
# cell type and supports the annotation.

DefaultAssay(seurat_obj) <- "RNA"

# For Seurat v5 objects, join RNA layers before DE/marker analysis if needed.
# If the object does not contain multiple layers, this step can be skipped.

seurat_obj <- JoinLayers(seurat_obj)

seurat_obj <- NormalizeData(
  object = seurat_obj,
  assay = "RNA",
  verbose = FALSE
)

Idents(seurat_obj) <- "cell_type"

major_celltype_markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

write.csv(
  major_celltype_markers,
  file = file.path(output_dir, "major_celltype_markers_all.csv"),
  row.names = FALSE
)

top20_major_celltype_markers <- major_celltype_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE)

write.csv(
  top20_major_celltype_markers,
  file = file.path(output_dir, "major_celltype_markers_top20.csv"),
  row.names = FALSE
)

# -----------------------------
# 11. Save annotated object
# -----------------------------

# This object is used as the input for downstream subset analysis,
# differential expression, pathway analysis, and cell-cell communication.

saveRDS(
  seurat_obj,
  file = file.path(output_dir, "seurat_integrated_annotated.rds")
)

# -----------------------------
# 12. Session information
# -----------------------------

sink(file.path(output_dir, "session_info_annotation.txt"))
sessionInfo()
sink()