# ============================================================
# Script: 09_cellchat_subcluster_analysis.R
# Project: Adipose snRNA-seq analysis
# Purpose: CellChat analysis at the subcluster level
# Author: Yi Hu
# ============================================================

# This script performs CellChat analysis using cell type subclusters.
# ASPC, endothelial cell, and macrophage subcluster annotations are
# combined with other major cell types before CellChat analysis.
#
# Input:
#   results/annotation/seurat_integrated_annotated.rds
#   results/aspc_subcluster/seurat_aspc_subcluster_annotated.rds
#   results/ec_subcluster/seurat_ec_subcluster_annotated.rds
#   results/macrophage_subcluster/seurat_macrophage_subcluster_annotated.rds
#
# Output:
#   results/cellchat_subcluster/cellchat_Chow_subcluster.rds
#   results/cellchat_subcluster/cellchat_HFSC_subcluster.rds
#   results/cellchat_subcluster/cellchat_merged_subcluster.rds
#   results/cellchat_subcluster/communication_tables/
#   results/cellchat_subcluster/pathway_tables/

# -----------------------------
# 1. Load packages
# -----------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(patchwork)
  library(future)
})

# -----------------------------
# 2. Define input and output paths
# -----------------------------

major_object_file <- "results/annotation/seurat_integrated_annotated.rds"
aspc_object_file <- "results/aspc_subcluster/seurat_aspc_subcluster_annotated.rds"
ec_object_file <- "results/ec_subcluster/seurat_ec_subcluster_annotated.rds"
macrophage_object_file <- "results/macrophage_subcluster/seurat_macrophage_subcluster_annotated.rds"

output_dir <- "results/cellchat_subcluster"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "communication_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "pathway_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3. Load annotated objects
# -----------------------------

seurat_major <- readRDS(major_object_file)
aspc_obj <- readRDS(aspc_object_file)
ec_obj <- readRDS(ec_object_file)
macrophage_obj <- readRDS(macrophage_object_file)

# -----------------------------
# 4. Merge subcluster annotations into major object
# -----------------------------

# Start with major cell type labels.
seurat_major$cellchat_group <- as.character(seurat_major$cell_type)

# Replace ASPC major labels with ASPC subcluster labels.
aspc_meta <- aspc_obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("barcode") %>%
  select(barcode, aspc_subcluster)

ec_meta <- ec_obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("barcode") %>%
  select(barcode, ec_subcluster)

macrophage_meta <- macrophage_obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("barcode") %>%
  select(barcode, macrophage_subcluster)

major_meta <- seurat_major@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("barcode") %>%
  left_join(aspc_meta, by = "barcode") %>%
  left_join(ec_meta, by = "barcode") %>%
  left_join(macrophage_meta, by = "barcode") %>%
  mutate(
    cellchat_group = case_when(
      !is.na(aspc_subcluster) ~ aspc_subcluster,
      !is.na(ec_subcluster) ~ ec_subcluster,
      !is.na(macrophage_subcluster) ~ macrophage_subcluster,
      TRUE ~ as.character(cell_type)
    )
  ) %>%
  column_to_rownames("barcode")

seurat_major$cellchat_group <- major_meta[
  rownames(seurat_major@meta.data),
  "cellchat_group"
]

# Remove remaining unassigned or unknown groups if present.
seurat_major <- subset(
  x = seurat_major,
  subset = !cellchat_group %in% c("Unassigned", "Unknown", NA)
)

seurat_major$cellchat_group <- factor(seurat_major$cellchat_group)

# Save group composition table.
cellchat_group_counts <- seurat_major@meta.data %>%
  count(condition, cellchat_group, name = "n_cells")

write.csv(
  cellchat_group_counts,
  file = file.path(output_dir, "cellchat_subcluster_group_counts_by_condition.csv"),
  row.names = FALSE
)

# -----------------------------
# 5. Define CellChat helper function
# -----------------------------

run_cellchat <- function(seurat_object, group_by, condition_name, output_dir) {
  
  message("Running CellChat for: ", condition_name)
  
  DefaultAssay(seurat_object) <- "RNA"
  
  data_input <- GetAssayData(
    object = seurat_object,
    assay = "RNA",
    slot = "counts"
  )
  
  meta <- seurat_object@meta.data
  meta <- meta[colnames(data_input), , drop = FALSE]
  
  meta_cellchat <- data.frame(
    row.names = rownames(meta),
    celltype = meta[[group_by]]
  )
  
  cellchat <- createCellChat(
    object = data_input,
    meta = meta_cellchat,
    group.by = "celltype"
  )
  
  cellchat@meta$condition <- condition_name
  
  CellChatDB <- CellChatDB.mouse
  
  CellChatDB.use <- subsetDB(
    CellChatDB,
    search = "Secreted Signaling",
    key = "annotation"
  )
  
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  
  future::plan("multisession", workers = 4)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(
    object = cellchat,
    type = "triMean"
  )
  
  cellchat <- filterCommunication(
    object = cellchat,
    min.cells = 10
  )
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  communication_table <- subsetCommunication(cellchat)
  
  write.csv(
    communication_table,
    file = file.path(
      output_dir,
      "communication_tables",
      paste0(condition_name, "_subcluster_LR_communication.csv")
    ),
    row.names = FALSE
  )
  
  pathway_table <- subsetCommunication(
    object = cellchat,
    slot.name = "netP"
  )
  
  write.csv(
    pathway_table,
    file = file.path(
      output_dir,
      "pathway_tables",
      paste0(condition_name, "_subcluster_pathway_communication.csv")
    ),
    row.names = FALSE
  )
  
  return(cellchat)
}

# -----------------------------
# 6. Split object by condition
# -----------------------------

seurat_Chow <- subset(
  x = seurat_major,
  subset = condition == "Chow"
)

seurat_HFSC <- subset(
  x = seurat_major,
  subset = condition == "HFSC"
)

# -----------------------------
# 7. Run CellChat for each condition
# -----------------------------

cellchat_Chow <- run_cellchat(
  seurat_object = seurat_Chow,
  group_by = "cellchat_group",
  condition_name = "Chow",
  output_dir = output_dir
)

cellchat_HFSC <- run_cellchat(
  seurat_object = seurat_HFSC,
  group_by = "cellchat_group",
  condition_name = "HFSC",
  output_dir = output_dir
)

# -----------------------------
# 8. Save individual CellChat objects
# -----------------------------

saveRDS(
  cellchat_Chow,
  file = file.path(output_dir, "cellchat_Chow_subcluster.rds")
)

saveRDS(
  cellchat_HFSC,
  file = file.path(output_dir, "cellchat_HFSC_subcluster.rds")
)

# -----------------------------
# 9. Merge CellChat objects for comparison
# -----------------------------

cellchat_list <- list(
  Chow = cellchat_Chow,
  HFSC = cellchat_HFSC
)

cellchat_merged <- mergeCellChat(
  object.list = cellchat_list,
  add.names = names(cellchat_list)
)

saveRDS(
  cellchat_merged,
  file = file.path(output_dir, "cellchat_merged_subcluster.rds")
)

# -----------------------------
# 10. Compare global interaction number and strength
# -----------------------------

pdf(file.path(output_dir, "plots", "subcluster_compare_interaction_number.pdf"), width = 5, height = 4)
compareInteractions(
  object = cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
)
dev.off()

pdf(file.path(output_dir, "plots", "subcluster_compare_interaction_strength.pdf"), width = 5, height = 4)
compareInteractions(
  object = cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)
dev.off()

# -----------------------------
# 11. Compare pathway information flow
# -----------------------------

pdf(file.path(output_dir, "plots", "subcluster_rankNet_comparison.pdf"), width = 8, height = 6)
rankNet(
  object = cellchat_merged,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)
dev.off()

# -----------------------------
# 12. Export selected pathway interactions
# -----------------------------

pathways_selected <- c(
  "COMPLEMENT",
  "CCL",
  "ANGPTL",
  "MIF",
  "TGFb"
)

for (pathway in pathways_selected) {
  
  message("Exporting pathway: ", pathway)
  
  chow_pathway <- subsetCommunication(
    object = cellchat_Chow,
    signaling = pathway
  )
  
  hfsc_pathway <- subsetCommunication(
    object = cellchat_HFSC,
    signaling = pathway
  )
  
  write.csv(
    chow_pathway,
    file = file.path(
      output_dir,
      "pathway_tables",
      paste0("Chow_", pathway, "_subcluster_LR_communication.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    hfsc_pathway,
    file = file.path(
      output_dir,
      "pathway_tables",
      paste0("HFSC_", pathway, "_subcluster_LR_communication.csv")
    ),
    row.names = FALSE
  )
}

# -----------------------------
# 13. Minimal pathway visualization
# -----------------------------

# These plots are workflow records only.
# Final manuscript figures were customized separately.

for (pathway in pathways_selected) {
  
  if (pathway %in% cellchat_Chow@netP$pathways) {
    pdf(
      file.path(output_dir, "plots", paste0("Chow_", pathway, "_subcluster_heatmap.pdf")),
      width = 6,
      height = 5
    )
    print(netVisual_heatmap(cellchat_Chow, signaling = pathway))
    dev.off()
  }
  
  if (pathway %in% cellchat_HFSC@netP$pathways) {
    pdf(
      file.path(output_dir, "plots", paste0("HFSC_", pathway, "_subcluster_heatmap.pdf")),
      width = 6,
      height = 5
    )
    print(netVisual_heatmap(cellchat_HFSC, signaling = pathway))
    dev.off()
  }
}

# -----------------------------
# 14. Session information
# -----------------------------

sink(file.path(output_dir, "session_info_cellchat_subcluster.txt"))
sessionInfo()
sink()