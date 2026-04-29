# ============================================================
# Script: 08_cellchat_major_celltype_analysis.R
# Project: Adipose scRNA-seq analysis
# Purpose: CellChat analysis at the major cell type level
# Author: Yi Hu
# ============================================================

# This script performs CellChat analysis using annotated major cell types.
# Chow and HFSC samples are analyzed separately and then compared.
#
# Input:
#   results/annotation/seurat_integrated_annotated.rds
#
# Output:
#   results/cellchat_major/cellchat_Chow_major.rds
#   results/cellchat_major/cellchat_HFSC_major.rds
#   results/cellchat_major/cellchat_merged_major.rds
#   results/cellchat_major/communication_tables/
#   results/cellchat_major/pathway_tables/

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

input_file <- "results/annotation/seurat_integrated_annotated.rds"
output_dir <- "results/cellchat_major"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "communication_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "pathway_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3. Load annotated Seurat object
# -----------------------------

seurat_obj <- readRDS(input_file)

# Major cell type annotation should be generated in:
# scripts/03_clustering_annotation.R

Idents(seurat_obj) <- "cell_type"

# -----------------------------
# 4. Define CellChat helper function
# -----------------------------

run_cellchat <- function(seurat_object, group_by, condition_name, output_dir) {
  
  message("Running CellChat for: ", condition_name)
  
  # Extract normalized or count matrix.
  # CellChat commonly uses normalized expression data, but count data can
  # also be used depending on workflow. Here we use RNA counts to match
  # the original analysis skeleton.
  
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
  
  # Create CellChat object.
  cellchat <- createCellChat(
    object = data_input,
    meta = meta_cellchat,
    group.by = "celltype"
  )
  
  cellchat@meta$condition <- condition_name
  
  # Set CellChat mouse ligand-receptor database.
  CellChatDB <- CellChatDB.mouse
  
  # Use secreted signaling interactions.
  CellChatDB.use <- subsetDB(
    CellChatDB,
    search = "Secreted Signaling",
    key = "annotation"
  )
  
  cellchat@DB <- CellChatDB.use
  
  # Preprocess expression data.
  cellchat <- subsetData(cellchat)
  
  future::plan("multisession", workers = 4)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Infer communication probability.
  cellchat <- computeCommunProb(
    object = cellchat,
    type = "triMean"
  )
  
  # Remove cell-cell interactions with too few cells.
  cellchat <- filterCommunication(
    object = cellchat,
    min.cells = 10
  )
  
  # Infer communication at the pathway level.
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate communication network.
  cellchat <- aggregateNet(cellchat)
  
  # Export ligand-receptor level communication table.
  communication_table <- subsetCommunication(cellchat)
  
  write.csv(
    communication_table,
    file = file.path(
      output_dir,
      "communication_tables",
      paste0(condition_name, "_major_celltype_LR_communication.csv")
    ),
    row.names = FALSE
  )
  
  # Export pathway-level communication table.
  pathway_table <- subsetCommunication(
    object = cellchat,
    slot.name = "netP"
  )
  
  write.csv(
    pathway_table,
    file = file.path(
      output_dir,
      "pathway_tables",
      paste0(condition_name, "_major_celltype_pathway_communication.csv")
    ),
    row.names = FALSE
  )
  
  return(cellchat)
}

# -----------------------------
# 5. Split object by condition
# -----------------------------

seurat_Chow <- subset(
  x = seurat_obj,
  subset = condition == "Chow"
)

seurat_HFSC <- subset(
  x = seurat_obj,
  subset = condition == "HFSC"
)

# -----------------------------
# 6. Run CellChat for each condition
# -----------------------------

cellchat_Chow <- run_cellchat(
  seurat_object = seurat_Chow,
  group_by = "cell_type",
  condition_name = "Chow",
  output_dir = output_dir
)

cellchat_HFSC <- run_cellchat(
  seurat_object = seurat_HFSC,
  group_by = "cell_type",
  condition_name = "HFSC",
  output_dir = output_dir
)

# -----------------------------
# 7. Save individual CellChat objects
# -----------------------------

saveRDS(
  cellchat_Chow,
  file = file.path(output_dir, "cellchat_Chow_major.rds")
)

saveRDS(
  cellchat_HFSC,
  file = file.path(output_dir, "cellchat_HFSC_major.rds")
)

# -----------------------------
# 8. Merge CellChat objects for comparison
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
  file = file.path(output_dir, "cellchat_merged_major.rds")
)

# -----------------------------
# 9. Compare global interaction number and strength
# -----------------------------

pdf(file.path(output_dir, "plots", "major_compare_interaction_number.pdf"), width = 5, height = 4)
compareInteractions(
  object = cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
)
dev.off()

pdf(file.path(output_dir, "plots", "major_compare_interaction_strength.pdf"), width = 5, height = 4)
compareInteractions(
  object = cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)
dev.off()

# -----------------------------
# 10. Compare pathway information flow
# -----------------------------

pdf(file.path(output_dir, "plots", "major_rankNet_comparison.pdf"), width = 8, height = 6)
rankNet(
  object = cellchat_merged,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)
dev.off()

# -----------------------------
# 11. Export selected pathway interactions
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
      paste0("Chow_", pathway, "_major_LR_communication.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    hfsc_pathway,
    file = file.path(
      output_dir,
      "pathway_tables",
      paste0("HFSC_", pathway, "_major_LR_communication.csv")
    ),
    row.names = FALSE
  )
}

# -----------------------------
# 12. Minimal pathway visualization
# -----------------------------

# These plots are included only as workflow records.
# Publication-quality CellChat figures were customized separately.

for (pathway in pathways_selected) {
  
  if (pathway %in% cellchat_Chow@netP$pathways) {
    pdf(
      file.path(output_dir, "plots", paste0("Chow_", pathway, "_major_heatmap.pdf")),
      width = 5,
      height = 4
    )
    print(netVisual_heatmap(cellchat_Chow, signaling = pathway))
    dev.off()
  }
  
  if (pathway %in% cellchat_HFSC@netP$pathways) {
    pdf(
      file.path(output_dir, "plots", paste0("HFSC_", pathway, "_major_heatmap.pdf")),
      width = 5,
      height = 4
    )
    print(netVisual_heatmap(cellchat_HFSC, signaling = pathway))
    dev.off()
  }
}

# -----------------------------
# 13. Session information
# -----------------------------

sink(file.path(output_dir, "session_info_cellchat_major.txt"))
sessionInfo()
sink()