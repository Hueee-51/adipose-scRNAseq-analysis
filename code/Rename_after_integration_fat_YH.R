rm(list = ls())
setwd('/Users/hueee/Lab/Feinberg\ Lab/Single_cell_adipose_tissue')
#install.packages("tidyverse")
##setwd("~/Downloads/scRNA_August")
##knitr::purl('DE_EC.Rmd')
#BiocManager::install('scCustomize')
#BiocManager::install('qs')
#BiocManager::install("ComplexHeatmap")
library(tidyverse)
library(Matrix)
library(RCurl)
library(cowplot)
library(Seurat)
library(multtest)
library(readr)
library(metap)
library(SingleCellExperiment)
library(AnnotationHub)
library(ensembldb)
library(patchwork)
library(qs)
library(scCustomize)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggsci)

output_dir = file.path('./3_Integration')
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Path already exists!")
}


seurat.combined.sct = readRDS('./3_Integration/data/2024-01-12/integration/seurat_combined_sct.rds', refhook = NULL)
seurat.combined.sct[[]] %>% head()
Idents(seurat.combined.sct) <- 'integrated_snn_res.0.5'

png(filename = './3_Integration/combined_cluster_annotation/umap_0502.png')
DimPlot(seurat.combined.sct, reduction = "umap", label = T)
dev.off()
# markers.to.plot <- c("Fbn1", "Msln", "Cdh5", "Steap4", "Myh11", "Mafb",  "Flt3", "Csf3r", "Ms4a1", "Cd3e", "Myocd","Lyz1","Sirpb1c","Acta2")

markers.to.plot <- c("Pecam1","Ptprc","Ter119","Itgb1","Cd34","Ly6a","Cfd") # ASPC:CD31-/CD45-/TER119-/CD29+/CD34+/SCA1+
markers.to.plot <- c("Tcf21","Pdgfra","Icam1","Dpp4","Aqp1") # PECAM-1 (CD31)
markers.to.plot <- c("Kcnj8", "Abcc9", "Vtn","Cspg4","Pdgfrb")
markers.to.plot <-"Acta2"


for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
    cat("### ",gene, "\n")
    png(filename = sprintf('./3_Integration/markers_for_clustering/markers_%s.png', gene), width=8,height=8,units="in",res=150)
    p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap",
                     order = TRUE,
                     min.cutoff = 'q1', 
                     label = TRUE)
    #split.by = "orig.ident")
    print(p)
    dev.off()
    cat("\n")
  }
  else{
    cat("### ",gene, " \n")
    cat("Not present in dataset \n")
  }
  cat("\n")
}

DT::datatable(as.matrix(table(seurat.combined.sct$seurat_clusters)), colnames = c("cluster", "number of cells"))

DimPlot(seurat.combined.sct, reduction = "umap", label =T)

# Renamed:
seurat.combined.sct <- NormalizeData(seurat.combined.sct, assay = "RNA")
DefaultAssay(seurat.combined.sct) <- "RNA"

assignment_df <- data.frame(
  cluster = c("10","15","0",
              "13",
              "1","2","9","11",
              "5","6","7","8","12","14","17",
              "16",
              "19",
              "4","18","20",
              "3"), 
  renamed_cluster = c(
                      rep("ASPC",3), 
                      rep("Mesothelial Cell",1),
                      rep("Endothelial Cell",4),
                      rep("Macrophage",7),
                      rep("Dendritic Cell", 1), 
                      rep("Neutrophil", 1), 
                      rep("B Cell",3),
                      rep("T Cell",1)
                    )
            )

print(assignment_df)



assignment_df$renamed_cluster = as.character(assignment_df$renamed_cluster)

#--------------------------------------------------------------------------------
temp_meta <- seurat.combined.sct@meta.data %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "barcode") %>% 
  left_join(assignment_df, by = c("seurat_clusters" = "cluster")) %>% 
  tibble::column_to_rownames(var = "barcode")

temp_meta <- temp_meta[,-33]
colnames(temp_meta)[32] <- "renamed_cluster"
temp_meta$renamed_cluster <- fct_explicit_na(temp_meta$renamed_cluster, "Unknown")

seurat.combined.sct$renamed_cluster <- temp_meta[rownames(seurat.combined.sct@meta.data),"renamed_cluster"]

seurat.combined.sct$seurat_clusters <- as.factor(seurat.combined.sct$seurat_clusters)
seurat.combined.sct$renamed_cluster <- as.factor(seurat.combined.sct$renamed_cluster)

DT::datatable(as.matrix(table(seurat.combined.sct$renamed_cluster)), colnames = c("cell type", "number of cells"))

# Generate a column for sample conditions
#-------------------------------------------------------------------------------------------
seurat.combined.sct$sample_simple <- stringi::stri_sub(seurat.combined.sct$sample,-4,-1)
#seurat.combined.sct@meta.data[seurat.combined.sct@meta.data$sample_simple %in% 'WD',]$sample_simple = 'WT-TCDD'
#seurat.combined.sct@meta.data[seurat.combined.sct@meta.data$sample_simple %in% 'KT',]$sample_simple = 'KO-TCDD'
seurat.combined.sct$DE_group <- paste0(seurat.combined.sct$renamed_cluster, "-", seurat.combined.sct$sample_simple)
seurat.combined.sct$sample_simple <- as.factor(seurat.combined.sct$sample_simple)
#seurat.combined.sct$sample_simple<- relevel(seurat.combined.sct$sample_simple, "WT-TCDD",after = 1)
#seurat.combined.sct$sample_simple<- relevel(seurat.combined.sct$sample_simple, ref = "WT")

DT::datatable(as.matrix(table(seurat.combined.sct$DE_group)), colnames = c("cell type", "number of cells"))

saveRDS(seurat.combined.sct, "./3_Integration/combined/renamed/seurat_combined_sct_umap_renamed.rds")
#seurat_preintegration= readRDS('./3_Integration/seurat_preintegration.rds', refhook = NULL)
seurat.combined.sct = readRDS('./3_Integration/combined/renamed/seurat_combined_sct_umap_renamed.rds', refhook = NULL)

# Plotting the reannotated umap
#-Umap plot--------------------------------------------------------------------------------------------
#mycol <- c( "#d5b491","#E6A819","#7E6148B2","#c9d0c5","#E64B35B2", "#F1CE83","#799385", "#4DBBD5B2",  "#748eb3" )
#mycol <- c("#B4DEE2","#cd6e1b","#a8afbb","#d7c5ad","#e7b24f","#a8af98","#e6d8cf","#866373","#43ae9e")
#mycol <- c("#91D1C2B2", "#4DBBD5B2", "#DC0000B2", "#e7b24f", "#8491B4B2", "#7E6148B2","#00A087B2", "#ad66a9", "#F39B7FB2")
mycol <- c("#8491B4B2", "#e7b24f", "#91D1C2B2", "#4DBBD5B2", "#F39B7FB2",  "#ad66a9", "#DC0000B2", "#00A087B2")
Idents(seurat.combined.sct) <- "renamed_cluster"
levels(seurat.combined.sct) <- c("Macrophage","Endothelial Cell", "ASPC", "B Cell","T Cell","Mesothelial Cell","Dendritic Cell","Neutrophil")
levels(seurat.combined.sct)

png(filename = './3_Integration/combined/renamed/cluter_renamed_label_splitsample.png', width=9,height=5,units="in",res=500)
DimPlot(seurat.combined.sct, label = T, split.by = 'sample_simple', 
        ncol = 2, cols = mycol, pt.size= 0.05, label.color= "#0E2432", label.size= 3) + theme(axis.text = element_text(size = 5),
                                                                                              axis.title = element_text(size = 7))
dev.off()
png(filename = './3_Integration/combined/renamed/cluter_renamed_label_splitsample_small.png', width=6,height=3,units="in",res=500)
DimPlot(seurat.combined.sct, label = T, split.by = 'sample_simple', 
        ncol = 2, cols = mycol, pt.size= 0.01, label.color= "#0E2432", label.size= 1.5) + theme(axis.text = element_text(size = 6),
                                                                                              axis.title = element_text(size = 8))
dev.off()
pdf('./3_Integration/combined/renamed/cluter_renamed_label_splitsamplep.pdf', width=9,height=5)
DimPlot(seurat.combined.sct, label = T, split.by = 'sample_simple', 
        ncol = 2, cols = mycol, pt.size= 0.05, label.color= "#0E2432",label.size= 3) + theme(axis.text = element_text(size = 5),
                                                                                             axis.title = element_text(size = 7))
dev.off()
png(filename = './3_Integration/combined/renamed/cluter_renamed_nolabel_splitsample.png', width=9,height=5,units="in",res=500)
DimPlot(seurat.combined.sct, label = F, split.by = 'sample_simple', 
        ncol = 2,cols = mycol, pt.size= 0.05, label.color= "#0E2432",label.size= 3) + theme(axis.text = element_text(size = 5),
                                                                                            axis.title = element_text(size = 7))
dev.off()
png(filename = './3_Integration/combined/renamed/cluter_renamed_label_umap.png', width=6,height=5,units="in",res=500)
DimPlot(seurat.combined.sct, label = T, cols = mycol, pt.size= 0.05, label.color= "#0E2432",label.size= 3) + 
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10))
dev.off()
png(filename = './3_Integration/combined/renamed/cluter_renamed_label_umap_small.png', width=4.5,height=3,units="in",res=500)
DimPlot(seurat.combined.sct, label = T, cols = mycol, pt.size= 0.01, label.color= "#0E2432",label.size= 1.5) + 
  theme(axis.text = element_text(size = 6),axis.title = element_text(size = 8))
dev.off()
pdf('./3_Integration/combined/renamed/cluter_renamed_label_umap.pdf', width=6,height=5)
DimPlot(seurat.combined.sct, label = T,cols = mycol, pt.size= 0.05, label.color= "#0E2432",label.size= 3) + 
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10))
dev.off()
png(filename = './3_Integration/combined/renamed/cluter_renamed_nolabel_umap.png', width=6,height=5,units="in",res=500)
DimPlot(seurat.combined.sct, label = F,cols = mycol, pt.size= 0.05, label.color= "#0E2432",label.size= 3) + 
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10))
dev.off()

# Highlight cell types
Idents(seurat.combined.sct) <- "renamed_cluster"

png(filename = sprintf('./3_Integration/combined/renamed/highlight macrophage.png'),width=10,height=5,units="in",res=500)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "Macrophage", highlight_color = "#8491B4B2",
                       background_color = "lightgray", split_seurat= T, split.by = "sample_simple",label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("")
dev.off()

png(filename = sprintf('./3_Integration/combined/renamed/highlight EC.png'),width=10,height=5,units="in",res=150)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "Endothelial Cell", highlight_color = "#e7b24f",
                       background_color = "lightgray", split.by = "sample_simple", split_seurat= T, label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("")
dev.off()

png(filename = sprintf('./3_Integration/combined/renamed/highlight ASPC.png'),width=10,height=5,units="in",res=150)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "ASPC", highlight_color = "#91D1C2B2",
                       background_color = "lightgray", split_seurat= T, split.by = "sample_simple",label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("")
dev.off()

#png(filename = sprintf('./3_Integration/combined/renamed/highlight pericyte.png'),width=10,height=5,units="in",res=150)
#Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "Pericyte", highlight_color = '#0d8a8c',
                       #background_color = "lightgray", split.by = "sample_simple", split_seurat= T, label = T,
                       #pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("")
#dev.off()

png(filename = sprintf('./3_Integration/combined/renamed/highlight B cell.png'),width=10,height=5,units="in",res=150)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "B Cell", highlight_color = "#ad66a9",
                       background_color = "lightgray", split.by = "sample_simple", split_seurat= T, label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("")
dev.off()

png(filename = sprintf('./3_Integration/combined/renamed/highlight T cell.png'),width=10,height=5,units="in",res=150)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "T Cell", highlight_color = "#F39B7FB2",
                       background_color = "lightgray", split.by = "sample_simple", split_seurat= T, label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("")
dev.off()




# No split
pdf('./3_Integration/combined/renamed/highlight Macrophage.pdf',width=5,height=5)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "Macrophage", highlight_color = "#8491B4B2",
                       background_color = "lightgray", label = T, 
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("") + theme(legend.position = "none",
                                                                                                  axis.text = element_text(size = 6),
                                                                                                  axis.title.x = element_text(size = 8), 
                                                                                                  axis.title.y = element_text(size = 8))
dev.off()

pdf('./3_Integration/combined/renamed/highlight EC.pdf',width=5,height=5)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "Endothelial Cell", highlight_color = "#e7b24f",
                       background_color = "lightgray", label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("") + theme(legend.position = "none",
                                                                                                  axis.text = element_text(size = 6),
                                                                                                  axis.title.x = element_text(size = 8), 
                                                                                                  axis.title.y = element_text(size = 8))
dev.off()

#pdf('./3_Integration/combined/renamed/highlight pericyte.pdf',width=5,height=5)
#Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "Pericyte", highlight_color = '#0d8a8c',
                       #background_color = "lightgray", label = T,
                       #pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("") + theme(legend.position = "none",
                                                                                                  #axis.text = element_text(size = 6),
                                                                                                  #axis.title.x = element_text(size = 8), 
                                                                                                  #axis.title.y = element_text(size = 8))
#dev.off()

pdf('./3_Integration/combined/renamed/highlight ASPC.pdf',width=5,height=5)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "ASPC", highlight_color = "#91D1C2B2",
                       background_color = "lightgray", label = T, 
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("") + theme(legend.position = "none",
                                                                                                  axis.text = element_text(size = 6),
                                                                                                  axis.title.x = element_text(size = 8), 
                                                                                                  axis.title.y = element_text(size = 8))
dev.off()

pdf('./3_Integration/combined/renamed/highlight B cell.pdf',width=5,height=5)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "B Cell", highlight_color = "#4DBBD5B2",
                       background_color = "lightgray", label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("") + theme(legend.position = "none",
                                                                                                  axis.text = element_text(size = 6),
                                                                                                  axis.title.x = element_text(size = 8), 
                                                                                                  axis.title.y = element_text(size = 8))
dev.off()

pdf('./3_Integration/combined/renamed/highlight T cell.pdf',width=5,height=5)
Cluster_Highlight_Plot(seurat_object = seurat.combined.sct, cluster_name = "T Cell", highlight_color = "#F39B7FB2",
                       background_color = "lightgray", label = T,
                       pt.size= 0.05, label.color= "#0E2432",label.size= 3) + ggtitle("") + theme(legend.position = "none",
                                                                                                  axis.text = element_text(size = 6),
                                                                                                  axis.title.x = element_text(size = 8), 
                                                                                                  axis.title.y = element_text(size = 8))
dev.off()


#-dot plot------------------------------------------------------------------------------------------------
seurat.combined.sct <- NormalizeData(seurat.combined.sct, assay = "RNA")
DefaultAssay(seurat.combined.sct) <- "RNA"

markers.to.plot <- c("Cd3e","Csf3r","Msln","Mafb","Cdh5","Flt3","Ms4a1","Fbn1")  

Idents(seurat.combined.sct) <- "renamed_cluster"
png(file=sprintf('./3_Integration/combined/renamed/Dot_plot_combined.png'),width=8,height=5,units="in",res=500)
p2 <- DotPlot(seurat.combined.sct, 
              features = markers.to.plot, 
              dot.scale = 8,
              cols = "RdYlBu") +
  labs(x = "Markers", y = "Cell Types", title = "")
              # split.by = "sample_simple") + 
  RotatedAxis()
print(p2)
cat("\n")
dev.off()


#-------marker------------
#markers = read.csv('./subset_FBS/newcutoff/FBS_marker_0.3.csv')
#top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#top10 <- markers %>% group_by(cluster) %>% top_n(n = -25,wt = p_val_adj)
DefaultAssay(seurat.combined.sct) <- "RNA" 
Idents(seurat.combined.sct) <- "renamed_cluster"
levels(seurat.combined.sct) <- c("Macrophage","Endothelial Cell", "ASPC", "B Cell","T Cell","Mesothelial Cell","Dendritic Cell","Neutrophil")
levels(seurat.combined.sct)
markers <- FindAllMarkers(seurat.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,file = sprintf('./3_Integration/combined/renamed/cell_marker.csv'),row.names = T)
markers = read.csv('./3_Integration/combined/renamed/cell_marker.csv')
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#top10$cluster <- factor(top10$cluster, levels = c("Macrophage", "Endothelial Cell", "ASPC", "B Cell", "T Cell",  "Mesothelial Cell", "Dendritic Cell", "Neutrophil"))
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20, file = sprintf('./3_Integration/combined/renamed/cell_marker_top20.csv'),row.names = T)

# top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = p_val)
#write.csv(top10,file = sprintf('./subset_FBS/newcutoff/FBS_marker_top10.csv'),row.names = T)

DefaultAssay(seurat.combined.sct) <- "integrated"
seurat.combined.sct <- NormalizeData(seurat.combined.sct, assay = "RNA")
seurat.combined.sct <- ScaleData(seurat.combined.sct)

pdf('./3_Integration/combined/renamed/featureheatmap_avglog2FC.pdf',width=10,height=6)
DoHeatmap(seurat.combined.sct, features = top10_sorted$gene, group.colors = mycol)
dev.off()

#png(filename = sprintf('./3_Integration/combined/renamed/featureheatmap_pval.png'),width=12,height=10,units="in",res=500)
#DoHeatmap(seurat.combined.sct, features = top10$gene, group.colors = mycol) #+ NoLegend()
#dev.off()

#For a more advanced take on the DotPlot scCustomize contains a separate function Clustered_DotPlot(). 
#This function allows for clustering of both gene expression patterns and identities in the final plot. 
library(scCustomize)
markers <- FindAllMarkers(object = seurat.combined.sct,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- Extract_Top_Markers(marker_dataframe = markers, num_genes = 5, named_vector = FALSE,
                                    make_unique = TRUE)


png(filename = sprintf('./3_Integration/combined/renamed/top5_marker.png'),width=8,height=12,units="in",res=150)
Clustered_DotPlot(seurat_object = seurat.combined.sct, features = top5_markers)
dev.off()
png(filename = sprintf('./3_Integration/combined/renamed/top5_marker_k means_number5.png'),width=8,height=12,units="in",res=150)
Clustered_DotPlot(seurat_object = seurat.combined.sct, features = top5_markers,k = 5)
dev.off()
png(filename = sprintf('./3_Integration/combined/renamed/top5_marker_k means_number5-1.png'),width=8,height=10,units="in",res=150)
Clustered_DotPlot(seurat_object = seurat.combined.sct, features = top5_markers,k = 5)
dev.off()
png(filename = sprintf('./3_Integration/combined/renamed/top5_marker_k means_1_number5.png'),width=8,height=12,units="in",res=150)
#png(filename = sprintf('./subset_FBS/newcutoff/renamed/highlight cluster8.png'),width=8,height=8,units="in",res=150)
Clustered_DotPlot(seurat_object = seurat.combined.sct, features = top5_markers, k = 5, exp_color_min = -1, exp_color_max = 2,
                  colors_use_exp = PurpleAndYellow())
dev.off()


#######################
### marker heatmap umap
library(patchwork)

# Marcophages
markers.to.plot <- c("Adgre1", "Cd68", "Mafb", "Trem2", 
                    "Csf1r", "Mrc1","C1qc","Folr2"
                     #"Cd163","Cbr2", "Folr2",
                     )
# EC
markers.to.plot <- c("Pecam1", "Tie1", "Cdh5", "Nos3",
                     "Adgrl4","Gpihbp1", "Cyyr1", "Sox17"
                     # "Klf2"
                     )

# ASPC
markers.to.plot <- c("Pdgfra","Pecam1","Ptprc",
                     "Itgb1","Cd34","Ly6a")

# B cell
markers.to.plot <- c("Cd79a","Ly6d","Cd19","Cd79b")

# T cell
markers.to.plot <- c("Cd3e","Cd4","Cd8a","Trbc2")

# Pericyte
# markers.to.plot <- c("Eng", "Nt5e", "Des","Itgam","Ptprc","Vtn", "Pecam1", "Nes","Dpp4", "Pdgfra",
                     # "Mcam", "Rgs5", "Pdgfrb", "Cspg4", "Ccl2", "Abcc9", "Notch3",
                     # "Steap4") # Desmin:Des; 


Idents(seurat.combined.sct) <- "renamed_cluster"

plots_list <- list()
for(gene in markers.to.plot) {
  if(gene %in% rownames(seurat.combined.sct)) {
    cat("### ", gene, "\n")
    png(filename = sprintf('./3_Integration/combined/renamed/markers/DEG_T_%s.png', gene), width = 3, height = 3, units = "in", res = 500)
    p <- FeaturePlot(seurat.combined.sct, 
                     features = gene, reduction = "umap",
                     order = TRUE,
                     label = TRUE,
                     max.cutoff = 50,
                     pt.size = 0.01,
                     label.size = 1.5,
                     min.cutoff = "q10", 
                     # cols = c("#e9ecf2","#8ea0cc"))+# Macrophage
                     # cols = c("#e9ecf2","#e7b24f"))+ # Endothelial Cell
                     # cols = c("#e9ecf2","#91D1C2B2"))+ # ASPC
                     # cols = c("#e9ecf2","#4DBBD5B2"))+ # B Cell
                      cols = c("#e9ecf2","#F39B7FB2"))+ # T Cell
                     # cols = c("#e9ecf2","#0d8a8c"))+ # Pericyte
      theme(axis.text = element_text(size = 6),
            axis.title = element_text(size = 6),
            axis.line = element_line(linewidth = 0.5),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2, 'mm'),
            title = element_text(size = 11))
    print(p)
    dev.off()
    plots_list[[gene]] <- p # 将图添加到列表中
    cat("\n")
  } else {
    cat("### ", gene, " \n")
    cat("Not present in dataset \n")
  }
  cat("\n")
}

final_plot <- wrap_plots(plots_list, ncol = 4)

# 使用 patchwork 拼接图形
pdf('./3_Integration/combined/renamed/markers/DEG_t_Merge.pdf', width = 12, height = 3)
print(final_plot)
dev.off()
png('./3_Integration/combined/renamed/markers/DEG_t_Merge.png', width = 12, height = 3, units = "in", res = 500)
print(final_plot)
dev.off()


# Violin plot--------------------------------------------------------------------------------------------
png(file=sprintf('./3_Integration/combined/Violin_plot_Ms4a1 .png'),width=20,height=7,units="in",res=300)
VlnPlot(seurat.combined.sct,
        features= "Ms4a1",
        group.by = "renamed_cluster",
        split.plot = FALSE,
        split.by = "sample_simple")
dev.off()