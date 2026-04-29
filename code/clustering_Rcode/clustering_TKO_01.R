rm(list = ls())
setwd('/Users/jc1012/Downloads/scRNA_August/')
#knitr::purl('clustering_s181_03_WT.Rmd')
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

params <-
list(filtered_seurat = "./data/seurat_TKO_01_filtered.rds", 
    sample_name = "TKO_01", output_dir = file.path('./clustering/TKO1'))
## ---------------------------------------------------------------------------------------------------------------
# Load the dataset
seurat <- readRDS(params$filtered_seurat)


## ---------------------------------------------------------------------------------------------------------------
seurat <- SCTransform(seurat, verbose = FALSE)
seurat <- RunPCA(seurat)


## ---------------------------------------------------------------------------------------------------------------
png(filename = './clustering/TKO1/fig1.png')
Seurat::ElbowPlot(seurat, ndims = 30)
dev.off()

## ---------------------------------------------------------------------------------------------------------------
n_pcs <- 25


## ---------------------------------------------------------------------------------------------------------------
seurat <- RunUMAP(seurat, dims = 1:n_pcs)


## ---- results = "hide"------------------------------------------------------------------------------------------
seurat <- FindNeighbors(seurat, dims = 1:n_pcs)
seurat <- FindClusters(seurat, resolution = seq(0,1,0.1))


## ---- results = 'asis'------------------------------------------------------------------------------------------
for(res in seq(0.1,1,0.1)){
  cat("#### res:",res, "\n")
  Idents(seurat) <-  paste0("SCT_snn_res.", res)
  png(filename = sprintf('./clustering/TKO1/res_%s.png', res))
  p2 <- DimPlot(seurat, label = T)
  print(p2)
  dev.off()
  cat("\n")
  cat("\n")

}


## ---------------------------------------------------------------------------------------------------------------
png(filename = './clustering/TKO1/fig2.png')
FeaturePlot(seurat, features = "percent.mt", min.cutoff = 'q10')
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './clustering/TKO1/fig3.png')
FeaturePlot(seurat, features = "nFeature_RNA", min.cutoff = 'q10')
dev.off()

## ---------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat) <- "SCT"


## ---------------------------------------------------------------------------------------------------------------
markers.to.plot <- c("Pi16", "Clec3b", "Ccl11", "C3", "Gsn", "Dcn", "Serping1", "Dpep1",
"Dpt", "Serpinf1", "Pcolce", "Igfbp6", "Scara5", "Tmem100", "Col14a1", "Fbln1",
"C1s1", "Htra3", "Nid1","Entpd2")


## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./clustering/TKO1/markers_%s.png', gene))  
    p <- FeaturePlot(seurat, features = gene, reduction = "umap")
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



## ---------------------------------------------------------------------------------------------------------------
markers.to.plot <- c("Bgn","Vcam1", "Mfap5", "Bmp4", "Fn1",
                     "Lox", "Pcolce2", "Fbln5", "Efemp1", "Thbs1", "Col4a3",
                     "Col5a2", "Icam1", "Col8a1", "Dcn")


## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./clustering/TKO1/markers_%s.png', gene))  
      p <- FeaturePlot(seurat, features = gene, reduction = "umap")
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



## ---------------------------------------------------------------------------------------------------------------
markers.to.plot <- c("Vcam1", "Fn1", "Thbs1", "Col4a3", "Icam1" ,"Col1a1", "Col1a2","Col3a1")


## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./clustering/TKO1/markers_%s.png', gene))  
     p <- FeaturePlot(seurat, features = gene, reduction = "umap")
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



## ---------------------------------------------------------------------------------------------------------------
markers.to.plot <- c("Gpihbp1", "Lpl", "Scarb1", "Mylip", "Angptl4", "Apoc1", "Pltp")


## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./clustering/TKO1/markers_%s.png', gene))   
     p <- FeaturePlot(seurat, features = gene, reduction = "umap")
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



## ---------------------------------------------------------------------------------------------------------------
cell_type_df <- data.frame(
  cell_type = c("Fibroblasts", rep("VSMC",2), "EC", rep("Monocyte", 3),
                rep("T/B_cells",2), rep("Macrophage/DC",2),
                rep("RBC",2)),
  gene = c("Pdgfra", "Cnn1", "Myh11", "Cdh5", "Cd68", "Cd115", 
                     "Lyz2", "Cd3d", "Cd19", "H2-Ab1", "Cd68", "Hbb-bt",
                     "Hbb-bs")
)


## ---------------------------------------------------------------------------------------------------------------
DT::datatable(cell_type_df %>% 
    group_by(cell_type) %>% 
    summarise(gene = paste(unique(gene), collapse = ', ')))



## ---- results = 'asis'------------------------------------------------------------------------------------------
for(cell_t in unique(cell_type_df$cell_type)){
  cat("### ",cell_t,"{.tabset} \n")
  df <- cell_type_df %>% dplyr::filter(cell_type == cell_t )
  
  for( gene in df$gene){
    if(gene %in% rownames(seurat)){
    cat("#### ",gene, "\n")
    #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
     #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
      png(filename = sprintf('./clustering/TKO1/markers_%s.png', gene))  
        p <- FeaturePlot(seurat, features = gene, reduction = "umap")
    print(p)
    dev.off()
    cat("\n")
    }
    else{
      cat("#### ",gene, " \n")
      cat("Not present in dataset \n")
    }
    cat("\n")
  }
  cat("\n")
}


## ---------------------------------------------------------------------------------------------------------------
saveRDS(seurat, "./clustering/TKO1/seurat_sctransformed_TKO1.rds")
#saveRDS(seurat, file = file.path(data_dir, paste0("seurat_sctransformed_", params$sample_name, ".rds")))

