rm(list = ls())
setwd('/Users/hueee/Lab/Feinberg\ Lab/Single_cell_adipose_tissue/')
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
list(filtered_seurat = "./data/seurat_HFSC_filtered.rds", 
    sample_name = "HFSC", output_dir = file.path('./2_Clustering/HFSC'))
## ---------------------------------------------------------------------------------------------------------------
# Load the dataset
seurat <- readRDS(params$filtered_seurat)


## ---------------------------------------------------------------------------------------------------------------
seurat <- SCTransform(seurat, verbose = FALSE)
seurat <- RunPCA(seurat)


## ---------------------------------------------------------------------------------------------------------------
png(filename = './2_Clustering/HFSC/fig1.png')
Seurat::ElbowPlot(seurat, ndims = 30)
dev.off()

## ---------------------------------------------------------------------------------------------------------------
n_pcs <- 20


## ---------------------------------------------------------------------------------------------------------------
seurat <- RunUMAP(seurat, dims = 1:n_pcs)


## ---- results = "hide"------------------------------------------------------------------------------------------
seurat <- FindNeighbors(seurat, dims = 1:n_pcs)
seurat <- FindClusters(seurat, resolution = seq(0,1,0.1))


## ---- results = 'asis'------------------------------------------------------------------------------------------
for(res in seq(0.1,1,0.1)){
  cat("#### res:",res, "\n")
  Idents(seurat) <-  paste0("SCT_snn_res.", res)
  png(filename = sprintf('./2_Clustering/HFSC/res_%s.png', res))
  p2 <- DimPlot(seurat, label = T)
  print(p2)
  dev.off()
  cat("\n")
  cat("\n")

}


## ---------------------------------------------------------------------------------------------------------------
png(filename = './2_Clustering/HFSC/fig2.png')
FeaturePlot(seurat, features = "percent.mt", min.cutoff = 'q10')
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './2_Clustering/HFSC/fig3.png')
FeaturePlot(seurat, features = "nFeature_RNA", min.cutoff = 'q10')
dev.off()

## ---------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat) <- "SCT"


## ---------------------------------------------------------------------------------------------------------------
markers.to.plot <- c("Adipoq","Pdgfra","Msln","Jam2","Cdh5","Pecam1","Prox1","Steap4","Myocd",
                     "Mafb","Cybb","Flt3","Cpa3","Csf3r","Ms4a1","Klrd1","Il7r","Dcdc2a","Erbb4",
                     "Vcam1", "Fn1", "Thbs1", "Col4a3", "Icam1" ,"Col1a1", "Col1a2","Col3a1",
                     "Pdgfra", "Cnn1", "Myh11", "Cdh5", "Cd68", "Cd115", 
                     "Lyz2", "Cd3d", "Cd19", "H2-Ab1", "Cd68", "Hbb-bt",
                     "Hbb-bs", "Gpihbp1", "Lpl", "Scarb1", "Mylip", "Angptl4", "Apoc1", "Pltp",
                     "Bgn","Vcam1", "Mfap5", "Bmp4", "Fn1",
                     "Lox", "Pcolce2", "Fbln5", "Efemp1", "Thbs1", "Col4a3",
                     "Col5a2", "Icam1", "Col8a1", "Dcn")

## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./2_Clustering/HFSC/markers_%s.png', gene))  
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
saveRDS(seurat, "./2_Clustering/HFSC/seurat_sctransformed_HFSC.rds")
#saveRDS(seurat, file = file.path(data_dir, paste0("seurat_sctransformed_", params$sample_name, ".rds")))

