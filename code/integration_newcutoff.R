rm(list = ls())
setwd('/Users/hueee/Lab/Feinberg\ Lab/Single_cell_adipose_tissue')

output_dir = file.path('./3_Integration')
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Path already exists!")
}

#setwd('/home/jic621/maomao/single_cell_August')
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
list(seurat_Chow_file = "./2_Clustering/Chow/seurat_sctransformed_Chow.rds", 
    seurat_HFSC_file = "./2_Clustering/HFSC/seurat_sctransformed_HFSC.rds", 
    sample_name = "integration", output_dir = file.path('./3_Integration'))

## ---------------------------------------------------------------------------------------------------------------
# Directory paths
data_dir <- file.path(params$output_dir, "data", Sys.Date(), params$sample_name)

invisible(mapply(
    FUN = dir.create,
    list(data_dir),
    MoreArgs = list(recursive = TRUE, showWarnings = FALSE)
))

## ---------------------------------------------------------------------------------------------------------------
# Load the datasets
seurat_Chow <- readRDS(params$seurat_Chow_file)
seurat_Chow$sample <- "Chow"
#seurat_Chow@meta.data = seurat_Chow@meta.data[!duplicated(rownames(seurat_Chow@meta.data)),]

seurat_HFSC <- readRDS(params$seurat_HFSC_file)
seurat_HFSC$sample <- "HFSC"
#seurat_HFSC@meta.data = seurat_HFSC@meta.data[!duplicated(rownames(seurat_HFSC@meta.data)),]

seurat <- list(seurat_Chow, seurat_HFSC)


## ---------------------------------------------------------------------------------------------------------------
seurat_Chow_RNA <- seurat_Chow
DefaultAssay(seurat_Chow_RNA) <- "RNA"
seurat_Chow_RNA <- DietSeurat(seurat_Chow_RNA, assays = "RNA")

seurat_HFSC_RNA <- seurat_HFSC
DefaultAssay(seurat_HFSC_RNA) <- "RNA"
seurat_HFSC_RNA <- DietSeurat(seurat_HFSC_RNA, assays = "RNA")

seurat_preintegration <- merge(seurat_Chow_RNA, seurat_HFSC_RNA, add.cell.ids = c('_Chow', '_HFSC'))
#seurat_preintegration2 <- merge(seurat_TKO_01_RNA, seurat_TKO_02_RNA, add.cells.ids = c('_tko1', '_tko2'))
#seurat_preintegration <- merge(seurat_preintegration1, seurat_preintegration2)
##Warning message:In CheckDuplicateCellNames(object.list = objects) :
##Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.

## ---------------------------------------------------------------------------------------------------------------
seurat_preintegration <- SCTransform(seurat_preintegration, verbose = F)

## ---------------------------------------------------------------------------------------------------------------
seurat_preintegration <- RunPCA(seurat_preintegration)

#seurat_preintegration= readRDS('./3_Integration/seurat_preintegration.rds', refhook = NULL)
## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/preintergration/fig1.png')
ElbowPlot(seurat_preintegration, ndim = 30)
dev.off()

## ---------------------------------------------------------------------------------------------------------------
n_pcs <- 20

## ---------------------------------------------------------------------------------------------------------------
seurat_preintegration <- RunUMAP(seurat_preintegration, dims = 1:n_pcs)

## ---- results = "hide"------------------------------------------------------------------------------------------
seurat_preintegration <- FindNeighbors(seurat_preintegration, dims = 1:n_pcs)
seurat_preintegration <- FindClusters(seurat_preintegration, resolution = 0.5)

saveRDS(seurat_preintegration, "./3_Integration/seurat_preintegration.rds")
## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/preintergration/fig2.png')
DimPlot(seurat_preintegration, reduction = "umap", label = T)
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/preintergration/fig3.png')
p1 <- DimPlot(seurat_preintegration, reduction = "umap", group.by = "sample",
              shuffle = T)
p1
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/preintergration/fig4.png',width = 800,height = 400)
p1 <- DimPlot(seurat_preintegration, reduction = "umap", split.by = "sample",
              ncol = 2 )
p1
dev.off()

# go to O2 to run the following
## ---------------------------------------------------------------------------------------------------------------
seurat.list <- seurat
if(!file.exists(file.path(data_dir,"seurat_combined_sct.rds"))){
  seurat.list <- lapply(X = seurat.list, FUN = SCTransform, verbose = FALSE) # Applies SCTransform to each Seurat object in `seurat.list` for normalization and variance stabilization.
  features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000) # Selects features (genes) for integration.
  seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features) # Prepares the Seurat objects for integration using the selected features.
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                           normalization.method = "SCT", 
                                           anchor.features = features)  # Identifies anchors for data integration across datasets.
  seurat.combined.sct <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT") # Integrates the data using the identified anchors.
  saveRDS(seurat.combined.sct, file = file.path(data_dir, "seurat_combined_sct.rds")) # Saves the integrated Seurat object to an RDS file for future use.
  }else{
    seurat.combined.sct <- readRDS(file.path(data_dir, "seurat_combined_sct.rds"))
}

#seurat.combined.sct = readRDS('./3_Integration/data/2024-01-12/integration/seurat_combined_sct.rds', refhook = NULL)
#seurat.combined.sct = readRDS('./3_Integration/newqc_combined/2021-10-07/integration/seurat_combined_sct.rds', refhook = NULL)
#seurat.combined.sct <- readRDS(file.path(data_dir, "seurat_combined_sct.rds"))
data_dir <- './3_Integration/data/2024-01-12/integration/'
metadata = seurat.combined.sct@meta.data 
png(filename = './3_Integration/cell_counts.png')
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

#### Victor's suggestion for doing plots
# seurat.combined.sct <- NormalizeData(seurat.combined.sct, assay = "RNA")
# DefaultAssay(seurat.combined.sct) <- "RNA" 
## ---------------------------------------------------------------------------------------------------------------
seurat.combined.sct <- RunPCA(seurat.combined.sct, verbose = FALSE)

## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/fig1.png')
ElbowPlot(seurat.combined.sct, ndim = 30)
dev.off()

## ---------------------------------------------------------------------------------------------------------------
ndims <- 20
seurat.combined.sct <- RunUMAP(seurat.combined.sct, reduction = "pca", dims = 1:ndims)

seurat.combined.sct <- FindNeighbors(seurat.combined.sct, dims = 1:ndims, verbose = FALSE)
seurat.combined.sct <- FindClusters(seurat.combined.sct, verbose = FALSE, resolution = seq(0.1,1,0.1))
seurat.combined.sct <- FindClusters(seurat.combined.sct, verbose = FALSE, resolution = 0.5)
## ---- results = 'asis'------------------------------------------------------------------------------------------
for(res in seq(0.1,1,0.1)){
  cat("#### res:",res, "\n")
  Idents(seurat.combined.sct) <-  paste0("integrated_snn_res.", res)
  png(filename = sprintf('./3_Integration/res_%s.png', res))
  p2 <- DimPlot(seurat.combined.sct, label = T)
  print(p2)
  dev.off()
  cat("\n")
  cat("\n")
}

seurat.combined.sct[[]] %>% head()
Idents(seurat.combined.sct) <- 'integrated_snn_res.0.5'
saveRDS(seurat.combined.sct, file = file.path(data_dir, "seurat_combined_sct.rds"))
## ---------------------------------------------------------------------------------------------------------------

#seurat.combined.sct = readRDS('./3_Integration/newqc_combined/renamed/seurat_combined_sct_umap_FBSmerged.rds')

# Explore heatmap of PCs
DefaultAssay(seurat.combined.sct) <- "SCT"
png(filename = './3_Integration/combined/Heatmap_1.png', width=8,height=8,units="in",res=300)
DimHeatmap(seurat.combined.sct, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()

#---------------------------------------------------------------------
seurat.combined.sct <- PrepSCTFindMarkers(seurat.combined.sct)
aorta.markers <- FindAllMarkers(seurat.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
aorta.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(aorta.markers,file = sprintf('./3_Integration/aorta_marker.csv'),row.names = T)
#aorta.markers = read.csv('../3_Integration/combined/aorta_marker.csv')
top10 <- aorta.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = sprintf('./3_Integration/aorta_marker_top10.csv'),row.names = T)
DefaultAssay(seurat.combined.sct) <- "integrated"
png(filename = sprintf('./3_Integration/featureheatmap1.png'),width=50,height=32,units="in",res=300)
DoHeatmap(seurat.combined.sct, features = top10$gene) #+ NoLegend()
dev.off()

# Printing out the most variable genes driving PCs
print(x = seurat.combined.sct[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
## ---------------------------------------------------------------------------------------------------------------
DT::datatable(as.data.frame.matrix(table(seurat.combined.sct$seurat_clusters, seurat.combined.sct$sample)))

DT::datatable(as.matrix(table(seurat.combined.sct$seurat_clusters)), colnames = c("cluster", "number of cells"))
## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/combined/umap.png')
p2 <- DimPlot(seurat.combined.sct, reduction = "umap", label =T)
p2
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/combined/uamp_sample.png')
p1 <- DimPlot(seurat.combined.sct, reduction = "umap", group.by = "sample", 
              shuffle = T)
p1
dev.off()


## ---------------------------------------------------------------------------------------------------------------
png(filename = './3_Integration/combined/fig4.png', width=800, height=400)
p1 <- DimPlot(seurat.combined.sct, reduction = "umap", split.by = "sample",
              ncol = 2 )
p1
dev.off()


## ---------------------------------------------------------------------------------------------------------------
DT::datatable(as.matrix(table(seurat.combined.sct$seurat_clusters)), colnames = c("cluster", "number of cells"))


## ---------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat.combined.sct) <- "SCT"


## before plotting normalization and change to RNA method---------------------------------------------------------------------------------------------------------------
### seurat.combined.sct <- NormalizeData(seurat.combined.sct, assay = "RNA")
seurat.combined.sct <- NormalizeData(seurat.combined.sct, assay = "RNA")
DefaultAssay(seurat.combined.sct) <- "RNA"
markers.to.plot <- c("Pi16", "Clec3b", "Ccl11", "C3", "Gsn", "Dcn", "Serping1", "Dpep1",
"Dpt", "Serpinf1", "Pcolce", "Igfbp6", "Scara5", "Tmem100", "Col14a1", "Fbln1",
"C1s1", "Htra3", "Nid1","Entpd2",
"Vcam1", "Fn1", "Thbs1", "Col4a3", "Icam1" ,"Col1a1", "Col1a2","Col3a1",
"Pdgfra", "Cnn1", "Myh11", "Cdh5", "Cd68", "Cd115", 
"Lyz2", "Cd3d", "Cd19", "H2-Ab1", "Cd68", "Hbb-bt",
"Hbb-bs", "Gpihbp1", "Lpl", "Scarb1", "Mylip", "Angptl4", "Apoc1", "Pltp",
"Bgn","Vcam1", "Mfap5", "Bmp4", "Fn1",
"Lox", "Pcolce2", "Fbln5", "Efemp1", "Thbs1", "Col4a3",
"Col5a2", "Icam1", "Col8a1", "Dcn",
'AI662270','Snhg12','Snhg18','Bach2os',
"Adipoq","Pdgfra","Msln","Jam2","Cdh5","Pecam1","Prox1","Steap4","Myocd",
"Mafb","Cybb","Flt3","Cpa3","Csf3r","Ms4a1","Klrd1","Il7r","Dcdc2a","Erbb4")

## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./3_Integration/combined/markers_%s.png', gene))
    p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap")
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
  if(gene %in% rownames(seurat.combined.sct)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./3_Integration/combined/markers_%s.png', gene))
      p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap")
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
markers.to.plot <- c("Vcam1", "Fn1", "Thbs1", "Col4a3", "Icam1" )


## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./3_Integration/combined/markers_%s.png', gene))
      p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap")
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
  if(gene %in% rownames(seurat.combined.sct)){
  cat("### ",gene, "\n")
  #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
   #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./3_Integration/combined/markers_%s.png', gene))
      p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap")
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
    if(gene %in% rownames(seurat.combined.sct)){
    cat("#### ",gene, "\n")
    #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
     #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
      png(filename = sprintf('./3_Integration/newqc_combined/markers_%s.png', gene))
        p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap", 
                         label = TRUE,order = TRUE,
                         min.cutoff = 'q10')
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
## ---- fibroblast activation------------------------------------------------------------------------------------------
seurat.combined.sct <- NormalizeData(seurat.combined.sct, assay = "RNA")
DefaultAssay(seurat.combined.sct) <- "RNA"

markers.to.plot <- c("Pdgfra","Ltbp2", "Col5a2", "Sparc", "Col1a2", 
                     "Col1a1", "Lgals1", "Timp1", "Tnc", "Fn1", 
                     "Tpm1", "Eln", "Thbs1", "Col3a1", "Fbn1", "Fstl1", "Acta2", "Col5a1", 
                     "Col12a1", "Wisp1", "Col4a2", "Serpinh1", "Ppic", "Actn1", "Tagln", "Col6a3",
                     "Fmod","Angtl1","Col8a1","Mmp2","Mmp3","Vim","Snhg18","Spp1")


## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
    cat("### ",gene, "\n")
    #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
    #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./3_Integration/newqc_combined/fibroblast_activation/markers_%s.png', gene), width=8,height=8,units="in",res=150)
    p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap",
                     order = TRUE,
                     min.cutoff = 'q10', 
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

## ---- T and NK cells------------------------------------------------------------------------------------------
markers.to.plot <- c("HBB-bs", "Hba-a1", "Hba-a2", "Calca", "Hbb-bt", "Tnfrsf4", "Cd74",
"Ccl1", "Scgb1a1", "Capg", "Rps27", "Ifi27l2a", "Klf2", "Itm2a", "Tnfrsf18","Klf4",
"Klf10", "Tnfrsf9", "Txnip", "S100a4", "Cirbp", "Pdcd1")
## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
    cat("### ",gene, "\n")
    #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
    #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./3_Integration/newqc_combined/T_NK cell/markers_%s.png', gene), width=8,height=8,units="in",res=150)
    p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap",
                     order = TRUE,
                     min.cutoff = 'q10', 
                     label = TRUE)
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

png(file=sprintf('./3_Integration/Foxp3'),width=12,height=10,units="in",res=300)
FeaturePlot(seurat.combined.sct, 
            reduction = "umap", 
            features = c("Foxp3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()

#######-renamed_cluster-------------------------------------------------------------------
Idents(seurat.combined.sct) <- "renamed_cluster"
#-long noncoding RNA of interest----------------------------------------------------------
markers.to.plot <- c("Snhg18","AI662270","AI427809","Gm6981","Gm10521",
                     "Rian","Snhg5","Snhg6","Bach2os","Gm39822","Cep83os","Gm4951",
                     "Gm10419","Lncbate1","Gm28437","Gm13910","Gm28438")
markers.to.plot <- c("Gm28439","Gm10222","Gm26809","AI838599","Gm10419", "Gm10282",
                     "AI661453","Gm10925","Gm9910","Gm45774","Gm15437",
                     "Gm15590","Gm26917","Gm24270","Gm6768","Gm10252",
                     "5031439G07Rik","4731419I09Rik","6330403L08Rik","Snhg3")
markers.to.plot <- c("Gm17315","Gm12216","Gm15222","A530010L16Rik")
markers.to.plot <- c("Snhg12","Carmn","Gm14461")
## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
    cat("### ",gene, "\n")
    #p <- FeaturePlot(seurat, features = gene, reduction = "umap",
    #                 cols = c("grey", "red")) & theme(legend.position = "bottom")
    png(filename = sprintf('./3_Integration/combined/lncRNA/markers_%s.png', gene), width=6,height=6,units="in",res=300)
    p <- FeaturePlot(seurat.combined.sct, features = gene, reduction = "umap",
                     order = TRUE,
                     min.cutoff = 'q10', 
                     label = TRUE)
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


## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
    cat("### ",gene, "\n")
    cat("**Contrast** \n")
      png(filename = sprintf('./3_Integration/newqc_combined/Violin/markers_%s.png', gene), width=12,height=8,units="in",res=300)
    p1 <- VlnPlot(seurat.combined.sct, 
                  features = gene,
                  pt.size = 0.1,
                  idents = params$contrast) + 
      theme(axis.text.x = element_text(angle = 90, hjust=1))
    print(p1)
    dev.off()
    cat("\n")
    cat("**All groups** \n")
    png(filename = sprintf('./3_Integration/newqc_combined//Violin/fig/markers_%s.png', gene), width=12,height=8,units="in",res=300)
    p2 <- VlnPlot(seurat.combined.sct, features = gene,
                  pt.size = 0.1) +
      theme(axis.text.x = element_text(angle = 90, hjust=1))
    print(p2)
    dev.off()
    cat("\n")
    cat("\n")
  }
  else{
    cat("### ",gene, " \n")
    cat("Not present in dataset \n")
  }
  cat("\n")
}


## ---------------------------------------------------------------------------------------------------------------
cols = gg_color_hue(length(unique(seurat.combined.sct$renamed_cluster)))
#cols = gg_color_hue(length(unique(seurat.combined.sct$seurat_clusters)))

## ---- results = 'asis'------------------------------------------------------------------------------------------
for( gene in markers.to.plot){
  if(gene %in% rownames(seurat.combined.sct)){
    cat("### ",gene, "{.tabset} \n")
    cat("#### Grouped by cluster, split by genotype {.tabset} \n")
    cat("\n")
    cat("##### Side by side \n")
    png(filename = sprintf('./3_Integration/newqc_combined/Violin/fig1/markers_%s.png', gene), width=12,height=8,units="in",res=300)
    p <- VlnPlot(seurat.combined.sct, features = gene,
                 pt.size = 0.1, group.by = "renamed_cluster") +
      theme(axis.text.x = element_text(angle = 90, hjust=1))
    print(p)
    dev.off()
    cat("\n")
    cat("\n")
    cat("##### Same violin \n")
    png(filename = sprintf('./3_Integration/newqc_combined/Violin/fig2/markers_%s.png', gene), width=12,height=8,units="in",res=300)
    p <- VlnPlot(seurat.combined.sct, features = gene,
                 pt.size = 0.1, group.by = "renamed_cluster") +
      theme(axis.text.x = element_text(angle = 90, hjust=1))
    print(p)
    dev.off()
    cat("\n")
    cat("\n")
    
    cat("#### Grouped by genotype, split by cluster\n")
    png(filename = sprintf('./3_Integration/newqc_combined/Violin/fig3/markers_%s.png', gene), width=12,height=8,units="in",res=300)    
    p <- VlnPlot(seurat.combined.sct, features = gene,
                 pt.size = 0.1, split.by = "renamed_cluster", cols = cols) +
      theme(axis.text.x = element_text(angle = 90, hjust=1),
      )
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


