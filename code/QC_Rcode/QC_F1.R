rm(list = ls())
setwd('/Users/hueee/Lab/Feinberg\ Lab/Single_cell_adipose_tissue/')
#knitr::purl('QC1.Rmd')
#setwd('/home/jic621/maomao/single_cell_August')
#load library
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

#create output folder
output_dir = file.path('./1_QC/F1')
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Path already exists!")
}

## setup filters
params <-
  list(matricesDir = "./data/F1/", 
       sample_name = "Chow", outputDir = ".", minGenes = 400L, 
       minUMIs = 600L, maxMitoPer = 20L, minNovelty = 0.8)
params <-
  list(matricesDir = "./data/F1/", 
       sample_name = "Chow", outputDir = ".", minGenes = 500L, 
       minUMIs = 800L, maxMitoPer = 18L, minNovelty = 0.8)
#sample_name = "Cre_01"
#minGenes = 400
#minUMIs = 600
#maxMitoPer = 20
#minNovelty = 0.8

# Load the dataset
# Initialize the Seurat object with the raw (non-normalized data).
seurat.data <- Read10X(data.dir = params$matricesDir)
seurat <- CreateSeuratObject(counts = seurat.data, project = "Visceral_Chow", min.cells = 3, min.features = 200)
seurat <- CreateSeuratObject(counts = seurat.data[!duplicated(row.names(seurat.data)),], 
                             project = "Visceral_Chow", 
                             min.cells = 3, 
                             min.features = 200)
  
  #hca_bm <- CreateSeuratObject(counts = hca_bm.counts[!which(duplicated(row.names(hca_bm.counts))),], project = "hca_bm", min.cells = 5, min.features =200)
seurat[["total_reads"]] <- Matrix::colSums(GetAssayData(seurat,slot = "counts"))
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
seurat[["novelty"]] <- log10(seurat$nFeature_RNA)/log10(seurat$nCount_RNA)
saveRDS(seurat,"data/seurat_Chow.rds")
#The initial dataset after cellranger filtering and seurat initial filtering (a gene must have expression in at least 3 cells and a cell has to have at least 200 genes) has `r ncol(seurat)` cells.

initial_cells <- ncol(seurat)

# QC Plots {.tabset}

## nGenes

#With the distribution of the number of genes we observe that most cells have a very low number of genes. 
#We set a low threshold (red line)of `r params$minGenes` number of Genes.
metadata = seurat@meta.data
#------------------------------------
png(filename = './1_QC/F1/fig_1.png')
ggplot(data = seurat[[]], 
       aes(x = nFeature_RNA)) + geom_histogram(binwidth = 50) + geom_vline(xintercept = params$minGenes, 
       col = "red", linetype = "dashed") + xlab("nGenes")
dev.off()
#------------------------------------
png(filename = './1_QC/F1/fig_2.png')
VlnPlot(seurat, features = "nFeature_RNA") + 
  geom_hline(yintercept = params$minGenes, 
             col = "red", linetype = "dashed")
dev.off()
## ---------------------------------------------------------------------------------------------------------------
nGenes_df <- seurat[[]] %>% as.data.frame() %>% tibble::rownames_to_column(var  = "barcode") %>% 
  dplyr::select(c(barcode, nFeature_RNA)) %>% dplyr::arrange(-nFeature_RNA)
nGenes_df$ranked_bc <- as.integer(rownames(nGenes_df))

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_3.png')
p <- ggplot() +
  geom_line(data = nGenes_df, aes(x = ranked_bc, y = nFeature_RNA), col = "red")

p  +
  scale_y_continuous(trans='log10') + geom_hline(yintercept = params$minGenes, col = "red", linetype = "dashed") + ylab("nGenes - log10") + xlab("barcode")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_4.png')
ggplot(data = seurat[[]], aes(x = nCount_RNA)) + geom_histogram(binwidth = 50)  + geom_vline(xintercept = params$minUMIs, col = "red", linetype = "dashed") + xlab("nUMIs") + xlim(0,5000)
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_5.png')
VlnPlot(seurat, features = "nCount_RNA") + geom_hline(yintercept = params$minUMIs, col = "red", linetype = "dashed") 
dev.off()

## ---------------------------------------------------------------------------------------------------------------
nUMIs_df <- seurat[[]] %>% as.data.frame() %>% tibble::rownames_to_column(var  = "barcode") %>% 
  dplyr::select(c(barcode, nCount_RNA)) %>% dplyr::arrange(-nCount_RNA)
nUMIs_df$ranked_bc <- as.integer(rownames(nUMIs_df))

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_6.png')
p <- ggplot() +
  geom_line(data = nUMIs_df, aes(x = ranked_bc, y = nCount_RNA), col = "red")

p  +
  scale_y_continuous(trans='log10') + geom_hline(yintercept = params$minUMIs, col = "red", linetype = "dashed") + ylab("nUMIs - log10") + xlab("barcode")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_7.png')
ggplot(data = seurat[[]], aes(x = percent.mt)) + geom_histogram(binwidth = 1) + xlim(0,100) + geom_vline(xintercept = params$maxMitoPer, col = "red", linetype = "dashed")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_8.png')
VlnPlot(seurat, features = "percent.mt") + geom_hline(yintercept = params$maxMitoPer, col = "red", linetype = "dashed")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_9.png')
seurat[[]] %>% 
  ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  geom_boxplot() + 
  ggtitle("NGenes")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_10.png')
seurat[[]] %>% 
  ggplot(aes(x = log10(nCount_RNA), y = log10(nFeature_RNA),color = percent.mt)) + scale_color_gradient(limits=c(0, 100)) +
  geom_point() + 
  stat_smooth(method = lm) + 
  facet_wrap(~orig.ident) + xlab("log10 nUMis") + ylab("log10 nGenes")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_11.png')
seurat[[]] %>%
  ggplot(aes(x=novelty))+
  geom_histogram(bins = 50) 
dev.off()

## ---------------------------------------------------------------------------------------------------------------
seurat <- subset(seurat, subset = nFeature_RNA > params$minGenes & percent.mt < params$maxMitoPer & nCount_RNA >params$minUMIs & novelty > params$minNovelty) 


## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_12.png')
ggplot(data = seurat[[]], aes(x = nFeature_RNA)) + geom_histogram(binwidth = 50) + xlab("nGenes")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_13.png')
ggplot(data = seurat[[]], aes(x = nCount_RNA)) + geom_histogram(binwidth = 50) + xlab("nUMIs")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_14.png')
ggplot(data = seurat[[]], aes(x = percent.mt)) + geom_histogram(binwidth = 1) + xlim(0,100) 
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_15.png')
seurat[[]] %>% 
  ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  geom_boxplot() + 
  ggtitle("NGenes")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_16.png')
seurat[[]] %>% 
  ggplot(aes(x = log10(nCount_RNA), y = log10(nFeature_RNA),color = percent.mt)) +  scale_color_gradient(limits=c(0, 100))+
  geom_point() + 
  stat_smooth(method = lm) + 
  facet_wrap(~orig.ident) + xlab("log10 nUMis") + ylab("log10 nGenes")
dev.off()

## ---------------------------------------------------------------------------------------------------------------
png(filename = './1_QC/F1/fig_17.png')
seurat[[]] %>%
  ggplot(aes(x=novelty))+
  geom_histogram(bins = 50) 
dev.off()

## ---------------------------------------------------------------------------------------------------------------
saveRDS(seurat,"data/seurat_Chow_filtered.rds")
# Conclusions QC
#From the `r initial_cells` initial cells, `r ncol(seurat)` passed the filtering 
#(`r round(ncol(seurat)*100/initial_cells,2)`%).
getwd()
seurat_Chow_filtered <- readRDS("data/seurat_Chow_filtered.rds")
summary(seurat_Chow_filtered$nCount_RNA)
summary(seurat_Chow_filtered$nFeature_RNA)
summary(seurat_Chow_filtered$percent.mt)
summary(seurat_Chow_filtered$novelty)
