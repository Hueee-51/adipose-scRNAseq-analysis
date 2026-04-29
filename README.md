# Adipose snRNA-seq Analysis

This repository contains the public analysis code skeleton for single-nucleus RNA-seq analysis of adipose tissue under diet-induced metabolic stress.

The scripts are organized to document the major computational steps used in the study, including quality control, sample integration, major cell type annotation, cell type-specific differential expression, subcluster analysis, and CellChat-based cell-cell communication analysis.

Large sequencing files, processed Seurat objects, and intermediate R objects are not included in this repository.

## Repository structure

```text
adipose-snRNAseq-analysis_public/
├── README.md
├── .gitignore
├── scripts/
│   ├── 01_quality_control.R
│   ├── 02_integration.R
│   ├── 03_clustering_annotation.R
│   ├── 04_major_celltype_DE_analysis.R
│   ├── 05_aspc_subcluster_analysis.R
│   ├── 06_endothelial_subcluster_analysis.R
│   ├── 07_macrophage_subcluster_analysis.R
│   ├── 08_cellchat_major_celltype_analysis.R
│   └── 09_cellchat_subcluster_analysis.R
├── config/
│   ├── sample_metadata_template.csv
│   ├── analysis_parameters.yaml
│   └── celltype_annotation_template.csv
├── data/
│   └── README.md
└── results/
    └── README.md