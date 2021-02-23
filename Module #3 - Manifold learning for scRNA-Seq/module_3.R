#### Module #3: Manifold Learning for scRNA-Seq #### 

### Part #0: Recap ###

## Load libraries ##
library(Seurat)
library(tidyverse)
library(ggpubr)
library(Matrix)

## Get the data ##
download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
              destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz",
              method = "curl")
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
mat <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(mat)

## Do QC ##

# Percent MT content
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC plots
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Make the cutoff
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Normalize and Scale ##

# Log-transform the counts
pbmc <- NormalizeData(pbmc)

# Find Variable Features
pbmc <- FindVariableFeatures(pbmc)

# Scale the data
pbmc <- ScaleData(pbmc)

## Use a pipe instead ##
pbmcraw <- CreateSeuratObject(mat)

# QC/Norm/scale pipeline
pbmcraw %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# QC/Norm/scale pipeline (save the results)
pbmc <- pbmcraw %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# With plotting
pbmcraw %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "percent.mt")

# With ggplot2
pbmcraw %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  pluck("meta.data") %>%
  ggplot(aes(y = percent.mt, x = nCount_RNA)) +
  geom_point() +
  theme_bw()

### Part #1: Principal Component Analysis ###

# Run PCA
pbmc <- RunPCA(pbmc)

# Plot PCs 1 and 2
DimPlot(pbmc)

# What are the genes that determine PC1 and PC2?
VizDimLoadings(pbmc, dims = 1:2)

# FeaturePlot of PC1 top features
FeaturePlot(pbmc, features = c("CST3", "MALAT1"))

# FeaturePlot of PC2 top features
FeaturePlot(pbmc, features = c("CD79A", "NKG7"))

# Select number of PCs
ElbowPlot(pbmc, ndims = 50)

## Activity: do the PCA on the GBM dataset ##
# Data page: https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_SC3v3_Human_Glioblastoma
# Get the GBM dataset
dataset_url <- "https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_filtered_feature_bc_matrix.h5"

# Download and untar
download.file(dataset_url, method = "curl",
          destfile = "Parent_SC3v3_Human_Glioblastoma_filtered_feature_bc_matrix.h5")

# Make Seurat object
gbm.mat <- Read10X_h5("Parent_SC3v3_Human_Glioblastoma_filtered_feature_bc_matrix.h5")
gbm <- CreateSeuratObject(gbm.mat)

# Run the PreProcessing pipeline
gbm <- gbm %>%
  _____(pattern = "^MT-", col.name = "percent.mt") %>%
  _____(subset = nFeature_RNA > 200 & nFeature_RNA < 8000 &
           percent.mt < 20) %>%
  _____() %>%
  _____() %>%
  _____() 

# Run PCA 

# Plot PCs 1 and 2

# What are the genes that determine PC1 and PC2?

# FeaturePlot of PC1 top features

# FeaturePlot of PC2 top features

# Make elbow plot and identify number of PCs to keep for downstream


### Part #2: Nearest neighbor graph, clustering, and embedding ###

# Find nearest neighbors and construct the graph
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find the clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Get the UMAP embedding
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot the UMAP with clustering
DimPlot(pbmc, reduction = "umap")

# Get the TSNE embedding
pbmc <- RunTSNE(pbmc, dims = 1:10)

# Plot the UMAP with clustering
DimPlot(pbmc, reduction = "tsne")

## Activity: Do the same for the GBM dataset ##

# Find nearest neighbors and construct the graph

# Find the clusters

# Get the UMAP embedding

# Plot the UMAP with clustering

# Get the TSNE embedding

# Plot the UMAP with clustering

### Part #3: Cluster Marker Identification ###

# UMAP for PBMC
DimPlot(pbmc, reduction = "umap")

# Find markers for cluster 3
cluster3.markers <- FindMarkers(pbmc, ident.1 = 3,
                                logfc.threshold = .5,
                                min.pct = 0.25)
cluster3.markers %>%
  top_n(5, -log10(p_val))

# Find markers of 0, 2, 4, and 6
cluster0246.markers <- FindMarkers(pbmc, ident.1 = c(0, 2, 4, 6),
                                   logfc.threshold = .5,
                                   min.pct = 0.25)
cluster0246.markers %>%
  top_n(5, -log10(p_val))

# Find markers that distinguish cluster 0 from 2, 4, and 6
cluster4.dist.markers <- FindMarkers(pbmc, ident.1 = 6,
                                     ident.2 = c(0, 2, 4),
                                     logfc.threshold = .5,
                                     min.pct = 0.25)
cluster4.dist.markers %>%
  top_n(5, -log10(p_val))

# Find all markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.5)
marker_summary <- pbmc.markers %>% 
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC)


# Make a violin plot for strong markers
VlnPlot(pbmc, features = c("CD79A", "CD8A", "CD14", "CCR7",
                           "PPBP", "NKG7"))

# Demonstrate subtypes with featureplots
FeaturePlot(pbmc, features = c("CD79A", "CD8A", "CD14", "CCR7",
                               "PPBP", "NKG7"))

# Assemble a marker list
# Derived from http://biocc.hrbmu.edu.cn/CellMarker/index.jsp
marker_list <- list(
  'T-Cells (CD4)' = c("NOSIP", "CCR7", "IL7R", "CD27"),
  'T-Cells (CD8)' = c("CD8A", "CD8B", "GZMK", "CD3E"),
  'B-Cells' = c("CD79A", "CD79B"),
  "NK Cells" = c("NKG7", "GZMB"),
  "Platelets" = "PPBP"
)

# Which are the CD4 T cells?
FeaturePlot(pbmc, features = marker_list$`T-Cells (CD4)`)

# Which are the CD8 T cells?
FeaturePlot(pbmc, features = marker_list$`T-Cells (CD8)`)

# Which are the B-cells?
FeaturePlot(pbmc, features = marker_list$`B-Cells`)

# Which are the NK-cells?
FeaturePlot(pbmc, features = marker_list$`NK Cells`)

# Which are the platelets?
FeaturePlot(pbmc, features = marker_list$Platelets)

## Do this on the GBM dataset. Which are the tumor cells? ##

# UMAP of GBM
DimPlot(gbm, reduction = "umap")

# Find all markers


# Markers derived from http://biocc.hrbmu.edu.cn/CellMarker/index.jsp
marker_list <- list(
  'Oligodendrocytes' = c("OLIG2", "FA2H", "UGT8", "CNP"),
  'Macrophages' = c("CD68", "CD74", "CD93", "ARG2"),
  'T-cells' = c("CXCR3", "CD8A", "IL7R", "CD3E"),
  'Endothelial Cells' = c("CD34", "ESAM", "PECAM1", "CDH5"),
  'GBM Cells' = c("SOX2", "PARP1", "NES", "EGFR"),
  'Astrocytes' = c("GFAP", "SOX9", "CBS", "ATP13A4")
)

# Using FeaturePlots, answer the following...
# 1. Which are the Oligodendrocytes?

# 2. Which are the Macrophages?

# 3. Which are the T cells?

# 4. Which are the Endothelial cells? 

# 5. Which are the GBM Tumor cells? 

# 6. Which are the Astrocytes? 







