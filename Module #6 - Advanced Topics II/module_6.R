#### Module #6: Advanced Topics II #### 

# Check for missing R packages, install if needed
list_of_pkgs <- c("hdf5r", "tidyverse", "ggpubr", "Seurat", "plotly", "phateR",
                  "cowplot", "patchwork", "VennDiagram", "devtools", "UpSetR",
                  "reticulate", "RColorBrewer", "pheatmap", "Rmagic")
install.packages(list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())])

# Check for missing BioC packages, install if needed
list_of_pkgs <- c("muscat")
bioc2install <- list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())]
if (length(bioc2install) > 0) {
  BiocManager::install(bioc2install)
}

### IMPORTANT: Download and install Python and python packages ###

# Method #1 -- through reticulate
reticulate::install_miniconda()
Rmagic::install.magic()
phateR::install.phate()

# Method #2 -- classic install
# Visit https://www.python.org/downloads/
# Download python and set it on your PATH
# Once you do this, restart rstudio
# This is not strictly necessary, but you will be able to follow along if you do it. 

## Try the test example -- MAKE SURE THIS WORKS BEFORE CONTINUING ##
# MAGIC
library(Rmagic)
data(magic_testdata)
MAGIC_data <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"))
ggplot2::ggplot(MAGIC_data) +
  ggplot2::geom_point(ggplot2::aes(x=VIM, y=CDH1, color=ZEB1))

# PHATE
library(phateR)
data(tree.data)
tree.phate <- phate(tree.data$data)
plot(tree.phate, col = tree.data$branches)


## Load libraries ##
library(Seurat)
library(VennDiagram)
library(tidyverse)
library(ggpubr)
library(phateR)
library(plotly)
library(Rmagic)
library(UpSetR)
library(cowplot)
library(SingleCellExperiment)
library(muscat)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(reticulate)

# Helper function for plotting enrichr results
WtOrRd_pal <- colorRampPalette(c("#FFFFFF", brewer.pal(9, "OrRd")))(100)
# Get enrichr databases
dbs <- listEnrichrDbs()
# Set random seed
set.seed(42) 

## Get the data ##
# From http://adsn.ddnetbio.com/
srt <- read_tsv("scRNA_rawCounts.tsv.gz") %>%
  column_to_rownames(var = "geneName") %>%
  CreateSeuratObject(meta.data = read_tsv("scRNA_metadata.tsv") %>%
                       column_to_rownames("sampleID")) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters(resolution = .8) %>%
  RunUMAP(dims = 1:20, n.neighbors = 60)

# UMAP plot
DimPlot(srt, group.by = "cellType") + DimPlot(srt, group.by = "batchCond") +
  DimPlot(srt, group.by = "sex") + DimPlot(srt, group.by = "patient")


### Part 1: Differential State Analysis ###

# Goal: Identify the DEGs between AD / CTR per each cell type

## 1A. Using Seurat ##

# Assign Idents as the cell type
Idents(srt) <- srt$cellType

# Run the DGE analysis using AD / CTR within Specific Cell Types
ad_vs_ctr_mg <- FindMarkers(srt, ident.1 = "AD", ident.2 = "ct",
                            group.by = "batchCond", subset.ident = "mg")

# Do the same with all cell types and collect in a list
resList <- lapply(levels(Idents(srt)), function(cellType_now) {
  FindMarkers(srt, ident.1 = "AD", ident.2 = "ct", only.pos = TRUE,
              group.by = "batchCond", subset.ident = cellType_now) %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < .05) %>%
    pull(gene)
})
names(resList) <- levels(Idents(srt))
resListFinal <- resList[-which(names(resList) %in% c("unID", "doublet", "endo"))]

# Use venn diagram to compare
venn.diagram(resListFinal, filename = "compared_celltypes_via_seurat.png", 
             margin = .1, fill = c("firebrick", "skyblue", "goldenrod",
                                   "forestgreen", "pink"))
ols <- calculate.overlap(resListFinal)
ols$a31

# FeaturePlots
VlnPlot(srt, group.by = "cellType", features = "LINGO1", split.by = "batchCond")
VlnPlot(srt, group.by = "cellType", features = "XIST", split.by = "batchCond")

# Hmm... XIST is curious. What is the effect of sex?
VlnPlot(srt, group.by = "cellType", features = "XIST", split.by = "sex")
# XIST gene: 
tb <- table(srt$sex, srt$batchCond)
tb / rowSums(tb) * 100

# Use an Upset Plot to compare
upset(fromList(resListFinal))

## 1B. Using Pseudobulk analysis with muscat ##
# Convert to sce
sce <- as.SingleCellExperiment(srt)

# Add columns for muscat
sce <- prepSCE(sce, 
               kid = "cellType", # subpopulation IDs (e.g., cell types)
               gid = "batchCond",  # condition (e.g., ctrl/AD)
               sid = "patient",   # sample IDs (e.g., patient ID)
               drop = TRUE)

# Aggregate data into pseudobulk
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
View(pb@assays@data$astro)

# MDS plot of pseudobulk counts
pbMDS(pb)

# run DS analysis
res <- pbDS(pb, verbose = TRUE)

# Get the results for the MicroGlia
mg_res <- res$table$ct$mg
View(mg_res)

# Get the significant results for all and find concordance between cell-types
resList <- lapply(res$table$ct, function(table_now) {
  table_now %>%
    filter(p_adj.loc < .05 & logFC < 0) %>%
    pull(gene)
})
resListFinal <- resList[-which(names(resList) %in% c("unID", "doublet", "endo"))]

# Use venn diagram to compare
venn.diagram(resListFinal, filename = "compared_celltypes_via_muscat.png", 
             margin = .1, fill = c("firebrick", "skyblue", "goldenrod",
                                   "forestgreen", "pink"))
ols <- calculate.overlap(resListFinal)
ols$a31

# Use an Upset Plot to compare
upset(fromList(resListFinal))

### Part 2: Graph Signal Processing ###

# Graph signal processing in scRNA-Seq is largely based on the work
# of the Krishnaswamy lab at yale: https://www.krishnaswamylab.org/

## 2A. MAGIC ##

# Why MAGIC?
FeatureScatter(srt, 
               group.by = "batchCond",
               feature1 = "XIST", 
               feature2 = "LINGO1")

# Run MAGIC
data_MAGIC <- magic(srt@assays$RNA@data, knn=15, seed = 42)

# Create Seurat with MAGIC data
srtMAGIC <- SetAssayData(srt, slot = "data", new.data = as.matrix(data_MAGIC$result))

# Recalculate with MAGIC
srtMAGIC <- srtMAGIC %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters(resolution = .8) %>%
  RunUMAP(dims = 1:20, n.neighbors = 60)

# DimPlot MAGIC vs Normal
DimPlot(srtMAGIC, group.by = "batchCond") + DimPlot(srt, group.by = "batchCond")

# FeaturePlot of gene expression compared
fs1 <- FeatureScatter(srtMAGIC, 
                      group.by = "batchCond",
                      feature1 = "XIST", 
                      feature2 = "LINGO1")
fs2 <- FeatureScatter(srt, 
                      group.by = "batchCond",
                      feature1 = "XIST", 
                      feature2 = "LINGO1")
fs1 + fs2

## 2B. PHATE ##

# Run PHATE
srt_phate <- phate(t(srt@assays$RNA@data))

# Get Embedding
emd <- srt_phate$embedding

# Add into SRT
srt$PHATE1 <- emd[,c(1)]
srt$PHATE2 <- emd[,c(2)]

# Plot the PHATE embedding
DimPlot(srt) + FeatureScatter(srt, feature1 = "PHATE1", 
                              plot.cor = FALSE, 
                              feature2 = "PHATE2") +
  DimPlot(srt, group.by = "batchCond") + FeatureScatter(srt, feature1 = "PHATE1", 
                                                        plot.cor = FALSE, group.by = "batchCond",
                                                        feature2 = "PHATE2")








