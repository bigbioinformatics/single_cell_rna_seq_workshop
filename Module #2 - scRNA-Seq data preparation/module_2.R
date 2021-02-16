#### scRNA-Seq Analysis Workshop Module #2: Data Prep & QC ####

# Check for missing packages, install if needed
list_of_pkgs <- c("hdf5r", "tidyverse", "ggpubr")
install.packages(list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())])

# Load libraries
library(Seurat)
library(tidyverse)
library(ggpubr)
library(Matrix)

### Part #1 - Loading Data from Different Sources ###

## 10X Genomics Website (Filtered Feature Barcode folder) ##

# We will be using a public 10X Genomics scRNA-Seq dataset in PBMC cells
# You can read more 10X public datasets: https://www.10xgenomics.com/resources/datasets/
# Download the dataset if it is not already present
if (! file.exists("filtered_feature_bc_matrix/barcodes.tsv.gz")) {
  pbmc_url <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz"
  download.file(pbmc_url,
                destfile = "pbmc_1k_v3_filtered_gene_bc_matrices.tar.gz")
  untar("pbmc_1k_v3_filtered_gene_bc_matrices.tar.gz")
} 
mat <- Read10X(data.dir = "filtered_feature_bc_matrix/")
srt <- CreateSeuratObject(mat)

## 10X Genomics Website (h5 file) ##

# We will be using a public 10X Genomics scRNA-Seq dataset in PBMC cells
# You can read more 10X public datasets: https://www.10xgenomics.com/resources/datasets/
# Download the dataset if it is not already present
if (! file.exists("pbmc_1k_v3_filtered_gene_bc_matrices.h5")) {
  pbmc_h5_url <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5"
  download.file(pbmc_h5_url, method = "curl",
                destfile = "pbmc_1k_v3_filtered_gene_bc_matrices.h5")
} 
mat <- Read10X_h5(filename = "pbmc_1k_v3_filtered_gene_bc_matrices.h5")
srt <- CreateSeuratObject(mat)

## From PangloaDB ##
# Download https://panglaodb.se/data_dl.php?sra=SRA553822&srs=SRS2119548&filetype=R&datatype=readcounts
# Add it to your project folder and then...
load("SRA553822_SRS2119548.sparse.RData")
srt <- CreateSeuratObject(sm)

## GEO ##

# GEO is the main source of all public HTS data
# An example scRNA-Seq study is here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146221

# Step #1: Download the dataset  
if (! file.exists("EWS_Cells_10X/GSM4368462_CHLA9_barcodes.tsv.gz")) {
  archive_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE146221&format=file"
  download.file(archive_url, destfile = "EWS_Cells.tar", method = "curl")  # Method = curl to avoid corruption
  untar("EWS_Cells.tar", exdir = "EWS_Cells_10X") 
}

# Step #2: Get the data into the matrix format in R
Read10X(data.dir = "EWS_Cells_10X/")
TC71_mat <- readMM(file = "EWS_Cells_10X/GSM4368464_TC71_matrix.mtx.gz")
genes <- read_tsv("EWS_Cells_10X/GSM4368464_TC71_features.tsv.gz",
                  col_names = FALSE)
rownames(TC71_mat) <- genes$X2
cellnames <- read_tsv("EWS_Cells_10X/GSM4368464_TC71_barcodes.tsv.gz",
                      col_names = FALSE)
colnames(TC71_mat) <- cellnames$X1
TC71_mat[1:10, 1:4]

# Step #3: Create the Seurat Object
srtTC71 <- CreateSeuratObject(TC71_mat, project = "TC71")

# Step #4: Repeat Steps 2-3 for CHLA9 
CHLA9_mat <- readMM(file = "EWS_Cells_10X/GSM4368462_CHLA9_matrix.mtx.gz")
genes <- read_tsv("EWS_Cells_10X/GSM4368462_CHLA9_features.tsv.gz",
                  col_names = FALSE)
rownames(CHLA9_mat) <- genes$X2
cellnames <- read_tsv("EWS_Cells_10X/GSM4368462_CHLA9_barcodes.tsv.gz",
                      col_names = FALSE)
colnames(CHLA9_mat) <- cellnames$X1
srtCHLA9 <- CreateSeuratObject(CHLA9_mat, project = "CHLA9")

# Step #5: Combine srtCHLA9 and srtTC71
CHLA9_TC71_srt <- merge(x = srtCHLA9, y = srtTC71)

# Activity: Add in CHLA10 to the Seurat object
CHLA10_mat <- readMM(file = "EWS_Cells_10X/GSM4368463_CHLA10_matrix.mtx.gz")
rownames(CHLA10_mat) <- read_tsv("EWS_Cells_10X/GSM4368463_CHLA10_features.tsv.gz",
                                 col_names = FALSE) %>%
  pull(X2)
cellnames <- read_tsv("EWS_Cells_10X/GSM4368463_CHLA10_barcodes.tsv.gz",
                      col_names = FALSE)
colnames(CHLA10_mat) <- cellnames$X1
srtCHLA10 <- CreateSeuratObject(CHLA10_mat, project = "CHLA10")
# Combine srtCHLA10 with srtCHLA9 and srtTC71
full_srt <- merge(x = CHLA9_TC71_srt, y = srtCHLA10)

# Bonus: Is there a more efficient way to do this?
cell_types <- c("CHLA9", "CHLA10", "TC71")
srtList <- lapply(cell_types, function(cell_now) {
  mtx <- list.files(path = "EWS_Cells_10X/", full.names = TRUE, 
                    pattern = paste0(cell_now, "_matrix"))
  features <- list.files(path = "EWS_Cells_10X/", full.names = TRUE,
                         pattern = paste0(cell_now, "_features"))
  barcodes <- list.files(path = "EWS_Cells_10X/", full.names = TRUE,
                         pattern = paste0(cell_now, "_barcodes"))
  
  sm <- readMM(mtx)
  rownames(sm) <- read_tsv(features, col_names = FALSE) %>% pull(X2)
  colnames(sm) <- read_tsv(barcodes, col_names = FALSE) %>% pull(X1)
  
  srt <- CreateSeuratObject(sm, project = cell_now)
  return(srt)
})
full_srt <- merge(x = srtList[[1]], y = srtList[-1])

## Conversion from other formats to Seurat ##

# An excellent vignette (https://satijalab.org/seurat/articles/conversion_vignette.html)
# is available and gives all the likely conversions. Let me know if you want me to 
# go over a specific one. 


### Part #2 - QC with scRNA-Seq datasets ###

rm(list=ls()); gc()  # Cleanup the environment to free up memory

## Step #1: Get the PBMC data into R ##

# We will be using a public 10X Genomics scRNA-Seq dataset in PBMC cells
# You can read more 10X public datasets: https://www.10xgenomics.com/resources/datasets/
# Download the dataset if it is not already present
if (! file.exists("filtered_gene_bc_matrices/hg19/barcodes.tsv")) {
  download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz")
  untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
} 

## Following the Seurat Vignette: 
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html 
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data)

## Step #2: Examine QC metrics ##
# Compare with online summary: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_web_summary.html

# Get the quality info from meta.data
View(pbmc@meta.data)

# Add in the Mitochondrial PCT% information
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# nCount_RNA is the number of UMI counts in a cell
hist(pbmc$nCount_RNA)

# nFeature_RNA is the number of different genes that had any reads
hist(pbmc$nFeature_RNA)

# percent.mt is the percent mitochondrial reads
hist(pbmc$percent.mt)

# Make a violin plot of the QC columns
plt <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3) 
ggsave(filename = "QC.png", plot = plt, width = 7, height = 3.5)

# Make scatter plots to compare
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

# Subset to remove outlier
dim(pbmc)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &
                 nFeature_RNA < 2500 & percent.mt < 5)
dim(pbmc)


## Activity: Prepare and perform the QC step on the dataset below
# Dataset info page: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/neuron_1k_v3
dataset_url <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_filtered_feature_bc_matrix.tar.gz"

# Download and untar
_________(______, method = "curl",
          destfile = ______)
untar(________)

# Make Seurat object
neurons.mat <- Read10X(________)
neurons <- _______(neurons.mat)

# Perform QC check (Hint: What kind of animal is this...?)
neurons$percent.mt <- _________(neurons, pattern = "^mt-")
VlnPlot(_______________) 
plot1 <- FeatureScatter(____________)
plot1
plot2 <- FeatureScatter(____________)
plot2

# Subset and remove outliers -- choose your cutoffs based on the plots
dim(neurons)
neurons <- subset(neurons, subset = ___________)
dim(neurons)

# How many cells did you filter out?

## Step #3: Normalize and Scale the data ##

# Log-transform the counts
pbmc <- NormalizeData(pbmc)

# Find Variable Features
pbmc <- FindVariableFeatures(pbmc)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

# Scale the data
pbmc <- ScaleData(pbmc)

# Compare the raw, log, and scaled count distributions
data.frame(
  counts = c(
    pbmc@assays$RNA@counts[,'AAACATACAACCAC-1'],
    pbmc@assays$RNA@counts[,'AAACATTGAGCTAC-1'],
    pbmc@assays$RNA@data[,'AAACATACAACCAC-1'],
    pbmc@assays$RNA@data[,'AAACATTGAGCTAC-1'],
    pbmc@assays$RNA@scale.data[,'AAACATACAACCAC-1'],
    pbmc@assays$RNA@scale.data[,'AAACATTGAGCTAC-1']
  ),
  cell = c(rep(rep(c('AAACATACAACCAC-1', 'AAACATTGAGCTAC-1'),
             each = length(rownames(pbmc@assays$RNA@counts))), 2),
           rep(c('AAACATACAACCAC-1', 'AAACATTGAGCTAC-1'),
               each = length(rownames(pbmc@assays$RNA@scale.data)))),
  type = factor(c(rep(c("Raw", "Log(x+1)"),
             each = 2*length(rownames(pbmc@assays$RNA@counts))),
             rep("Scaled", 2*length(rownames(pbmc@assays$RNA@scale.data)))),
             levels = c("Raw", "Log(x+1)", "Scaled"))
) %>% 
  ggplot(aes(x = cell, y = counts)) +
  geom_boxplot() +
  theme_bw() +
  rotate_x_text(45) +
  facet_wrap(facets = "type", scales = "free_y")

## Activity: Normalize and scale your neurons dataset as in the above example
# No hints this time -- do your best. 

## Step #4: Perform PCA

# Run PCA
pbmc <- RunPCA(pbmc)

# Plot the PCA
DimPlot(pbmc, reduction = "pca")

# Plot the loadings
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# Choose the number of principle components to keep
ElbowPlot(pbmc)





