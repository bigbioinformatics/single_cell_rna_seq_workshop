### scRNA-Seq Analysis Workshop Module #1: Loading data into Seurat ###

# We will be using a public 10X Genomics scRNA-Seq dataset in PBMC cells
# You can read more 10X public datasets: https://www.10xgenomics.com/resources/datasets/
# Download the dataset if it is not already present
if (! file.exists("filtered_gene_bc_matrices/hg19/barcodes.tsv")) {
  download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz")
  untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
} 

## Following the Seurat Vignette: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html ##
# Load the Seurat library
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# What kind of object is this?
pbmc.data

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data)

# What kind of object is this?
pbmc

## That's all for today! ##
## Activity: Try finding a new single cell dataset and loading it into R ##
## Datasets can be found here: https://www.10xgenomics.com/resources/datasets/ ##


