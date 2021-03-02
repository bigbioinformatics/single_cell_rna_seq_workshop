### Wes Wilon's T-cell Pipeline, Minor Edits for Henry's scRNA-seq Class 20210202
#01101000 01100001 01100011 01101011  01110100 01101000 01100101  01110000 01101100 01100001 01101110 01100101 01110100
#__    __    ___  _____ __      ______         __    ___  _      _      _____
#|  T__T  T  /  _]/ ___/T  |    |      T       /  ]  /  _]| T    | T    / ___/
#|  |  |  | /  [_(   \_ l_ |    |      |      /  /  /  [_ | |    | |   (   \_ 
#|  |  |  |Y    _]\__  T  \l    l_j  l_j     /  /  Y    _]| l___ | l___ \__  T
#l  `  '  !|   [_ /  \ |          |  |      /   \_ |   [_ |     T|     T/  \ |
# \      / |     T\    |          |  |      \     ||     T|     ||     |\    |
#  \_/\_/  l_____j \___j          l__j       \____jl_____jl_____jl_____j \___j
#01101000 01100001 01100011 01101011  01110100 01101000 01100101  01110000 01101100 01100001 01101110 01100101 01110100


# Load Libraries

# scRNA-seq
library("DropletUtils")
library("SingleCellExperiment")
library("scater")
library("scran")

# Plotting
library("cowplot")
library("UpSetR")
library("grid")

# Tidyverse
library("tidyverse")

#custom Functions
source("load.R")
source("annotate.R")
source("output.R")
source("plotting.R")

#Set Cores
bpparam <- BiocParallel::MulticoreParam(workers = 24)


# This Reads in the RAW counts
raw <- read10xCounts("../2020_TcellRepeats/Repeat4/outs/raw_feature_bc_matrix/")


# Here is what cellranger thinks is a cell
# Read Filtered barcodes
filt_barcodes <- readr::read_lines(file.path("../2020_TcellRepeats/Repeat4/outs/filtered_feature_bc_matrix", "barcodes.tsv.gz"))

colData(raw)$CellRangerFilt <- colData(raw)$Barcode %in% filt_barcodes

summary(colData(raw)$CellRangerFilt)


default1 <- defaultDrops(counts(raw))
colData(raw)$DefaultFilt <- c(default1)  #This is a list becasue I normally process multiple samples at once


# Remove GEM's with zero genes & remove genes with 0 reads
raw <- raw[rowSums(counts(raw)) > 0, colSums(counts(raw)) > 0]


#@ Creatting Rank Plot , Cell selection QC
empty_thresh <- 100

bc_ranks <- barcodeRanks(counts(raw), lower = empty_thresh)

colData(raw)$BarcodeRank   <- bc_ranks$rank
colData(raw)$BarcodeTotal  <- bc_ranks$total
colData(raw)$BarcodeFitted <- bc_ranks$fitted

bc_data <- colData(raw) %>%
  as.data.frame() %>%
  select(Barcode, Kept = CellRangerFilt, Rank = BarcodeRank,
         Total = BarcodeTotal, Fitted = BarcodeFitted) %>%
  arrange(Rank)


# Generating Plot, this is just a ggplot, and this isn't a ggplot tutorial so will leave it at that.
# I stole the colours off of Luke Zappia's 2018 plots, although where the data is located and stored
# I have put in manually to work with the newest versions of Bioconductor packages, and this works as of 20210302
ggplot(bc_data, aes(x = Rank, y = Total)) +
  geom_point(shape = 1, aes(colour = Kept)) +
  geom_line(aes(y = Fitted), colour = "red") +
  geom_hline(yintercept = bc_ranks@metadata$knee, colour = "dodgerblue", linetype = "dashed") +
  annotate("text", x = 0, y = bc_ranks@metadata$knee, label = "Knee", colour = "dodgerblue", hjust = 0, vjust = -1) +
  geom_hline(yintercept = bc_ranks@metadata$inflection, colour = "forestgreen", linetype = "dashed") +
  annotate("text", x = 0, y = bc_ranks@metadata$inflection, label = "Inflection", colour = "forestgreen", hjust = 0, vjust = -1) +
  geom_hline(yintercept = empty_thresh, colour = "darkorchid", linetype = "dashed") +
  annotate("text", x = 0, y = empty_thresh, label = "Empty threshold", colour = "darkorchid", hjust = 0, vjust = -1) +
  scale_x_log10(labels = scales::number) +
  scale_y_log10(labels = scales::number) +
  scale_colour_manual(values = c("black", "violet")) +
  ylab("Total counts") +
  theme_minimal()


### Remove Empties /w EmptyDrops (this is gold standard approach and found in EmptyDrops vignette)
set.seed(1337)
emp_iters <- 30000
emp_fdr <- 0.01

emp_drops <- emptyDrops(counts(raw), lower = empty_thresh, niters = emp_iters, test.ambient = TRUE, BPPARAM = bpparam)

is_cell <- emp_drops$FDR <= emp_fdr
is_cell[is.na(is_cell)] <- FALSE

colData(raw)$EmpDropsLogProb <- emp_drops$LogProb
colData(raw)$EmpDropsPValue  <- emp_drops$PValue
colData(raw)$EmpDropsLimited <- emp_drops$Limited
colData(raw)$EmpDropsFDR     <- emp_drops$FDR
colData(raw)$EmpDropsFilt    <- is_cell

table(Limited = emp_drops$Limited, Significant = is_cell)


#### Comparison Visualization Idea stolen from Luke Zappia 2018, Updated to work with modern SingleCellExperiment on 10X data by Wes Wilson 2020

plot_data <- list(
  "Cell Ranger"    = colData(raw)$Barcode[colData(raw)$DefaultFilt],
  "Cell Ranger v5" = colData(raw)$Barcode[colData(raw)$CellRangerFilt],
  "Drop Method"     = colData(raw)$Barcode[colData(raw)$EmpDropsFilt]
)

upset(fromList(plot_data), order.by = "freq", matrix.color="darkslategrey", sets.bar.color= c("coral3","cornflowerblue","coral"),
      sets.x.label = "Number of cells", text.scale = c(2, 1.2, 2, 1.2, 2, 3))


#Use "cells" from BOTH methods to start analysis

selected <- raw[, colData(raw)$CellRangerFilt | colData(raw)$EmpDropsFilt]
selected <- selected[Matrix::rowSums(counts(selected)) > 0, ]

colData(selected)$SelMethod <- "Both"
colData(selected)$SelMethod[!colData(selected)$CellRangerFilt] <- "emptyDrops"
colData(selected)$SelMethod[!colData(selected)$EmpDropsFilt] <- "CellRanger"


### Data Clean Up

# scRNA-seq
library("scater")
library("scran")

# RNA-seq
library("edgeR")

# Plotting
library("cowplot")
library("gridExtra")

# ensembl changed list function past this version, students use different method
# I recomend library(AnnotationHub) 

#Change RowNames to GeneSybols
rownames(selected) <- rowData(selected)$Symbol

selected <- annotateSCE(selected,
                       org        = "human",
                       add_anno   = FALSE,
                       host       = "jan2019.archive.ensembl.org",
                       calc_qc    = TRUE,
                       calc_cpm   = TRUE,
                       cell_cycle = TRUE,
                       BPPARAM    = bpparam,
                       verbose = TRUE)


glimpse(as.data.frame(colData(selected)))

set.seed(1)

## Plot Highest Expressed Genes, this is a good QC method
plotHighestExprs(selected)


## Normalize library sizes
sizeFactors(selected) <- librarySizeFactors(selected)

#Log2
selected <- logNormCounts(selected)

selected <- runPCA(selected)
selected <- runTSNE(selected)
selected <- runUMAP(selected)


## Detect Doublets
dbl.dens <- doubletCells(selected)
summary(dbl.dens)


#Graph Cell Cycle
library(ggplot2)

cell_data <- as.data.frame(colData(selected))

ggplot(cell_data, aes(x = G1Score, y = G2MScore, colour = CellCycle, shape = Sample)) +
  geom_point() +
  xlab("G1 score") +
  ylab("G2/M score") +
  theme_minimal()



#SingleR by Aaron Lun, Jared M. Andrews1, Friederike Dündar2 and Daniel Bunis
#Checkout scMatch by Al Forrest's Group too...

library(SingleR)
set.seed(100)
#Try HumanPrimaryCellAtlasData() as a database, its all about what cell types you expect
ref1 <- BlueprintEncodeData()
pred <- SingleR(test = selected, ref = ref1, labels = ref1$label.main, BPPARAM= bpparam)
table(pred$labels)
plotScoreHeatmap(pred)

# In real life I use library("irlba") for my PCA method , and continue with normal analysis stuff


etc etc etc


### YOU CAN DO EVERYTHING IN Bioconductor / cran for you're entire analysis BUT you can jump into Seurat anytime you want
## dont use as.Seurat , ITS A TRAP!!!!
library(Seurat)
library(patchwork)
tcounts <- as.matrix(counts(selected))
colnames(tcounts) <- selected$Barcode
tcells.so <- CreateSeuratObject(tcounts, min.cells = 3, min.features = 200, names.field = selected$Barcode) 


### OLD WAY Still works
mt.genes <- rownames(tcells.so)[grep("^MT-",rownames(tcells.so))]
C<-GetAssayData(object = tcells.so, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
tcells.so <- AddMetaData(tcells.so, percent.mito, col.name = "percent.mito")

rb.genes <- rownames(tcells.so)[grep("^RP[SL]",rownames(tcells.so))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
tcells.so <- AddMetaData(tcells.so, percent.ribo, col.name = "percent.ribo")

### NEW Seurat Way , USE THIS KIDS!

tcells.so[["percent.mt"]] <- PercentageFeatureSet(tcells.so, pattern = "^MT-")
tcells.so[["percent.rb"]] <- PercentageFeatureSet(tcells.so, pattern = "^RP[SL]")


VlnPlot(tcells.so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)


tcells.so <- NormalizeData(tcells.so, normalization.method = "LogNormalize", scale.factor = 10000)

tcells.so <- FindVariableFeatures(tcells.so, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tcells.so), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(tcells.so)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


tcells.so <- subset(tcells.so, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

tcells.so <- SCTransform(tcells.so, verbose = FALSE)

tcells.so <- RunPCA(tcells.so)
tcells.so <- RunUMAP(tcells.so, dims = 1:30)

tcells.so <- FindNeighbors(tcells.so, dims = 1:30)
tcells.so <- FindClusters(tcells.so)

library("Nebulosa")

p4 <- plot_density(tcells.so, c("CD8A", "CCR7"), joint = TRUE)
p4 + plot_layout(ncol = 1)

