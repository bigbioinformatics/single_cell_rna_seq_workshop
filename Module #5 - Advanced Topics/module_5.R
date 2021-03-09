#### Module #5: Advanced Topics #### 

### Part #0: Recap ###

# Check for missing R packages, install if needed
list_of_pkgs <- c("hdf5r", "tidyverse", "ggpubr", "Seurat", "enrichR",
                  "cowplot", "patchwork", "VennDiagram", "devtools",
                  "RColorBrewer", "pheatmap", "msigdbr")
install.packages(list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())])

# Check for missing BioC packages, install if needed
list_of_pkgs <- c("limma", "org.Hs.eg.db", "SingleR", "batchelor")
bioc2install <- list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())]
if (length(bioc2install) > 0) {
  BiocManager::install(bioc2install)
}

# Install SeuratWrappers (optional)
devtools::install_github('satijalab/seurat-wrappers')

# Install monocle3 (optional)
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")

## Load libraries ##
library(Seurat)
library(tidyverse)
library(ggpubr)
library(limma)
library(SingleR)
library(enrichR)
library(cowplot)
library(patchwork)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
library(monocle3)
library(msigdbr)
# Helper function for plotting enrichr results
plot_eres <- function(eres_name, eres_list, n = 10) {
  eres_list[[eres_name]] %>%
    top_n(n = n, wt = -log10(Adjusted.P.value)) %>%
    arrange(-log10(Adjusted.P.value)) %>%
    mutate(Term = factor(Term, levels = Term)) %>%
    ggplot(mapping = aes(x = Term, y = -log10(Adjusted.P.value), fill = Combined.Score)) +
    geom_bar(stat = "identity") +
    ggpubr::rotate() +
    theme_bw(base_size = 16) +
    rremove("ylab") +
    labs(title = eres_name)
}
# Define color palettes
WtOrRd_pal <- colorRampPalette(c("#FFFFFF", brewer.pal(9, "OrRd")))(100)
# Get enrichr databases
dbs <- listEnrichrDbs()
# Set random seed
set.seed(42) 

## Get the data ##
download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
              destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz",
              method = "curl")
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")

# Pipeline to run all previous steps
pbmc <- CreateSeuratObject(Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")) %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>% 
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:10) %>%
  RenameIdents(c(
    "0" = "CD4 + T Cells",
    "1" = "CD14+ Monocytes",
    "2" = "CD4 + T Cells",
    "3" = "B Cells",
    "4" = "CD8+ T Cells",
    "5" = "CD16+ Monocytes",
    "6" = "NK Cells",
    "7" = "Dendritic Cells",
    "8" = "Platelets"
  )) 
pbmc$cell_type <- Idents(pbmc)
Idents(pbmc) <- pbmc$seurat_clusters

# Dim plot of UMAP
DimPlot(pbmc, group.by = "seurat_clusters") + DimPlot(pbmc, group.by = "cell_type")

## Find cluster markers ##
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.5)
sig_markers <- pbmc.markers %>% 
  filter(p_val_adj < .05) %>%
  write_csv(file = "pbmc_sig_markers.csv")

# Cluster 0 enrichr: https://maayanlab.cloud/Enrichr/enrich?dataset=b33d5997396c6721c3bf4bf0e67b8fbe
to_check <- c("Human_Gene_Atlas", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
              "KEGG_2019_Human", "MSigDB_Hallmark_2020")
# Pull the Cluster 0 genes as a vector
cluster0_genes <- sig_markers %>%
  filter(cluster == "0") %>%
  pull(gene)
# Run through enrichr
cluster0_eresList <- enrichr(cluster0_genes, databases = to_check)
# Get the top hits from each and plot them
plotList <- lapply(names(cluster0_eresList), plot_eres, eres_list = cluster0_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 0 Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster0_pbmc_enrichr_plots.png", width = 20, height = 10)

## Find markers of cluster 0 vs cluster 2 ##
cluster0vs2.markers <- FindMarkers(pbmc,
                                   ident.1 = 0, 
                                   ident.2 = 2,
                                   logfc.threshold = .25,
                                   min.pct = 0.25)
# Pull the Cluster 0 vs 2 genes as a vector
cluster0_genes <- cluster0vs2.markers %>%
  filter(avg_log2FC > 0 & p_val_adj < .05) %>%
  rownames_to_column("gene") %>%
  pull(gene)
# Run through enrichr
cluster0_eresList <- enrichr(cluster0_genes, databases = to_check)
# Get the top hits from each and plot them
plotList <- lapply(names(cluster0_eresList), plot_eres, eres_list = cluster0_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 0 vs 2 Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster0_vs2_pbmc_enrichr_plots.png", width = 20, height = 10)


### Part 0.5: Brain dataset marker analysis ###
brain.raw <- read.csv("../Module #4 - Biological Interpretation/Brain_Non-Myeloid-counts.small.csv.gz", 
                      row.names = 1)
brain <- Seurat::CreateSeuratObject(brain.raw) %>%
  subset(subset = nFeature_RNA > 600 & nFeature_RNA < 5000 &
           nCount_RNA < 2.5E6) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution= .5) %>%
  RunUMAP(dims = 1:20) %>%
  RenameIdents(c(
    "0" = "Oligodendrocytes",
    "1" = "Endothelial Cells",
    "2" = "Astrocytes",
    "3" = "Unknown",
    "4" = "OPCs|aNSCs",
    "5" = "Neurons",
    "6" = "Pericytes",
    "7" = "Neurons"
  )) 
brain$cell_type <- Idents(brain)
Idents(brain) <- brain$seurat_clusters

# Dim plot of UMAP
DimPlot(brain, group.by = "seurat_clusters") + DimPlot(brain, group.by = "cell_type")

## Goal #1: Find the gene sets and pathways that distinguish Cluster 0 ##
brain.markers <- FindAllMarkers(brain, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.5)
sig_markers <- brain.markers %>% 
  filter(p_val_adj < .05) %>%
  write_csv(file = "brain_sig_markers.csv")

# Cluster 0 enrichr: https://maayanlab.cloud/Enrichr/enrich?dataset=b33d5997396c6721c3bf4bf0e67b8fbe
to_check <- c("Human_Gene_Atlas", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
              "KEGG_2019_Human", "MSigDB_Hallmark_2020")
# Pull the Cluster 0 genes as a vector
cluster0_genes <- sig_markers %>%
  filter(cluster == "0") %>%
  pull(gene)
# Run through enrichr
cluster0_eresList <- enrichr(cluster0_genes, databases = to_check)
# Get the top hits from each and plot them
plotList <- lapply(names(cluster0_eresList), plot_eres, eres_list = cluster0_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 0 Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster0_brain_enrichr_plots.png", width = 20, height = 10)
 
## Goal #2: Find the genes which differentiate clusters 7 and 5 ##
cluster7vs5.markers <- FindMarkers(brain,
                                   ident.1 = 7, 
                                   ident.2 = 5,
                                   logfc.threshold = .25,
                                   min.pct = 0.25)
# Pull the Cluster 7 vs 5 genes as a vector
cluster7_genes <- cluster7vs5.markers %>%
  filter(avg_log2FC > 0 & p_val_adj < .05) %>%
  rownames_to_column("gene") %>%
  pull(gene)
# Run through enrichr
cluster7_eresList <- enrichr(cluster7_genes, databases = to_check)
# Get the top hits from each and plot them
plotList <- lapply(names(cluster7_eresList), plot_eres, eres_list = cluster7_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 7 vs 5 Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster7_vs5_enrichr_plots.png", width = 20, height = 10)

### Part #1: Data Integration ###

## Data integration in EWS Cell Lines ##

# Get Ewing Sarcoma Cell Lines CHLA9 and CHLA10 
# From: https://www.mdpi.com/2072-6694/12/4/948
CHLA10_mat <- read.csv("CHLA10_matrix.csv.gz", row.names = 1)
CHLA9_mat <- read.csv("CHLA9_matrix.csv.gz", row.names = 1)
cells <- merge(x = CreateSeuratObject(CHLA9_mat, project = "CHLA9"),
                  y = CreateSeuratObject(CHLA10_mat, project = "CHLA10")) %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(nFeature_RNA > 2500 & nCount_RNA > 12000 & percent.mt < 18) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunUMAP(dims = 1:50)

# Dim Plots
DimPlot(cells, group.by = "orig.ident") 
DimPlot(cells, group.by = "orig.ident") + DimPlot(cells, group.by = "seurat_clusters") 

# Compare Batch and Cluster ID
compTable <- table(cells$orig.ident, cells$seurat_clusters)
compTable <- (compTable / rowSums(compTable)) * 100
pheatmap(compTable, color = WtOrRd_pal)

# Split the seurat object and integrate with CCA
cellsList <- SplitObject(cells, split.by = "orig.ident")

# Normalize and identify variable features for each dataset independently
cellsList <- lapply(X = cellsList, SCTransform) 
features <- SelectIntegrationFeatures(object.list = cellsList, 
                                      nfeatures = 3000)
cellsList <- PrepSCTIntegration(object.list = cellsList, 
                                   anchor.features = features)
cells.anchors <- FindIntegrationAnchors(object.list = cellsList,
                                           normalization.method = "SCT",
                                           anchor.features = features)
cells.int <- IntegrateData(anchorset = cells.anchors, 
                              normalization.method = "SCT")

# Re-run the pipeline on the integrated data
cells.int <- cells.int %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30) 

# Dim Plots
DimPlot(cells.int, group.by = "orig.ident") + DimPlot(cells, group.by = "orig.ident")
DimPlot(cells.int, group.by = "orig.ident") + DimPlot(cells.int, group.by = "seurat_clusters")

# Compare Batch and Cluster ID
compTable.int <- table(cells.int$orig.ident, cells.int$seurat_clusters)
compTable.int <- (compTable.int / rowSums(compTable.int)) * 100
pheatmap(compTable.int, color = WtOrRd_pal)

## Activity #1: Cell cycle and Integration on PDX batches ##

# Get the PDX data
filesNow <- list.files("Delattre_EWS_PDX", pattern = "mat\\.csv\\.gz", full.names = TRUE)
pdxIDs <- gsub(pattern = ".+/(.+)_mat\\.csv\\.gz", x = filesNow, replacement = "\\1")
names(filesNow) <- pdxIDs
srtList <- lapply(names(filesNow), FUN = function(pdx_now) {
  mat_now <- read.csv(filesNow[pdx_now], row.names = 1)
  CreateSeuratObject(mat_now, project = pdx_now)
})
pdx <- merge(x = srtList[[1]], y = srtList[-1]) %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(nFeature_RNA > 800 & nCount_RNA > 4000 & 
           percent.mt < 20 & nCount_RNA < 60000) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

# Dim Plot to see the batches
DimPlot(pdx, group.by = "orig.ident")
DimPlot(pdx, group.by = "orig.ident") + DimPlot(pdx, group.by = "seurat_clusters")

# Compare Batch and Cluster ID
compTable <- table(pdx$orig.ident, pdx$seurat_clusters)
compTable <- (compTable / rowSums(compTable)) * 100
pheatmap(compTable, color = WtOrRd_pal)

# 1. Integrate the 5 batches together

# 2. Run the pipeline on the integrated dataset

# 3. Plot the Integrated and Original datasets side-by-side in UMAP

# 4. Compare Batch and Cluster ID

### Part #2: Module Scoring ###

## Cell cycle scoring ##
# Get the cell cycle scores
cells.int <- CellCycleScoring(cells.int, search = TRUE,
                              s.features = cc.genes.updated.2019$s.genes,
                              g2m.features = cc.genes.updated.2019$g2m.genes)
# Plot in UMAP
DimPlot(cells.int, group.by = "Phase") + DimPlot(cells.int, group.by = "orig.ident")
# Plot with clusters
compTable.int <- table(cells.int$orig.ident, cells.int$Phase)
compTable.int <- (compTable.int / rowSums(compTable.int)) * 100
pheatmap(compTable.int, color = WtOrRd_pal, breaks = 1:100)
# FeaturePlot for showing mod score
FeaturePlot(cells.int, features = c("S.Score", "G2M.Score"))

# Module Scoring and Comparison:

# 1. Get the "HALLMARK_HYPOXIA" gene set
h_gene_sets <-msigdbr(category = c("H"))
g2m_genes <- h_gene_sets %>%
  filter(gs_name == "HALLMARK_G2M_CHECKPOINT") %>%
  pull(gene_symbol)

# 2. Calculate the module score
cells.int <- AddModuleScore(cells.int, name = "g2m_score", assay = "RNA", 
                          features = list(g2m_genes),
                          search = TRUE)
FeaturePlot(cells.int, features = "g2m_score1") 

# 3. Compare to clusters using VlnPlot and DimPlot
DimPlot(cells.int, group.by = "seurat_clusters") + VlnPlot(cells.int, features = "g2m_score1")

# 4. Scatter plot comparing G2M module score with the G2M score found by CellCycleScoring
FeatureScatter(cells.int, feature1 = "g2m_score1", 
               feature2 = "G2M.Score",
               group.by = "Phase")

## Activity #2: Score Hypoxia and Cell Cycle in PDX ##

# Cell cycle scoring:

# 1. Do the Cell Cycle Scoring

# 2. Plot Cell Cycle Phases in UMAP

# 3. Plot Phases and clusters in a heatmap
compTable.int <- table(_____, ______)
compTable.int <- (compTable.int / rowSums(compTable.int)) * 100
_______(compTable.int, color = WtOrRd_pal, breaks = 1:100)

# Module Scoring and Comparison:

# 1. Get the "HALLMARK_HYPOXIA" gene set
h_gene_sets <-msigdbr(category = c("H"))
hypoxia_genes <- h_gene_sets %>%
  filter(gs_name == "HALLMARK_HYPOXIA") %>%
  pull(gene_symbol)

# 2. Calculate the Hypoxia module score

# 3. Compare to clusters for Hypoxia Score using VlnPlot

# Do the same, but with Metastasis:
# 1. Get the "TAVAZOIE_METASTASIS" gene set
cp_gene_sets <-msigdbr(category = c("C2"))
emt_genes <- cp_gene_sets %>%
  filter(gs_name == "TAVAZOIE_METASTASIS") %>%
  pull(gene_symbol)

# 2. Calculate the Metastasis module score

# 3. Compare to clusters for Metastasis Score using VlnPlot

# 4. Bonus: Do a feature scatter to compare hypoxia and Metastasis scores in cluster 2 and 5

### Part #3: Psuedotime ###

## Pseudotime example ##

# Convert to monocle3 cell_data_set
cells.cds <- as.cell_data_set(cells.int)

# Re-cluster (required for trajectories in monocle3)
cells.cds <- cluster_cells(cells.cds)

# Learn the trajectory
cells.cds <- learn_graph(cells.cds)

# Plot the trajectory
plot_cells(cells.cds, color_cells_by = "seurat_clusters", show_trajectory_graph = TRUE)

# Find the pseudotime order
DimPlot(cells.int, group.by = "Phase")
cells.cds <- order_cells(cells.cds)

# Plot pseudotime
plot_cells(cells.cds, color_cells_by = "pseudotime")

# Add pseudotime score back into the dataset
cells.int$pseudotime_score <- cells.cds@principal_graph_aux@listData$UMAP$pseudotime
cells.int$pseudotime_score[cells.int$pseudotime_score == Inf] <- NA

# FeaturePlot with Pseudotime
FeaturePlot(cells.int, features = c("pseudotime_score"))

## Activity 3 (optional): Pseudotime analysis ##

# Convert to monocle3 cell_data_set
pdx.cds <- as.cell_data_set(pdx.int)

# Re-cluster (required for trajectories in monocle3)
pdx.cds <- cluster_cells(pdx.cds)

# Learn the trajectory
pdx.cds <- learn_graph(pdx.cds)

# Plot the trajectory
plot_cells(pdx.cds, color_cells_by = "seurat_clusters", show_trajectory_graph = TRUE)

# Find the pseudotime order
DimPlot(pdx.int, group.by = "Phase")
pdx.cds <- order_pdx(pdx.cds)

# Plot pseudotime
plot_pdx(pdx.cds, color_pdx_by = "pseudotime")

# Add pseudotime score back into the dataset
pdx.int$pseudotime_score <- pdx.cds@principal_graph_aux@listData$UMAP$pseudotime
pdx.int$pseudotime_score[pdx.int$pseudotime_score == Inf] <- NA

# FeaturePlot with Pseudotime
FeaturePlot(pdx.int, features = c("pseudotime_score"))


### Next time: Optional Project assignments and Final discussion of missed topics ###
### Send me your requests for any topics I didn't cover that you wanted ###



