#### Module #4: Biological Interpretation #### 

### Part #0: Recap ###

# Check for missing R packages, install if needed
list_of_pkgs <- c("hdf5r", "tidyverse", "ggpubr", "Seurat", "enrichR", "cowplot", "patchwork", "VennDiagram")
install.packages(list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())])

# Check for missing BioC packages, install if needed
list_of_pkgs <- c("limma", "org.Hs.eg.db", "celldex", "SingleR")
bioc2install <- list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())]
if (length(bioc2install) > 0) {
  BiocManager::install(bioc2install)
}


## Load libraries ##
library(Seurat)
library(tidyverse)
library(ggpubr)
library(limma)
library(celldex)
library(SingleR)
library(enrichR)
library(cowplot)
library(patchwork)
library(VennDiagram)

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
  RunUMAP(dims = 1:10)


### Part #1: Cluster Marker Identification and cell type classification ###

## Part 1a: Finding Cluster Markers ##

# UMAP for PBMC
DimPlot(pbmc, reduction = "umap")

# Find markers for cluster 3
cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, 
                                only.pos = TRUE,
                                logfc.threshold = .5,
                                min.pct = 0.25)
cluster3.markers %>%
  top_n(5, -log10(p_val))

# Feature plot for top 2 Cluster 3 markers
FeaturePlot(pbmc, features = c("MS4A1", "CD79A"))
# Violin top 2 Plot
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# Find markers of 0, 2, 4, and 6
cluster0246.markers <- FindMarkers(pbmc, ident.1 = c(0, 2, 4, 6),
                                   logfc.threshold = .5, only.pos = TRUE,
                                   min.pct = 0.25)
cluster0246.markers %>%
  top_n(5, -log10(p_val))
# Feature plot for 0, 2, 4, 6 markers
FeaturePlot(pbmc, features = c("IL32", "CD3D"))
# Violin Plot
VlnPlot(pbmc, features = c("IL32", "CD3D"))


## Activity 1a: Find Cluster Markers in Brain Dataset ##
# Data comes from here: https://www.nature.com/articles/s41586-018-0590-4
# FigShare: https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687

# # Processed the data to make it smaller...
# brain.raw <- read.csv("Brain_Non-Myeloid-counts.csv.gz", row.names = 1)
# set.seed(42); brain.raw <- brain.raw[,sample(colnames(brain.raw), 2400)]
# write.csv(brain.raw, file = "Brain_Non-Myeloid-counts.small.csv")

brain.raw <- read.csv("Brain_Non-Myeloid-counts.small.csv.gz", row.names = 1)
brain <- Seurat::CreateSeuratObject(brain.raw) %>%
  subset(subset = nFeature_RNA > 600 & nFeature_RNA < 5000 &
           nCount_RNA < 2.5E6) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution= .5) %>%
  RunUMAP(dims = 1:20)


# UMAP of brain
DimPlot(brain, reduction = "umap")

# Find the markers for cluster 0

# Feature plot for top 2 Cluster 0 markers

# Violin Plot for top 2 Cluster 0 markers

# What kind of cells are likely in this cluster?

## Part 1b: Manual Cell Type Classification from Prior Knowledge ##

# Marker list from: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
marker_list <- list(
  'Naive CD4+ T' = c("IL7R", "CCR7"),
  'Memory CD4+' = c('IL7R', 'S100A4'),
  'B Cells' = c("MS4A1", "CD79A"),
  "CD8+ T" = c("CD8A", "GZMK"),
  'CD16+ Mono' = c("FCGR3A", "MS4A7"),
  'NK Cells' = c("GNLY", "NKG7"),
  "Dendritic Cells" = c("FCER1A", "CST3"),
  'CD14+ Mono' = c("CD14", "LYZ"),
  "Platelets" = c("PPBP")
)

# Use FeaturePlots/VlnPlots to find...
# B Cells
FeaturePlot(pbmc, features = marker_list$`B Cells`)
VlnPlot(pbmc, features = marker_list$`B Cells`)
# Naive CD4+ T
FeaturePlot(pbmc, features = marker_list$`Naive CD4+ T`)
VlnPlot(pbmc, features = marker_list$`Naive CD4+ T`)
# CD14+ Mono
FeaturePlot(pbmc, features = marker_list$`CD16+ Mono`)
VlnPlot(pbmc, features = marker_list$`CD16+ Mono`)
# NK Cells
FeaturePlot(pbmc, features = marker_list$`NK Cells`)
VlnPlot(pbmc, features = marker_list$`NK Cells`)
# Memory CD4+
FeaturePlot(pbmc, features = marker_list$`Memory CD4+`)
VlnPlot(pbmc, features = marker_list$`Memory CD4+`)
# CD8+ T
FeaturePlot(pbmc, features = marker_list$`CD8+ T`)
VlnPlot(pbmc, features = marker_list$`CD8+ T`)
# CD16+ Mono
FeaturePlot(pbmc, features = marker_list$`CD16+ Mono`)
VlnPlot(pbmc, features = marker_list$`CD16+ Mono`)
# Dendritic Cells
FeaturePlot(pbmc, features = marker_list$`Dendritic Cells`)
VlnPlot(pbmc, features = marker_list$`Dendritic Cells`) 

# Do all VlnPlots at once automatically and save
marker_vln <- function(srt, marker_list, marker_now) {
  vln <- VlnPlot(srt, features = marker_list[[marker_now]]) +
    patchwork::plot_annotation(title = marker_now,
                               theme = theme(title = element_text(size = 22)))
  return(vln)
}
vlnList <- lapply(names(marker_list), marker_vln, 
                  srt = pbmc, marker_list = marker_list)
ga <- ggarrange(plotlist = vlnList)
ggsave(ga, filename = "arranged_pbmc_violins.png", height = 10, width = 20)


## Activity 1b: Prior knowledge manual analysis of brain clusters ##
# Markers derived from http://biocc.hrbmu.edu.cn/CellMarker/index.jsp
# and also https://tabula-muris-senis.ds.czbiohub.org/brain-non-myeloid/facs/
marker_list_brain <- list(
  'Oligodendrocytes' = c("Plp1", "Mag", "Cnp"),
  'OPCs' = c('Pdgfra', 'C1ql1', 'Cntn1'),
  'Endothelial Cells' = c("Ly6c1", "Esam", "Pecam1"),
  "Neurons" = c("Meg3", "Snap25", "Syt1"),
  'Astrocytes' = c("Gfap", "Sox9", "Slc1a2"),
  'Pericyte' = c("Mcam", "Pdgfrb"),
  "Neuroepithelial Cells" = c("Krt5", "Krt15")
)

# Using FeaturePlots and/or VlnPlots and find each population from the marker list
# Note: You can use the lapply function here as well...


## Part 1c: Automated classification using celldex and singleR ##
# Predict using the MonacoImmuneData from celldex
mimd.se <- MonacoImmuneData()
# What is the structure of mimd.se?
mimd.se
View(mimd.se)

# Predict using the Broad labels
sceP <- as.SingleCellExperiment(pbmc)
pred.mimd <- SingleR(test = sceP, 
                     ref = mimd.se, 
                     labels = mimd.se$label.main)
plotScoreHeatmap(pred.mimd,
                 max.labels = 8,
                 clusters = pbmc$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)

# Predict using the Fine-grained labels
pred.mimd.fine <- SingleR(test = sceP, 
                          ref = mimd.se, 
                          labels = mimd.se$label.fine)
plotScoreHeatmap(pred.mimd.fine,
                 max.labels = 14,
                 clusters = pbmc$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)

# Combine the labels for improved accuracy
pred.comb <- combineCommonResults(
  list(
    "Broad" = pred.mimd,
    "Fine" = pred.mimd.fine
  )
)
plotScoreHeatmap(pred.comb,
                 max.labels = 14,
                 clusters = pbmc$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)

# Add the predicted labels and compare to clusters
# You are looking for clusters which have unambiguous labels
pbmc$predicted_id <- pred.comb$pruned.labels
table(pbmc$predicted_id, pbmc$seurat_clusters)

# Select final IDs
new_ids <- c(
  "0" = "CD4 + T Cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD4 + T Cells",
  "3" = "B Cells",
  "4" = "CD8+ T Cells",
  "5" = "CD16+ Monocytes",
  "6" = "NK Cells",
  "7" = "Dendritic Cells",
  "8" = "Platelets"
)
pbmc <- RenameIdents(pbmc, new_ids)
pbmc$cell_type <- Idents(pbmc)
DimPlot(pbmc)
DimPlot(pbmc, label.size = 5, label = TRUE) + NoLegend() 


## Activity 1c: Do automated cell type classification on the Brain dataset ##
# Do automated cell type classification using SingleR 
# Then, select final labels, add them to the Seurat object, and plot the UMAP
mouse_cells <- celldex::MouseRNAseqData()

# Use the mouse_cells reference with SingleR with broad labels

# Plot the resulting heatmap

# Use the mouse_cells reference with SingleR with fine labels

# Plot the heatmap of cluster assignment results

# Combine the results and plot the combined heatmap

# Make the cell type assignments on the Seurat object

# Plot the UMAP with cell types labeled


## Part 1 Bonus: Semi-automated Cell Type Classification from Prior Knowledge ##
# Find all markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.5)
marker_summary <- pbmc.markers %>% 
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC)
# Violin plot to discover cluster identity
VlnPlot(pbmc, features = marker_summary$gene, ncol = 3)

# Get the cell marker list for PBMCs
# http://biocc.hrbmu.edu.cn/CellMarker/search.jsp?quickSearchInfo=Blood
pbmc_markerlist <- read_csv("CellMarker_PBMC.csv")
# Clean the list and compile markers from scRNA-Seq experiments
pbmc_markerlist <- pbmc_markerlist %>%
  filter(Species == "Human" & Cancer == "Normal") %>%
  select(`Cell Type`, `Cell Marker`) %>%
  separate_rows(`Cell Marker`, sep = ", ") %>%
  mutate(gene = alias2SymbolTable(`Cell Marker`)) %>%
  filter(! is.na(gene)) %>%
  select( -`Cell Marker`)

# Join together with the results of FindMarkers
pbmc_markresults <- pbmc.markers %>% 
  inner_join(pbmc_markerlist, by = "gene") %>%
  group_by(gene) %>%
  mutate(cellTypes = paste0(unique(`Cell Type`), collapse = ", ")) %>%
  select(-`Cell Type`) %>%
  dplyr::distinct() %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = -log10(p_val)) 


### Part #2: Cluster Marker Analysis ###

# First, re-assign idents as the clusters
Idents(pbmc)
Idents(pbmc) <- pbmc$seurat_clusters
Idents(pbmc)

# Find all markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.5)
sig_markers <- pbmc.markers %>% 
  filter(p_val_adj < .05)

## Part 2a: Pathway enrichment online in enrichr: https://maayanlab.cloud/Enrichr/ ##
write_csv(sig_markers, file = "pbmc_sig_markers.csv")

# Cluster 0 enrichr: https://maayanlab.cloud/Enrichr/enrich?dataset=b33d5997396c6721c3bf4bf0e67b8fbe
dbs <- listEnrichrDbs()
to_check <- c("Human_Gene_Atlas", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
              "KEGG_2019_Human", "MSigDB_Hallmark_2020")

# Pull the Cluster 0 genes as a vector
cluster0_genes <- sig_markers %>%
  filter(cluster == "0") %>%
  pull(gene)
# Run through enrichr
cluster0_eresList <- enrichr(cluster0_genes, databases = to_check)
# Plotting function for cluster enrichment results
plot_eres <- function(eres_name, eres_list, n = 10) {
  eres <- eres_list[[eres_name]]
  eres %>%
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
# Get the top hits from each and plot them
plotList <- lapply(names(cluster0_eresList), plot_eres, eres_list = cluster0_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 0 Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster0_enrichr_plots.png", width = 20, height = 10)

## Part 2b: Finding differences using genes ##
DimPlot(pbmc) + DimPlot(pbmc, group.by = "cell_type")

# Find differences in marker genes
cluster_gene_list <- lapply(c("0", "2"), function(cluster_now) {
  sig_markers %>%
    filter(cluster == cluster_now) %>%
    pull(gene)
})
names(cluster_gene_list) <- c("Cluster: 0", "Cluster: 2")

# Use Venn Diagram to Compare
venn.diagram(cluster_gene_list, filename = "clusters_0_2_compared.png", 
             fill = c("firebrick", "skyblue"), margin = .05)
overlap_list <- calculate.overlap(cluster_gene_list)

## Part 2c: Use a differential model to find differences in biology ##
# Find markers in Cluster 0 VS Cluster 2
cluster0vs2.markers <- FindMarkers(pbmc,
                                   ident.1 = 0, 
                                   ident.2 = 2,
                                   logfc.threshold = .25,
                                   min.pct = 0.25)

# Pull the Cluster 0 vs 2 genes as a vector
cluster0_genes <- cluster0vs2.markers %>%
  filter(avg_log2FC > 0) %>%
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
  ggsave(filename = "cluster0_vs2_enrichr_plots.png", width = 20, height = 10)

# Pull the Cluster 2 vs 0 genes as a vector
cluster2_genes <- cluster0vs2.markers %>%
  filter(avg_log2FC < 0) %>%
  rownames_to_column("gene") %>%
  pull(gene)
# Run through enrichr
cluster2_eresList <- enrichr(cluster2_genes, databases = to_check)
# Get the top hits from each and plot them
plotList <- lapply(names(cluster2_eresList), plot_eres, eres_list = cluster2_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 2 vs 0 Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster2_vs0_enrichr_plots.png", width = 20, height = 10)

## Activity do 2a, 2b, and 2c on the brain dataset ##
# Goal #1: Find the gene sets and pathways that distinguish Cluster 0
# Goal #2: Find the genes which differentiate clusters 5 and 7


