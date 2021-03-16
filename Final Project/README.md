# Final Project

This project will get you experienced with a real-world dataset while you uncover
metastatic subpopulations within triple-negative breast cancers. Please answer the
following questions and format your results either in RMarkdown (bonus points) or
in a PowerPoint presentation -- then send to Henry (millerh1@livemail.uthscsa.edu).

Data: `counts_rsem.csv.gz`

Meta Data: `meta_data.csv`

Metastasis Genesets: `metastasis_genesets/`

Data Sources (including gene sets): [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118389), [Paper](https://www.nature.com/articles/s41467-018-06052-0), and [GitHub](https://github.com/Michorlab/tnbc_scrnaseq)

## Questions to answer (answer these using figures / tables generated in analysis):

1. What quality issues were present in the data? How did you deal with them?
2. What cell types exist in the dataset? How did you identify the tumor cells?
3. What biologically-relevant tumor subpopulations exist?
4. Were you able to identify metastatic subpopulations? How?
5. Is there a correlation between metastasis signature and any other genes / processes...?

## Suggested steps -- not required:

1. Download the dataset and import into `Seurat`
2. Perform downstream QC steps to filter the data
3. Do the steps necessary for getting a UMAP of the data and clusters
4. Assign cell types to each cell and identify the tumor cells
5. Do a cluster marker analysis to identify biologically-meaningful tumor subpopulations
6. Use Module Scoring to identify metastatic cells. 
7. Find the correlation between metastasis and other genes (might use MAGIC first).
