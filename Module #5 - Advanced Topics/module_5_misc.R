# Code to get the PDX data
download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130025/suppl/GSE130025_RAW.tar", 
              destfile = "Delattre_EWS_PDX.tar", method = "curl")
untar("Delattre_EWS_PDX.tar", exdir = "Delattre_EWS_PDX")
rmFiles <- list.files("Delattre_EWS_PDX/", pattern = "\\.txt|\\.bed", full.names = T)
file.remove(rmFiles)
filesNow <- list.files("Delattre_EWS_PDX", pattern = "mtx", full.names = T)
pdxIDs <- gsub(pattern = ".+_(.+)_matrix.+", x = filesNow, replacement = "\\1")
names(filesNow) <- pdxIDs
srtList <- lapply(names(filesNow), FUN = function(idNow, filesNow) {
  mtxNow <- filesNow[idNow]
  barcodesNow <- gsub(x = mtxNow, pattern = "matrix.mtx", replacement = "barcodes.tsv")
  genesNow <- gsub(x = mtxNow, pattern = "matrix.mtx", replacement = "genes.tsv")
  sm <- Matrix::readMM(mtxNow)
  genesNow <- read.table(genesNow)
  rownames(sm) <- gsub(genesNow$V2, pattern = "GRCh38_(.+)", replacement = "\\1")
  barcodesNow <- read.table(barcodesNow)
  colnames(sm) <- barcodesNow$V1
  sm <- sm[, sample(colnames(sm), 600)]
  write.csv(as.data.frame(sm), file = paste0("Delattre_EWS_PDX/", idNow, "_mat.csv"))
  srtTemp <- CreateSeuratObject(sm, project = idNow)
  return(srtTemp)
}, filesNow = filesNow)
DELATTRE <- merge(x = srtList[[1]], y = c(srtList[c(2:length(srtList))]))

# Code to get the cell line data
if (! file.exists("EWS_Cells_10X/GSM4368462_CHLA9_barcodes.tsv.gz")) {
  archive_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE146221&format=file"
  download.file(archive_url, destfile = "EWS_Cells.tar", method = "curl")  # Method = curl to avoid corruption
  untar("EWS_Cells.tar", exdir = "EWS_Cells_10X") 
}
cell_types <- c("CHLA9", "CHLA10")
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
  
  sm <- sm[, sample(colnames(sm), 1200)]
  write.csv(as.data.frame(sm), file = paste0(cell_now, "_matrix.csv"))
  
  srt <- CreateSeuratObject(sm, project = cell_now)
  return(srt)
})
full_srt <- merge(x = srtList[[1]], y = srtList[-1])

