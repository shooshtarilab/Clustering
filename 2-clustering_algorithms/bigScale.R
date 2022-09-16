#devtools::install_github("iaconogi/bigSCale2")
library("bigSCale")
library("scater")
library(SingleCellExperiment)

pre_clean <- function(sce1) {
  ave.counts <- rowMeans(counts(sce1))
  keep <- ave.counts >= 1
  sce1 <- sce1[keep,]
  sce1 <- calculateQCMetrics(sce1)
  sce1 <- scater::normalize(sce1)
  
  return(sce1)
}

main <- function(filename) {
  
  load(filename)
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/bigScale.csv")
  start = Sys.time()
  
  sce1 <- sce
  sce1=bigscale(sce1,speed.preset='fast', modality = "pca", clustering='recursive')
  sce1=setDistances(sce1)
  sce1=setClusters(sce1)
  sce1=storeTransformed(sce1)
  end = Sys.time()
  time = end - start
  res <- data.frame(cell=colnames(sce1), bigScale=getClusters(sce1))
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " bigscale"), file="time.txt", append=TRUE, sep = "\n")
  
}


