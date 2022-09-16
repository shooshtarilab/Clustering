# BiocManager::install("SIMLR")
# install.packages("RSpectra")
# library("Rspectra")
# library(SIMLR)
# library("RcppEigen")
library(SingleCellExperiment)

suppressPackageStartupMessages(library(SIMLR))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(RSpectra))
suppressPackageStartupMessages(library(RcppEigen))

main <- function(filename) {
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/simlr.csv")
  
  load(filename)
  
  start = Sys.time()
  
  set.seed(345454654)
  cluster_num <- SIMLR_Estimate_Number_of_Clusters(X = as.matrix(counts(sce)), NUMC = 2:20, cores.ratio = 1)
  cluster_num <- (min(which.min(cluster_num$K1), which.min(cluster_num$K2))+1)
  output <-  SIMLR_Large_Scale(X = as.matrix(counts(sce)), c = cluster_num, k = 10, kk = 100)
   
  end = Sys.time()
  time = end - start
  #colData(sce)$SIMLR <- output$y$cluster
  res <- data.frame(row.names(colData(sce)),output$y$cluster)
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = TRUE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " simlr"), file="time.txt", append=TRUE, sep = "\n")
  
}

datasets <- c("tirosh_melanoma_allCells_counts.RData")

for (i in datasets) {
  main(i)
}


