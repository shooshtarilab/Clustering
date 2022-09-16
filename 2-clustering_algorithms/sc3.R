suppressPackageStartupMessages(library(rngtools))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(SC3))

pre_clean <- function(sce) {
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  sce1 <- calculateQCMetrics(sce1)
  sce1 <- scater::normalize(sce1)
  
  return(sce1)
}

main <- function(filename) {
  
  load(filename)
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/sc3.csv")
  #sce1 <- pre_clean(sce)
  sce1 <- pre_clean(sce)
  rowData(sce1)$feature_symbol <- rowData(sce1)$genes
  set.seed(4758747)
  start = Sys.time()
  
  sce1 <- sce1[!duplicated(rowData(sce1)$feature_symbol), ]
  #isSpike(sce1, "ERCC") <- FALSE
  
  #IF USING ALREADY NORMALIZED DATA
  #logcounts(sce1) <- assay(sce1)
  sce1 <- sc3_estimate_k(sce1)
  k_est <- sce1@metadata$sc3$k_estimation
  sce1 <- sc3(sce1, ks = k_est, biology = FALSE, n_cores=4,  k_estimator = FALSE, rand_seed=0, gene_filter = FALSE)
  
  test <- function() {
    sce1 <- sc3_run_svm(sce1, ks = 53)
  }
  
  end = Sys.time()
  eval(parse(text=paste0("colData(sce)$SC3 <- colData(sce1)$sc3_", k_est, "_clusters")))
  res <-data.frame(row.names(colData(sce)), colData(sce)$SC3)
  colnames(res) <- c("cell", "sc3")
  
  time = end - start
  write.csv(res, file=folder, quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " sc3"), file="time.txt", append=TRUE, sep = "\n")

}

