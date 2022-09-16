library(cidr)
library(SingleCellExperiment)

pre_clean <- function(sce) {
  
  ave.counts <- rowMeans(as.matrix(counts(sce)))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  
  return(sce1)
}

main <- function(filename) {
  
  load(filename)
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/cidr.csv")
  
  sce1 <- pre_clean(sce)
  #View(counts(sce1))
  
  start = Sys.time() 
  singleCell <- scDataConstructor(assay(sce), tagType = "raw")
  initial = colData(sce1)
  
  singleCell <- determineDropoutCandidates(singleCell, min1 = 3, min2 = 8,
  N = 2000, alpha = 0.1, fast = TRUE, zerosOnly = FALSE,
  bw_adjust = 1)

  singleCell <- wThreshold(singleCell, cutoff = 0.5, plotTornado = FALSE)
  singleCell <- scDissim(singleCell, correction = FALSE, threads = 0,
  useStepFunction = TRUE)
  singleCell <- scPCA(singleCell, plotPC = TRUE)
  singleCell <- nPC(singleCell)
  set.seed(89549585)

  singleCell <- scCluster(singleCell, n = NULL, nCluster = NULL,
  nPC = NULL, cMethod = "ward.D2")
  
  end = Sys.time()
  time = end-start

  res <- data.frame(row.names(colData(sce)), singleCell@clusters)
  colnames(res) <- c("cell", "cidr")
 
  write.csv(res, file=folder, quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " cidr"), file="time.txt", append=TRUE, sep = "\n")
  
}



