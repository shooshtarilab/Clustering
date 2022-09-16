suppressPackageStartupMessages(library(CountClust))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(optparse))

library(CountClust)

main <- function(filename) {

  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/countclust.csv")
  
  load(filename)
  
  start = Sys.time()

  input <- as.matrix(counts(sce))
  outs <- FitGoM(t(input), K = 8, tol = 0.1, path_rda = NULL)
  
  end = Sys.time()
  #colData(sce)$countClust <- unlist(apply(outs$fit$omega, 1,  which.max))
  res <- colData(sce)
  
  res<-data.frame(cell = row.names(colData(sce)), countClust = unlist(apply(outs$fit$omega, 1,  which.max)))
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = TRUE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " countClust"), file="time.txt", append=TRUE, sep = "\n")

}