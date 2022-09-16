#devtools::install_github("JustinaZ/pcaReduce")
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
library("pcaMethods")
library("pcaReduce")

main <- function (filename) {
  load(filename)
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/pcaReduce.csv")
  counts(sce) <- log2(counts(sce)+1)
  
  start = Sys.time()
  
  Output_M <- PCAreduce(counts(sce), nbt=100, q=30, method='M')
  
  end = Sys.time()
  res <- data.frame(cells = colData(sce), pcaReduce = res <- Output_M[[1]][,1])
  time = end - start
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = TRUE, col.names=TRUE, row.names = FALSE)
  cat(paste0(time, " ", dataset, " pcaReduce"), file="time.txt", append=TRUE, sep = "\n")
  
}