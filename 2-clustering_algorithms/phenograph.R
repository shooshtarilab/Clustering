#devtools::install_github("JinmiaoChenLab/Rphenograph")
library("Rphenograph")
library(SingleCellExperiment)

main <- function(filename) {
  
  load(filename)
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/phenograph.csv")
  
  
  data <- t(counts(sce))
  start = Sys.time()
  Rphenograph_out <- Rphenograph(data, k = 45)
  #modularity(Rphenograph_out[[2]])
  membership(Rphenograph_out[[2]])
  
  end = Sys.time()
  time = end - start
  res <- data.frame(cell=colnames(sce), phenograph=factor(membership(Rphenograph_out[[2]])))
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " phenograph"), file="time.txt", append=TRUE, sep = "\n")
  
}

datasets <- c("lambrechts_lung_nonTumour_TPM.RData")

for (i in datasets) {
  main(i)
}

