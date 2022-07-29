library(SingleCellExperiment)
library(hopach)
setwd("/Users/alainamahalanabis/Documents/clustering/GEO_datasets")
main <- function(filename) {
  
  load(filename)
  mat <- counts(sce)
  print("number of cells in rData file")
  print(ncol(mat))
  print (filename)
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/altAnalyze.csv")
  print(folder)
  start = Sys.time()
  
  vars<-apply(mat,1,var)
  subset<-vars>quantile(vars,(nrow(mat)-200)/nrow(mat))
  mat.subset<-mat[subset,]
  gnames.subset<-genes[subset,]
  array.hobj<-hopach(t(mat.subset),d="euclid")
  print(length(array.hobj$clustering$labels))
  print(ncol(mat))
  print(array.hobj$clustering$k)
  res <- data.frame(cell = colnames(mat), altAnalyze = array.hobj$clustering$labels)
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = FALSE, col.names=TRUE)
  
}

datasets <- c("peng_pancreatic_malignant_counts.RData", "vangalen_AML_malignant_counts.RData")

for (i in datasets) {
  main(i)
}



