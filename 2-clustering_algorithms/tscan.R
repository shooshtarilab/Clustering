#BiocManager::install("TSCAN")
#source("https://bioconductor.org/biocLite.R")
#biocLite("TSCAN", ref="1.19.0", dependencies=F)
# library(TSCAN)
# library(SingleCellExperiment)

suppressPackageStartupMessages(library(TSCAN))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(optparse))

main <- function(filename) {

  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2])
  folder = paste0("./", dataset, "/tscan.csv")
  
  load(filename)
  
  start = Sys.time()

  procdata <- preprocess(as.matrix(counts(sce)), minexpr_percent = 0.01, 
                         clusternum = NULL, takelog = FALSE,
                         pseudocount = 1, minexpr_value = 1,
                         cvcutoff = 1)
  set.seed(4343646)
  lpsmclust <- exprmclust(procdata, clusternum = 2:20, modelNames = "VVV", reduce = T)
  res <- data.frame(row.names(colData(sce)), lpsmclust$clusterid)
  
  end = Sys.time()
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = TRUE, col.names=TRUE)
  #cat(paste0(time, " ", dataset, " tscan"), file="time.txt", append=TRUE, sep = "\n")

}

datasets <- c("li_crc_counts.RData")

for (i in datasets) {
  main(i)
}


