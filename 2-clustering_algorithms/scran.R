
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dynamicTreeCut))

pre_clean <- function(sce) {
  
  ave.counts <- rowMeans(as.matrix(counts(sce)))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  
  return(sce1)
}

main <- function(filename) {
  
  # option_list <- list(
  #   make_option(c("-i", "--input"), default="NA",
  #               help="A RData file containing a SingleCellExperiment object"),
  #   #
  #   make_option(c("-o", "--outdir"), default="NA",
  #               help="A path/name for the results directory"),
  # )
  # 
  # Input<- opt$input
  # Outdir<- opt$outdir
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/scran.csv")
  
  load(filename)
  
  sce1 <- pre_clean(sce)
  
  start = Sys.time()

  clusters <- scran::quickCluster(counts(sce1), use.ranks=FALSE)
  sce1 <- scran::computeSumFactors(sce1, cluster=clusters)
  sce1 <- normalize(sce1)
  my.dist <- dist(t(counts(sce1)))
  set.seed(13115645)
  my.tree <- hclust(my.dist, method = "ward.D2")
  my.clusters <- unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0))
  
  end = Sys.time()
  
  res <- data.frame(row.names(colData(sce)), my.clusters)

  time = end - start
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = TRUE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " scran"), file="time.txt", append=TRUE, sep = "\n")

}
