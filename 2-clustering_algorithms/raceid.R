library(SingleCellExperiment)
library(RaceID)

main <- function(filename) {
  
  print(filename)
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/raceid.csv")
  
  load(filename)
  
  mat <- data.frame(counts(sce))
  start = Sys.time()
  sc <- SCseq(mat)
  #sc <- SCseq(as.data.frame(as.matrix(counts(sce))))
  sc <- filterdata(sc, mintotal=1, minexpr = 1, minnumber = 5,
                   LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                   bmode = "RaceID")
  
  sc <- compdist(sc, metric="euclidean", FSelect = FALSE, knn = NULL)
  sc <- clustexp(sc, sat = TRUE, samp = NULL, cln = NULL, clustnr = 30,
                 bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
  sc <- findoutliers(sc, probthr = 0.001, outminc = 5, outlg = 2,
                     outdistquant = 0.95)
  
  res <- data.frame(colnames(mat), sc@cpart)
  res$RaceID3 <- sc@cpart
  length(unique(res$RaceID3))
  colnames(res) <- c("cell", "raceid")
  
  end = Sys.time()
  
  time = end - start
  write.csv(res, file=folder, quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " raceid"), file="time.txt", append=TRUE, sep = "\n")
  
}

