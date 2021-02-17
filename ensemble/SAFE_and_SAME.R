#INSTALL PACKAGES
#devtools::install_github("yycunc/SAMEclustering")
#devtools::install_github("yycunc/SAFEclustering")

#LOAD PACKAGES
library("SAFEclustering")
library("SAMEclustering")

run_SAFE <- function(input_dataframe, tools, dir) {
  final <- input_dataframe
  res_SAFE <- data.frame(cell = final$cell)
  final$cell=NULL
  
  final <- t(final)
  dims <- dim(final)
  safe_input <- t(apply(final, 1, as.numeric))
  
  cluster.SAFE <- SAFE(cluster_results = safe_input, program.dir = "/Users/alainamahalanabis/Documents/clustering", 
                       MCLA = TRUE, CSPA = FALSE, HGPA = TRUE, SEED = 123)
  

  res_SAFE <- cbind(res_SAFE, cluster = cluster.SAFE$optimal_clustering)

  colnames(res_SAFE)[ncol(res_SAFE)] <- "SAFE_5"
  dir = paste0(dir, "SAFE_5.csv")
  
  write.csv(res_SAFE, file = dir, row.names = FALSE)
  
}

run_SAME <- function(input_dataframe, tools, dir) {
  
  final <- input_dataframe
  res_SAME <- data.frame(cell = final$cell)
  final$cell=NULL
  
  final <- t(final)
  dims <- dim(final)
  safe_input <- t(apply(final, 1, as.numeric))
  
  same_input <- t(safe_input)
  rownames(same_input) = NULL
  colnames(same_input) = NULL
  
  cluster.SAME <- SAMEclustering(Y = same_input, rep = 3, SEED = 123)
  
  res_SAME <- data.frame(cell=res_SAME$cell, SAME_AIC_3 = cluster.SAME$AICcluster ,SAME_BIC_3 = cluster.SAME$BICcluster)
  
  filename1=paste(dir, "SAME_AIC_3.csv")
  filename2=paste(dir, "SAME_BIC_3.csv")
  
  write.csv(data.frame(cell=res_SAME$cell, SAME_AIC_3 = res_SAME$SAME_AIC), file = filename1, row.names = FALSE)
  write.csv(data.frame(cell = res_SAME$cell, SAME_BIC_3 = res_SAME$SAME_BIC), file = filename2, row.names = FALSE)
  
}

all.datasets <- c("li_crc_nonTumour", "tirosh_melanoma_nonTumour", "chung_breast_nonTumour")

for (dataset in all.datasets) {
  filename = paste0("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/Clustering_Results/", dataset, ".csv")
  final <- read.csv(filename)
  
  print("DATSET")
  print (dataset)
  
  cols <- intersect(colnames(final), c("cell", "cidr", "monocle", "bigScale", "raceid", "rca"))
  final <- final[,cols]
  for (col in colnames(final)) {
    if (col != "cell") {
      final <- transform(final,col=as.numeric(factor(final[,col])))
      final[,col] <- final[,"col"]
      final[,"col"] = NULL
      final[,col][final[,col] == 0] <- max(unique(final[,col])) + 1
    }
  }
  dir = paste0("/Users/alainamahalanabis/Documents/clustering/GEO_datasets/", dataset, "/")
  run_SAFE(final, cols, dir)
}
