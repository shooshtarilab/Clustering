#INSTALL PACKAGES
#devtools::install_github("yycunc/SAMEclustering")
#devtools::install_github("yycunc/SAFEclustering")

#LOAD PACKAGES
library("SAFEclustering")
library("SAMEclustering")

run_SAFE <- function(input_dataframe, tools) {
  final <- input_dataframe
  res_SAFE <- data.frame(cell = final$cell)
  final <- final[,tools]
  final$cell=NULL
  
  final <- t(final)
  dims <- dim(final)
  safe_input <- t(apply(final, 1, as.numeric))
  
  cluster.SAFE <- SAFE(cluster_results = safe_input, program.dir = "/Users/alainamahalanabis/Documents/clustering", 
                       MCLA = TRUE, CSPA = FALSE, HGPA = TRUE, SEED = 123)
  

  res_SAFE <- cbind(res_SAFE, cluster = cluster.SAFE$optimal_clustering)

  #colnames(res_SAFE)[ncol(res_SAFE)] <- "SAFE_3"
  
  return (res_SAFE)
  #dir = paste0(dir, "SAFE_3.csv")
  
  #write.csv(res_SAFE, file = dir, row.names = FALSE)
  
}

all_alg <- c("seurat", "bigScale", "cellranger", "pcaReduce", "ascend", "cidr", "countClust", 
             "monocle", "sc3", "rca", "tscan", "simlr", "sincera", "raceid", "scran", 
             "phenograph", "altanalyze")

all.datasets <- c("JA_melanoma_allCells", "tirosh_melanoma_allCells",
                   "li_crc_allCells")

setwd("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/Clustering_Results")

for (dataset in all.datasets) {
  print(dataset)
  clustering_results <- read.csv(paste0(dataset, ".csv"))
  alg <- intersect(all_alg, colnames(clustering_results))
  res_dataset <- data.frame(cell = clustering_results$cell, truth = clustering_results$truth)
  for (i in 3:length(alg)) {
    res <- run_SAFE(clustering_results, alg[1:i])
    colnames(res) <- c("cell", paste0("SAFE_", i))
    res_dataset <- cbind(res_dataset, res[ncol(res)])
  }
  out_dir <- paste0("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/ensemble_results/",
                    dataset, "_SAFE.csv")
  write.csv(file=out_dir, res_dataset, row.names=FALSE)
}


run_SAME <- function(input_dataframe, tools) {
  
  print(tools)
  final <- input_dataframe
  res_SAME <- data.frame(cell = final$cell)
  final <- final[,tools]
  final$cell=NULL
  
  final <- t(final)
  dims <- dim(final)
  safe_input <- t(apply(final, 1, as.numeric))
  
  same_input <- t(safe_input)
  rownames(same_input) = NULL
  colnames(same_input) = NULL
  
  cluster.SAME <- SAMEclustering(Y = same_input, rep = 3, SEED = 123)
  
  res_SAME <- data.frame(cell=res_SAME$cell, SAME_AIC = cluster.SAME$AICcluster ,SAME_BIC = cluster.SAME$BICcluster)
  
  return (res_SAME)
  
  #filename1=paste(dir, "SAME_AIC_3.csv")
  #filename2=paste(dir, "SAME_BIC_3.csv")
  
  #write.csv(data.frame(cell=res_SAME$cell, SAME_AIC_3 = res_SAME$SAME_AIC), file = filename1, row.names = FALSE)
  #write.csv(data.frame(cell = res_SAME$cell, SAME_BIC_3 = res_SAME$SAME_BIC), file = filename2, row.names = FALSE)
  
}

setwd("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/Clustering_Results")
all_alg <- c("seurat", "bigScale", "cellranger", "pcaReduce", "ascend", "cidr", "countClust", 
             "monocle", "sc3", "rca", "tscan", "simlr", "sincera", "raceid", "scran", 
             "phenograph", "altanalyze")

all.datasets <- c( "chung_breast_allCells",
                  "darmanis_glioblastoma_nonTumour", "JA_melanoma_nonTumour", "li_crc_nonTumour", 
                  "tirosh_melanoma_nonTumour", "chung_breast_nonTumour",
                  "darmanis_glioblastoma_allCells", "JA_melanoma_allCells", "tirosh_melanoma_allCells",
                  "li_crc_allCells")

# for (dataset in all.datasets) {
#   print (dataset)
#   clustering_results <- read.csv(paste0(dataset, ".csv"))
#   alg <- intersect(all_alg, colnames(clustering_results))
#   res_dataset_AIC <- data.frame(clustering_results$cell, clustering_results$truth)
#   res_dataset_BIC <- data.frame(clustering_results$cell, clustering_results$truth) d
#   for (col in colnames(clustering_results)) {
#     if (col != "cell") {
#       clustering_results <- transform(clustering_results,col=as.numeric(factor(clustering_results[,col])))
#       clustering_results[,col] <- clustering_results[,"col"]
#       clustering_results[,"col"] = NULL
#       clustering_results[,col][clustering_results[,col] == 0] <- max(unique(clustering_results[,col])) + 1
#     }
#   }
#   for (i in 3:length(alg)) {
#     print("SAME")
#     res <- run_SAME(clustering_results, alg[1:i])
#     colnames(res) <- c("cell", paste0("SAME_AIC_", i), paste0("SAME_BIC_", i))
#     res_dataset_AIC <- data.frame(cell = clustering_results$cell, truth = clustering_results$truth)
#     res_dataset_BIC <- data.frame(cell = clustering_results$cell, truth = clustering_results$truth)
#   }
#   out_dir_AIC <- paste0("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/ensemble_results/",
#                     dataset, "_SAME_AIC.csv")
#   out_dir_BIC <- paste0("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/ensemble_results/",
#                         dataset, "_SAME_BIC.csv")
#   write.csv(file=out_dir_AIC, res_dataset_AIC, row.names=FALSE)
#   write.csv(file=out_dir_BIC, res_dataset_BIC, row.names=FALSE)
# }



