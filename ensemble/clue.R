#PACKAGES TO BE INSTALLED
install.packages("clue")
install.packages("rlist")
install.packages("relations")

#LOAD AND SOURCE THESE PACKAGES
library("clue")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/AAA.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/partition.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/ensemble.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/consensus.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/dissimilarity.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/agreement.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/utilities.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/membership.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/registration.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/objects.R")
library("relations")
library(rlist)


run_clue <- function(input_dataframe, tools) {
  final <- input_dataframe
  res <- data.frame(cell = final$cell)
  final <- final[tools]
  clue_input <- list()
  for (col in colnames(final)) {
    print (col)
    tmp <-list(cell = final$cell , class_ids=final[,col])
    class(tmp) <- "cluster"
    clue_input <- list.append(clue_input, tmp)
  }
  
  final$cell = NULL
  names(clue_input) <- colnames(final)
  
  ens<-cl_ensemble(list = clue_input)
  
  consensus <- cl_consensus(ens, "HE") 
  
  res <- cbind(res, cluster = c(cl_class_ids(consensus)))
  colnames(res)[ncol(res)] = "clue"
  #colnames(res) <- c("cell", "clue")
  return (res)
  #write.csv(res, file="clue.csv", row.names=FALSE)
  
}

setwd("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/Clustering_Results")
all_alg <- c("seurat", "bigScale", "cellranger", "pcaReduce", "ascend", "cidr", "countClust", 
         "monocle", "sc3", "rca", "tscan", "simlr", "sincera", "raceid", "scran", 
         "phenograph", "altanalyze")

all.datasets <- c("darmanis_glioblastoma_allCells", "JA_melanoma_allCells", "li_crc_allCells", 
                  "tirosh_melanoma_allCells", "chung_breast_allCells",
                  "darmanis_glioblastoma_nonTumour", "JA_melanoma_nonTumour", "li_crc_nonTumour", 
                  "tirosh_melanoma_nonTumour", "chung_breast_nonTumour")

for (dataset in all.datasets) {
  clustering_results <- read.csv(paste0(dataset, ".csv"))
  alg <- intersect(all_alg, colnames(clustering_results))
  res_dataset <- data.frame(cell = clustering_results$cell, truth = clustering_results$truth)
  for (i in 3:length(alg)) {
    res <- run_clue(clustering_results, alg[1:i])
    colnames(res) <- c("cell", paste0("clue_", i))
    res_dataset <- cbind(res_dataset, res[ncol(res)])
  }
  out_dir <- paste0("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/ensemble_results/",
                    dataset, "_clue.csv")
  write.csv(file=out_dir, res_dataset, row.names=FALSE)
}
