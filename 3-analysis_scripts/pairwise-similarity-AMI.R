rm(list = ls())
Work.Dir <- "/Users/pAMIsa shooshtAMI/Documents/Projects/scRNAseq-Clustering"
source("Functions.R")


all.datasets <- c("lambrechts_lung_allCells", "lambrechts_lung_nonTumour",
                  "darmanis_glioblastoma_allCells", "JA_melanoma_allCells", "li_crc_allCells",
                  "tirosh_melanoma_allCells", "chung_breast_allCells",
                  "darmanis_glioblastoma_nonTumour", "JA_melanoma_nonTumour", "li_crc_nonTumour",
                  "tirosh_melanoma_nonTumour", "chung_breast_nonTumour",
                  "peng_pancreatic_allCells", "peng_pancreatic_nonTumour",
                  "vangalen_AML_allCells", "vangalen_AML_nonTumour")

# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
all.algorithms <- c("altAnalyze",	"ascend", "bigScale",
                    "cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",
                    "phenograph",	"raceid", "rca",	"sc3",	"scran",	"seurat",
                    "simlr", "sincera",	"tscan")



low.CI <- high.CI <- AMI <- cov.data <- cor.data <- AMI.avg <- list()
for (dataset in all.datasets) {
  AMI[[dataset]] <- cor.data[[dataset]] <- cov.data[[dataset]] <- AMI.avg[[dataset]] <- matrix(0, nrow=length(all.algorithms), ncol=length(all.algorithms), dimnames=list(all.algorithms, all.algorithms))
  
  print(dataset)
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  
  for (alg1 in intersect(all.algorithms, colnames(data))) {
    for (alg2 in intersect(all.algorithms, colnames(data))) {
      print (alg1)
      print (alg2)
      C <- data[, alg1]
      K <- data[, alg2]
      
      # Identify Mean AMI
      AMI[[dataset]][alg1, alg2] <- AMI.Func(C, K)
    }
  }
  
  # Measuring average AMI
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      AMI.avg[[dataset]][alg1, alg2] <- round((AMI[[dataset]][alg1, alg2] + AMI[[dataset]][alg2, alg1])/2, 2)
    }
  }
}


AMI.total <- (AMI.avg[[1]][all.algorithms, all.algorithms] +
                AMI.avg[[2]][all.algorithms, all.algorithms] +
                AMI.avg[[3]][all.algorithms, all.algorithms] +
                AMI.avg[[4]][all.algorithms, all.algorithms] +
                AMI.avg[[5]][all.algorithms, all.algorithms] +
                AMI.avg[[6]][all.algorithms, all.algorithms] +
                AMI.avg[[7]][all.algorithms, all.algorithms] +
                AMI.avg[[8]][all.algorithms, all.algorithms] +
                AMI.avg[[9]][all.algorithms, all.algorithms] +
                AMI.avg[[10]][all.algorithms, all.algorithms] +
                AMI.avg[[11]][all.algorithms, all.algorithms] +
                AMI.avg[[12]][all.algorithms, all.algorithms] +
                AMI.avg[[13]][all.algorithms, all.algorithms] +
                AMI.avg[[14]][all.algorithms, all.algorithms] +
                AMI.avg[[15]][all.algorithms, all.algorithms] +
                AMI.avg[[16]][all.algorithms, all.algorithms])

num_entries <- matrix(0, nrow=length(all.algorithms), ncol=length(all.algorithms), dimnames=list(all.algorithms, all.algorithms))
for (dataset in all.datasets) {
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      if (AMI.avg[[dataset]][alg1, alg2] != 0) {
        num_entries[alg1, alg2] = num_entries[alg1, alg2] + 1
      }
    }
  }
}

for (alg1 in all.algorithms) {
  for (alg2 in all.algorithms) {
    print(alg1)
    print(alg2)
    AMI.total[alg1, alg2] <- AMI.total[alg1, alg2] / num_entries[alg1, alg2]
  }
}

write.csv(AMI.total, file = "./pairwise_results/AMI.csv")


pdf(file=paste0("./Figures/AMI-Similar-Algorithms-Heatmap.pdf"), width=8, height=8)
heatmap(AMI.total, keep.dendro=T, scale="none", margins = c(7, 7))
dev.off()


f.data <- 1-AMI.total
hc <- hclust(d=dist(f.data), method = 'complete')

pdf(file=paste0("./Figures/AMI-Similar-Algorithms-Hierarchical.pdf"), width=8, height=8)
plot(hc)
dev.off()
memb <- cutree(hc, k = 8)
