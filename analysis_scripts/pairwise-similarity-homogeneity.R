rm(list = ls())
Work.Dir <- "/Users/pHomogeneitysa shooshtHomogeneity/Documents/Projects/scRNAseq-Clustering"
source("Functions.R")


##################################################################
#######                 Setting Parameters                ########
##################################################################
ref <- "truth"
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


low.CI <- high.CI <- Homogeneity <- cov.data <- cor.data <- Homogeneity.avg <- list()
for (dataset in all.datasets) {
  low.CI[[dataset]] <- high.CI[[dataset]] <- Homogeneity[[dataset]] <- cor.data[[dataset]] <- cov.data[[dataset]] <- Homogeneity.avg[[dataset]] <- matrix(0, nrow=length(all.algorithms), ncol=length(all.algorithms), dimnames=list(all.algorithms, all.algorithms))
  
  print(dataset)
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  
  for (alg1 in intersect(all.algorithms, colnames(data))) {
    for (alg2 in intersect(all.algorithms, colnames(data))) {
      print (alg1)
      print (alg2)
      C <- data[, alg1]
      K <- data[, alg2]
      
      # Identify Mean F measures
      Homogeneity[[dataset]][alg1, alg2] <- Homogeneity.Func(C, K)
    }
  }
  
  # Measuring average F-Measure
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      Homogeneity.avg[[dataset]][alg1, alg2] <- round((Homogeneity[[dataset]][alg1, alg2] + Homogeneity[[dataset]][alg2, alg1])/2, 2)
    }
  }
}

Homogeneity.total <- (Homogeneity.avg[[1]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[2]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[3]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[4]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[5]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[6]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[7]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[8]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[9]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[10]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[11]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[12]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[13]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[14]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[15]][all.algorithms, all.algorithms] +
                Homogeneity.avg[[16]][all.algorithms, all.algorithms])


num_entries <- matrix(0, nrow=length(all.algorithms), ncol=length(all.algorithms), dimnames=list(all.algorithms, all.algorithms))
for (dataset in all.datasets) {
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      if (Homogeneity.avg[[dataset]][alg1, alg2] != 0) {
        num_entries[alg1, alg2] = num_entries[alg1, alg2] + 1
      }
    }
  }
}

for (alg1 in all.algorithms) {
  for (alg2 in all.algorithms) {
    print(alg1)
    print(alg2)
    Homogeneity.total[alg1, alg2] <- Homogeneity.total[alg1, alg2] / num_entries[alg1, alg2]
  }
}

write.csv(Homogeneity.total, file = "./pairwise_results/Homogeneity.csv")



pdf(file=paste0("./Figures/Homogeneity-Similar-Algorithms-Heatmap.pdf"), width=8, height=8)
heatmap(Homogeneity.total, keep.dendro=T, scale="none", margins = c(7, 7))
dev.off()


f.data <- 1-Homogeneity.total
hc <- hclust(d=dist(f.data), method = 'complete')

pdf(file=paste0("./Figures/Homogeneity-Similar-Algorithms-Hierarchical.pdf"), width=8, height=8)
plot(hc)
dev.off()
memb <- cutree(hc, k = 8)
