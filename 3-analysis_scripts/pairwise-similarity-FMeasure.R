rm(list = ls())
Work.Dir <- "/Users/parisa shooshtari/Documents/Projects/scRNAseq-Clustering"
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



low.CI <- high.CI <- f.measure <- cov.data <- cor.data <- f.measure.avg <- list()
for (dataset in all.datasets) {
  low.CI[[dataset]] <- high.CI[[dataset]] <- f.measure[[dataset]] <- cor.data[[dataset]] <- cov.data[[dataset]] <- f.measure.avg[[dataset]] <- matrix(0, nrow=length(all.algorithms), ncol=length(all.algorithms), dimnames=list(all.algorithms, all.algorithms))
  
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
      F.Measure.List <- F.Measure.Func(C, K)
      
      # Measure 95% confidence intervals using bootstrapping 
      F.vector <- F.Measure.List$F.vector
      CI <- quantile(F.Bootstrap.Func(F.vector), probs = c(0.025, 0.5, 0.975))
      low.CI[[dataset]][alg1, alg2] <- round(CI[1], 2)
      high.CI[[dataset]][alg1, alg2] <- round(CI[3], 2)
      f.measure[[dataset]][alg1, alg2] <- F.Measure.List$F.measure # round(CI[2], 2) #
      
      # cov.data[[dataset]][alg1, alg2] <- cov(C,K)
      # cor.data[[dataset]][alg1, alg2] <- cor(C,K, method = "spearman")
    }
  }
  
  # Measuring average F-Measure
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      f.measure.avg[[dataset]][alg1, alg2] <- round((f.measure[[dataset]][alg1, alg2] + f.measure[[dataset]][alg2, alg1])/2, 2)
    }
  }
}


f.measure.total <- (f.measure.avg[[1]][all.algorithms, all.algorithms] +
                      f.measure.avg[[2]][all.algorithms, all.algorithms] +
                      f.measure.avg[[3]][all.algorithms, all.algorithms] +
                      f.measure.avg[[4]][all.algorithms, all.algorithms] +
                      f.measure.avg[[5]][all.algorithms, all.algorithms])/5


pdf(file=paste0("./Figures/F-Measure-Similar-Algorithms-Heatmap.pdf"), width=8, height=8)
heatmap(f.measure.total, keep.dendro=T, scale="none", margins = c(7, 7))
dev.off()


f.data <- 1-f.measure.total
hc <- hclust(d=dist(f.data), method = 'complete')

pdf(file=paste0("./Figures/F-Measure-Similar-Algorithms-Hierarchical.pdf"), width=8, height=8)
plot(hc)
dev.off()
memb <- cutree(hc, k = 8)
