rm(list = ls())
Work.Dir <- "/Users/pVIsa shooshtVI/Documents/Projects/scRNAseq-Clustering"
source("Functions.R")

all.datasets <- c("darmanis_glioblastoma_allCells", "JA_melanoma_allCells", "li_crc_allCells", "tirosh_melanoma_allCells", "chung_breast_allCells")

# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
all.algorithms <- c("altAnalyze",	"ascend",	"cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",	"phenograph",	"rca",	"sc3",	"scran",	"seurat",	"sincera",	"tscan")


low.CI <- high.CI <- VI <- cov.data <- cor.data <- VI.avg <- list()
for (dataset in all.datasets) {
  VI[[dataset]] <- cor.data[[dataset]] <- cov.data[[dataset]] <- VI.avg[[dataset]] <- matrix(0, nrow=length(all.algorithms), ncol=length(all.algorithms), dimnames=list(all.algorithms, all.algorithms))
  
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
      VI[[dataset]][alg1, alg2] <- VI.Func(C, K)
    }
  }
  
  # Measuring average F-Measure
  for (alg1 in all.algorithms) {
    for (alg2 in all.algorithms) {
      VI.avg[[dataset]][alg1, alg2] <- round((VI[[dataset]][alg1, alg2] + VI[[dataset]][alg2, alg1])/2, 2)
    }
  }
}


VI.total <- (VI.avg[[1]][all.algorithms, all.algorithms] +
                VI.avg[[2]][all.algorithms, all.algorithms] +
                VI.avg[[3]][all.algorithms, all.algorithms] +
                VI.avg[[4]][all.algorithms, all.algorithms] +
                VI.avg[[5]][all.algorithms, all.algorithms])/5


pdf(file=paste0("./Figures/VI-Similar-Algorithms-Heatmap.pdf"), width=8, height=8)
heatmap(VI.total, keep.dendro=T, scale="none", margins = c(7, 7))
dev.off()


f.data <- 1-VI.total
hc <- hclust(d=dist(f.data), method = 'complete')

pdf(file=paste0("./Figures/VI-Similar-Algorithms-Hierarchical.pdf"), width=8, height=8)
plot(hc)
dev.off()
memb <- cutree(hc, k = 8)
