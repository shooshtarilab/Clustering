source("Functions.R")
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(gplots)

setwd("/Users/alainamahalanabis/Documents/clustering/analysis_scripts")

all.datasets <- c("peng_pancreatic_malignant",
                  "vangalen_AML_malignant", 
                  "lambrechts_lung_malignant",
                  "darmanis_glioblastoma_malignant", 
                  "JA_melanoma_malignant", "li_crc_malignant",
                  "tirosh_melanoma_malignant", "chung_breast_malignant")

all.algorithms <- c("altAnalyze",	"ascend", "bigScale",
                    "cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",
                    "phenograph",	"raceid", "rca",	"sc3",	"scran",	"seurat",
                    "sincera",	"tscan", "truth")

all.datasets <- c("peng_pancreatic_malignant",
                  "vangalen_AML_malignant", 
                  "lambrechts_lung_malignant",
                  "tirosh_melanoma_malignant",
                  "JA_melanoma_malignant")

#NUMBER OF CLUSTERS - ABSOLUTE VALUES

cluster_heatmap <- matrix(0, nrow=length(all.datasets), ncol=length(all.algorithms), dimnames = list(all.datasets, all.algorithms))
truth_heatmap <- matrix(0, nrow=length(all.datasets), ncol=2, dimnames = list(all.datasets))

for (dataset in all.datasets) {
  print (dataset)
  truth_heatmap[dataset,1] <- length(unique(data[,"truth"]))
  
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  if ("patient" %in% colnames(data)) {
    truth_heatmap[dataset,2] <- length(unique(data[,"patient"]))
  } else {
    truth_heatmap[dataset,2] <- NA
  }
  
  for (alg in intersect(all.algorithms, colnames(data))) {
    num_clusters <- length(unique(data[,alg]))
    cluster_heatmap[dataset, alg] <- num_clusters
  }
}

cluster_heatmap <- cbind(truth_heatmap, cluster_heatmap)
colnames(cluster_heatmap)[1] <- "aaaaaaa"
colnames(cluster_heatmap)[2] <- "aaaaaab"
cluster_heatmap <- cluster_heatmap[order(rownames(cluster_heatmap)) ,order(colnames(cluster_heatmap))]

'''
rownames(cluster_heatmap) <- c("Breast : all cells", "Breast : non-tumour", "Glioblastoma : all cells",
"Glioblastoma : non-tumour", "Melanoma : all cells", "Melanoma : non-tumour",
"Lung : all cells", "Lung : non-tumour", "Colorectal : all cells", 
"Colorectal : non-tumour", "Pancreatic : all cells", "Pancreatic : non-tumour",
"Metastatic melanoma : all cells", "Metastatic melanoma : non-tumour", 
"AML : all cells", "AML : non-tumour")
'''

rownames(cluster_heatmap) <- c("Breast : malignant", "Glioblastoma : malignant",
                               "Melanoma : malignant",
                               "Lung : malignant", "Colorectal : malignant", 
                               "Pancreatic : malignant",
                               "Metastatic melanoma : malignant",
                               "AML : malignant")
'''
rownames(cluster_heatmap) <- c("Melanoma : malignant",
"Lung : malignant",
"Pancreatic : malignant",
"Metastatic melanoma : malignant",
"AML : malignant")
'''

cluster_heatmap <- cluster_heatmap[order(rownames(cluster_heatmap)) ,order(colnames(cluster_heatmap))]

colnames(cluster_heatmap) <- c("truth", "patient", "AltAnalyze", "Ascend", "bigSCale", "Cell Ranger",
                               "CIDR", "CountClust", "Monocle", "pcaReduce", "PhenoGraph",
                               "RaceID", "RCA", "SC3", "Scran", "Seurat", "SIMLR", "SINCERA",
                               "TSCAN","truth")


cluster_heatmap <- t(cluster_heatmap)

cluster_heatmap[cluster_heatmap == 0] <- NA
mat2 <- cluster_heatmap
mat2[is.na(mat2)] <- ""

col1 <- brewer.pal(n=9, name="Blues")[3:9]
col2 <- brewer.pal(n=9, name="Purples")[3:9]

pdf("num_clusters_malignant_patient.pdf", width=8, height=8)
pheatmap(cluster_heatmap, cluster_rows=FALSE, cluster_cols = FALSE, scale="none", main = "Number of clusters
         for each malignant dataset and alg", color=c(col1, col2), breaks = c(0, 3, 6, 9, 12,
                                                                              15, 18, 20, 100, 200, 300, 400, 500), display_numbers = mat2, number_color="white",
         gaps_row = 2)

dev.off()
#NUMBER OF CLUSTERS - DIFFERENCE VALUES

#BiocManager::install("ComplexHeatmap")

cluster_heatmap <- matrix(0, nrow=length(all.datasets), ncol=length(all.algorithms), dimnames = list(all.datasets, all.algorithms))
truth_heatmap <- matrix(0, nrow=length(all.datasets), dimnames = list(all.datasets))

for (dataset in all.datasets) {
  print (dataset)
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  truth_heatmap[dataset,1] <- length(unique(data[,"truth"]))
  for (alg in intersect(all.algorithms, colnames(data))) {
    num_clusters <- length(unique(data[,alg]))
    cluster_heatmap[dataset, alg] <- num_clusters - truth_heatmap[dataset,1]
  }
}

cluster_heatmap <- cbind(truth_heatmap, cluster_heatmap)
colnames(cluster_heatmap)[1] <- "aaaaaaa"
cluster_heatmap <- cluster_heatmap[order(rownames(cluster_heatmap)) ,order(colnames(cluster_heatmap))]

'''
rownames(cluster_heatmap) <- c("Breast : all cells", "Breast : non-tumour", "Glioblastoma : all cells",
"Glioblastoma : non-tumour", "Melanoma : all cells", "Melanoma : non-tumour",
"Lung : all cells", "Lung : non-tumour", "Colorectal : all cells", 
"Colorectal : non-tumour", "Pancreatic : all cells", "Pancreatic : non-tumour",
"Metastatic melanoma : all cells", "Metastatic melanoma : non-tumour", 
"AML : all cells", "AML : non-tumour")
'''

rownames(cluster_heatmap) <- c("Breast : malignant", "Glioblastoma : malignant",
                               "Melanoma : malignant",
                               "Lung : malignant", "Colorectal : malignant", 
                               "Pancreatic : malignant",
                               "Metastatic melanoma : malignant",
                               "AML : malignant")


cluster_heatmap <- cluster_heatmap[order(rownames(cluster_heatmap)) ,order(colnames(cluster_heatmap))]

colnames(cluster_heatmap) <- c("truth", "AltAnalyze", "Ascend", "bigSCale", "Cell Ranger",
                               "CIDR", "CountClust", "Monocle", "pcaReduce", "PhenoGraph",
                               "RaceID", "RCA", "SC3", "Scran", "Seurat", "SIMLR", "SINCERA",
                               "TSCAN")


cluster_heatmap <- t(cluster_heatmap)
cluster_heatmap[cluster_heatmap == 0] <- NA
mat2 <- cluster_heatmap
mat2[is.na(mat2)] <- ""
median <- rowMedians(cluster_heatmap, na.rm = TRUE)
cluster_heatmap <- cbind(median, cluster_heatmap)
#annotations = data.frame(truth = cluster_heatmap[1,])
#anno_colors = data.frame(truth = cluster_heatmap[1,])
#anno_colors$truth <- colorRampPalette(rev(brewer.pal(9, "Oranges")) )(17)
mat2 <- cluster_heatmap
mat2[is.na(mat2)] <- ""
#cluster_heatmap <- cluster_heatmap[-1,]
#mat2 <- cluster_heatmap
#mat2[is.na(mat2)] <- ""

col1 <- rev(brewer.pal(n=9, name="Blues")[3:9])
col2 <- brewer.pal(n=9, name="Greens")[3:9]

col1 <- c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C")
col2 <- c("#00441B", "#006D2C",  "#238B45", "#41AB5D", "#74C476", "#A1D99B")


pdf("diff_clusters.pdf", width=8, height=8)
pheatmap(cluster_heatmap, cluster_rows=FALSE, cluster_cols = FALSE, scale="none", main = "cluster - truth malignant", color=c(col1, col2), breaks = c(-30, -25, -20, -15, -10, -5,
                                                                                                                                                      0, 5, 10, 15, 20, 25, 30), display_numbers = mat2, number_color="white",
         gaps_col=1, gaps_row=1)

dev.off()

library(ggplot2)

ggplot(annotations, aes(x = rownames(annotations), y = 1, fill = truth)) + 
  geom_tile() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  scale_color_brewer(palette="Dark2")
