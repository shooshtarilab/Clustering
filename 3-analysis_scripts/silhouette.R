suppressPackageStartupMessages(library(SingleCellExperiment))

pre_clean <- function(sce) {

  ave.counts <- rowMeans(as.matrix(counts(sce)))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]
  #del_gns <- c(which(rowData(sce1)$is_feature_control_rbp), which(rowData(sce1)$is_feature_control_mt))
  #sce1 <- sce1[-del_gns,]  

  return(sce1)
}

all.datasets <- c("lambrechts_lung_allCells", "lambrechts_lung_nonTumour",
                  "darmanis_glioblastoma_allCells", "JA_melanoma_allCells", "li_crc_allCells", 
                  "tirosh_melanoma_allCells", "chung_breast_allCells",
                  "darmanis_glioblastoma_nonTumour", "JA_melanoma_nonTumour", "li_crc_nonTumour", 
                  "tirosh_melanoma_nonTumour", "chung_breast_nonTumour",
                  "peng_pancreatic_allCells", "peng_pancreatic_nonTumour",
                  "vangalen_AML_allCells", "vangalen_AML_nonTumour")


# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
all.algorithms <- c("SAME_AIC_3", "SAME_AIC_5", "SAME_BIC_3", "SAME_BIC_5",
                    "clue_3", "clue_5", "SAFE_3", "SAFE_5",
                    "altAnalyze",       "ascend", "bigScale",
                    "cellranger",       "cidr", "countClust",   "monocle",      "pcaReduce",
                    "phenograph",       "raceid", "rca",        "sc3",  "scran",        "seurat",
                    "simlr", "sincera", "tscan")

silhouette_matrix <- matrix(0, nrow=length(all.datasets), ncol=length(all.algorithms), dimnames=list(all.datasets, all.algorithms))

for (dataset in all.datasets) {
  
  filename = paste0(dataset, "_TPM.RData")
  print(filename)
  load(filename)
  data.path <- paste0("Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T, stringsAsFactors = FALSE)
  sce1 <- pre_clean(sce)  
  mat <- counts(sce1)[,colnames(counts(sce1)) %in% data[,"cell"]]
  rownames(mat)=NULL
  data <- data[match(colnames(mat), data$cell),]
  mat <- log2(mat + 1)
  mat <- t(mat)
  mat <- mat[ , which(apply(mat, 2, var) != 0)]
  matrix.pca <- prcomp(mat, center = TRUE,scale. = TRUE) 
  D = dist(matrix.pca$x[,1:10])
  
  for (alg in intersect(all.algorithms, colnames(data))) {
    print (alg)
    sil <- cluster::silhouette(data[,alg], D)
    silhouette_matrix[dataset, alg] <- mean(sil[,3])
    write.csv(silhouette_matrix, file="dataframe_results/silhouette.csv")
  }
}

write.csv(silhouette_matrix, file="silhouette.csv")
