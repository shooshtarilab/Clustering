source("Functions.R")
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(gplots)

all.algorithms <- c("altAnalyze",	"ascend", "bigScale",
                    "cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",
                    "phenograph",	"raceid", "rca",	"sc3",	"scran",	"seurat",
                    "sincera",	"tscan", "truth")

all.datasets <- c("peng_pancreatic_malignant",
                  "vangalen_AML_malignant", 
                  "lambrechts_lung_malignant",
                  "tirosh_melanoma_malignant",
                  "JA_melanoma_malignant")

setwd("/Users/alainamahalanabis/Documents/clustering/analysis_scripts")
num=0
total = 0
df <- c()
lowest_percent <- matrix(0, nrow=length(all.datasets), ncol=length(all.algorithms), dimnames = list(all.datasets, all.algorithms))
res <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("method", "f1_score", "truth", "dataset"))))
for (dataset in all.datasets) {
  print(dataset)
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  num_categories <- length(unique(data[,"truth"]))
  rare_cells = {}
  for (cell in unique(data[,"truth"])) {
    num_cells = nrow(data[data$truth==cell,])
    pct = num_cells / nrow(data)
    print (cell)
    print (pct)
    print (num_cells)
    if (num_cells < 500 && pct < 0.2 || pct < 0.05) {
      rare_cells <- c(rare_cells, cell)
    }
    
  }
  category_F_measure <- matrix(0, nrow=length(intersect(all.algorithms, colnames(data))), ncol=num_categories, dimnames=list(intersect(all.algorithms, colnames(data)), unique(data[,"truth"])))
  for (alg in intersect(all.algorithms, colnames(data))) {
    print (alg)
    data[,alg] <- as.numeric(as.factor(data[,alg]))
    num_clusters <- length(unique(data[,alg]))
    all_categories <- unique(data[,"truth"])
    category <- matrix(0, nrow=num_clusters, ncol=num_categories, dimnames=list(sort(unique(data[,alg])), unique(data[,"truth"])))
    F_measure <-F.Measure.Func(data[,"truth"], data[,alg])[2]
    F_measure <- cbind(c(1:num_clusters), F_measure[[1]])
    for (i in 1:num_clusters) {
      total_cells <- length(which(data[,alg] == i))
      for (j in all_categories) {
        cells_in_category <- length(intersect(which(data[,alg] == i), which(data[,"truth"] == j)))
        category[i, j] = cells_in_category/total_cells
      }
      F_measure[i,1] <- colnames(category)[which(category[i,]==max(category[i,]))[1]]
    }
    lowest_percent[dataset, alg] <- min(apply(category, 1, max))
    total = total + 1
    if (lowest_percent[dataset, alg] < 0.5) {
      num = num + 1
    }
    for (j in rare_cells) {
      category_F_measure[alg, j] <- median(as.numeric(F_measure[F_measure[,1] == j,2]))
      if (!is.na(category_F_measure[alg, j])) {
        res <- rbind(res, data.frame(method=alg, F1_score=category_F_measure[alg, j], category = j, dataset=dataset))
      } else {
        res <- rbind(res, data.frame(method=alg, F1_score=0, category = j, dataset=dataset))
      }
    }
  }
}

print(num)
print(total)
write.csv(res, file="bigdf.csv", sep=",", row.names=FALSE)

#CREATE A CATEGORY HEATMAP FOR EACH DATASET

bigdf <- read.csv("bigdf.csv", sep=',', stringsAsFactors = F)

categories = unique(bigdf$category)
all_methods = unique(bigdf$method)
datasets = unique(bigdf$dataset)

median_hm <- matrix(nrow = length(datasets), ncol=length(categories), dimnames = list(datasets, categories))

#dataset <- "peng_pancreatic_nonTumour"

for (dataset in datasets) {
  print(dataset)
  methods <- intersect(all_methods, bigdf[bigdf$dataset == dataset, "method"])
  hm <- matrix(nrow = length(categories), ncol=length(methods), dimnames = list(categories, methods))
  for (method in methods) {
    print (method)
    for (category in categories) {
      print(category)
      f_score <- bigdf[bigdf$category == category & bigdf$method == method & 
                         bigdf$dataset == dataset, "F1_score"]
      if (length(f_score) > 0) {
        hm[category, method] <- f_score
      } else {
        hm[category, method] <- NA  
      }
    }
  }
  hm <- hm[apply(hm,1,sum, na.rm = TRUE) > 0,]
  hm <- hm[order(apply(hm,1, median, na.rm = TRUE), decreasing=T),]
  median <- apply(hm,1, median, na.rm = TRUE)
  median <- median[order(median, decreasing=T)]
  
  hm <- cbind(median, hm)
  pdf(paste0(dataset, ".pdf"), width=15, height=15)
  pheatmap(hm, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,
           color = brewer.pal(n=9, name="Greens")[3:9], number_color="white",
           angle_col = 45, number_format = "%.2f", gaps_col = 1, cellheight = 50, main=dataset)
  
  dev.off()
  
  for (category in categories) {
    print(category)
    if (category %in% rownames(hm)) {
      median_hm[dataset, category] <- hm[category, "median"]
    } else {
      median_hm[dataset, category] <- NA
    }
  }
}


#COMBINED HEATMAP FOR ALL DATASETS

median_hm <- t(median_hm)
median_hm<- median_hm[,order(colnames(median_hm))]

colnames(median_hm) <- c("Breast : all cells", "Glioblastoma : all cells",
                         "Melanoma : all cells", 
                         "Lung : all cells", "Colorectal : all cells", 
                         "Pancreatic : all cells", 
                         "Metastatic melanoma : all cells",
                         "AML : all cells")

median_hm <- median_hm[order(apply(median_hm,1,median, na.rm = TRUE), decreasing=T),]

median <- apply(median_hm,1, median, na.rm = TRUE)

#median_hm <- median_hm[,order(apply(median_hm,2,sum, na.rm = TRUE), decreasing=T)]

median_hm <- cbind(median, median_hm)
mat2 <- median_hm
mat2 <- round(mat2, digits=2)
mat2[is.na(mat2)] <- ""


pdf("rare_cells.pdf", width=20, height=40)
pheatmap(median_hm, display_numbers = mat2, cluster_rows = FALSE, cluster_cols = FALSE,
         color = brewer.pal(n=9, name="Blues")[3:9], number_color="white",
         angle_col = 45, number_format = "%.2f", gaps_col = 1, cellheight = 40, main="All datasets and categories",
         fontsize_number = 12, cellwidth = 40)
dev.off()


#COMBINED RARE CELL TYPE

res <- read.csv("bigdf.csv", sep=',', stringsAsFactors = F)
res$dataset[res$dataset == "darmanis_glioblastoma_allCells"] <- "Glioblastoma : all cells"
res$dataset[res$dataset == "lambrechts_lung_allCells"] <- "Lung : all cells"
res$dataset[res$dataset == "JA_melanoma_allCells"] <- "Melanoma : all cells"
res$dataset[res$dataset == "li_crc_allCells"] <- "Colorectal : all cells"
res$dataset[res$dataset == "tirosh_melanoma_allCells"] <- "Metastatic melanoma : all cells"
res$dataset[res$dataset == "chung_breast_allCells"] <- "Breast : all cells"
res$dataset[res$dataset == "peng_pancreatic_allCells"] <- "Pancreatic : all cells"
res$dataset[res$dataset == "vangalen_AML_allCells"] <- "AML : all cells"


res$rare <- paste0(res$dataset, "_", res$category)

library("reshape2")
test <- dcast(res, method~rare, value.var="F1_score")
rownames(test) <- test$method
test$method=NULL
test <- as.matrix(test)
test <- test[,!apply(is.na(test), 2, any)]
median <- apply(test,1, median, na.rm = TRUE)
test <- cbind(median, test)
test <- test[order(-median),]
test <- test[-13,]

mat2 <- test
mat2 <- round(mat2, digits=2)
mat2[is.na(mat2)] <- ""



pdf("rare_cells.pdf", width=32, height=15)
pheatmap(test, display_numbers = mat2, cluster_rows = FALSE, cluster_cols = FALSE,
         color = brewer.pal(n=9, name="Blues")[3:9], number_color="white",
         angle_col = 45, number_format = "%.2f", gaps_col = 1, cellheight = 40, main="All datasets and rare cell types",
         fontsize_number = 12, cellwidth = 40)

dev.off()
