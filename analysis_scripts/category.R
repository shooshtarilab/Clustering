# FOR ALL DATASETS AND ALL ALGORITHMS
# CALCULATE F SCORE FOR EACH CLUSTER   
# ASSIGN EACH CLUSTER A CATEGORY BASED ON MAJORITY CELL TYPE
# TAKE MEDIAN F SCORE ACCROSS ALL CLUSTERS FOR A CATEGORY TO GET F SCORE FOR THAT CATEGORY
# FINAL DATAFRAME STORED IN 'bigdf'

source("Functions.R")

all.datasets <- c("peng_pancreatic_allCells", "peng_pancreatic_nonTumour",
                  "vangalen_AML_allCells", "vangalen_AML_nonTumour",
                  "lambrechts_lung_allCells", "lambrechts_lung_nonTumour",
                  "darmanis_glioblastoma_allCells", "JA_melanoma_allCells", "li_crc_allCells",
                  "tirosh_melanoma_allCells", "chung_breast_allCells",
                  "darmanis_glioblastoma_nonTumour", "JA_melanoma_nonTumour", "li_crc_nonTumour",
                  "tirosh_melanoma_nonTumour", "chung_breast_nonTumour")

all.algorithms <- c("SAME_AIC_3", "SAME_AIC_5", "SAME_BIC_3", "SAME_BIC_5",
                    "clue_3", "clue_5", "SAFE_3", "SAFE_5",
                    "altAnalyze",	"ascend", "bigScale",
                    "cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",
                    "phenograph",	"raceid", "rca",	"sc3",	"scran",	"seurat",
                    "simlr", "sincera",	"tscan")


num =0
total = 0
df <- c()
lowest_percent <- matrix(0, nrow=length(all.datasets), ncol=length(all.algorithms), dimnames = list(all.datasets, all.algorithms))
res <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("method", "f1_score", "category", "dataset"))))
for (dataset in all.datasets) {
  print(dataset)
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  num_categories <- length(unique(data[,"category"]))
  category_F_measure <- matrix(0, nrow=length(intersect(all.algorithms, colnames(data))), ncol=num_categories, dimnames=list(intersect(all.algorithms, colnames(data)), unique(data[,"category"])))
  for (alg in intersect(all.algorithms, colnames(data))) {
    print (alg)
    data[,alg] <- as.numeric(as.factor(data[,alg]))
    num_clusters <- length(unique(data[,alg]))
    all_categories <- unique(data[,"category"])
    category <- matrix(0, nrow=num_clusters, ncol=num_categories, dimnames=list(sort(unique(data[,alg])), unique(data[,"category"])))
    F_measure <-F.Measure.Func(data[,"truth"], data[,alg])[2]
    F_measure <- cbind(c(1:num_clusters), F_measure[[1]])
    for (i in 1:num_clusters) {
      total_cells <- length(which(data[,alg] == i))
      for (j in all_categories) {
        cells_in_category <- length(intersect(which(data[,alg] == i), which(data[,"category"] == j)))
        category[i, j] = cells_in_category/total_cells
      }
      F_measure[i,1] <- colnames(category)[which(category[i,]==max(category[i,]))[1]]
    }
    lowest_percent[dataset, alg] <- min(apply(category, 1, max))
    total = total + 1
    if (lowest_percent[dataset, alg] < 0.5) {
      num = num + 1
    }
    for (j in as.character(all_categories)) {
      category_F_measure[alg, j] <- median(as.numeric(F_measure[F_measure[,1] == j,2]))
      if (!is.na(category_F_measure[alg, j])) {
        res <- rbind(res, data.frame(method=alg, F1_score=category_F_measure[alg, j], category = j, dataset=dataset))
      }
    }
  }
}

print(num)
print(total)
write.csv(res, file="bigdf.csv", sep=",", row.names=FALSE)