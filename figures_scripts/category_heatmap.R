#CREATE A CATEGORY HEATMAP FOR EACH DATASET

bigdf <- read.csv("./category/bigdf.tsv", sep=',')

categories = unique(bigdf$category)
methods = unique(bigdf$method)
datasets = unique(bigdf$dataset)

mean_hm <- matrix(nrow = length(datasets), ncol=length(categories), dimnames = list(datasets, categories))

for (dataset in datasets) {
  hm <- matrix(nrow = length(categories), ncol=length(methods), dimnames = list(categories, methods))
  for (method in methods) {
    for (category in categories) {
      f_score <- bigdf[bigdf$category == category & bigdf$method == method & 
                                    bigdf$dataset == dataset, "f1_score"]
      if (length(f_score) > 0) {
        hm[category, method] <- f_score
      } else {
        hm[category, method] <- 0  
      }
    }
  }
  hm <- hm[rowSums(hm) > 0,]
  hm <- hm[order(rowMeans(hm), decreasing=T),]
  mean <- rowMeans(hm)
  mean <- mean[order(mean, decreasing=T)]
  
  hm <- cbind(mean, hm)
  pdf(paste0(dataset, ".pdf"), width=15, height=4)
  pheatmap(hm, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,
           color = brewer.pal(n=9, name="Greens")[3:9], number_color="white",
           angle_col = 45, number_format = "%.2f", gaps_col = 1, cellheight = 50, main=dataset)
  
  dev.off()
  
  for (category in categories) {
    if (category %in% rownames(hm)) {
      mean_hm[dataset, category] <- hm[category, "mean"]
    } else {
      mean_hm[dataset, category] <- 0
    }
  }
}


#COMBINED HEATMAP FOR ALL DATASETS

mean_hm <- mean_hm[,order(colMeans(mean_hm), decreasing=T)]
mean_hm <- mean_hm[,colSums(mean_hm) > 0]
mean <- rowMeans(mean_hm)

mean_hm <- cbind(mean, mean_hm)
pdf("category_all_datasets.pdf", width=20, height=20)
pheatmap(mean_hm, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,
         color = brewer.pal(n=9, name="Blues")[3:9], number_color="white",
         angle_col = 45, number_format = "%.2f", gaps_col = 1, cellheight = 40, main="All datasets and categories",
         fontsize_number = 12, cellwidth = 40)
dev.off()
