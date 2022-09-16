#CREATE A CATEGORY HEATMAP FOR EACH DATASET

bigdf <- read.csv("bigdf.csv", sep=',', stringsAsFactors = F)

categories = unique(bigdf$category)
all_methods = unique(bigdf$method)
datasets = unique(bigdf$dataset)

median_hm <- matrix(nrow = length(datasets), ncol=length(categories), dimnames = list(datasets, categories))

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
  pdf(paste0(dataset, ".pdf"), width=15, height=4)
  pheatmap(hm, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,
           color = brewer.pal(n=9, name="Greens")[3:9], number_color="white",
           angle_col = 45, number_format = "%.2f", gaps_col = 1, cellheight = 50, main=dataset)
  
  dev.off()
  
  for (category in categories) {
    if (category %in% rownames(hm)) {
      median_hm[dataset, category] <- hm[category, "median"]
    } else {
      median_hm[dataset, category] <- NA
    }
  }
}


#COMBINED HEATMAP FOR ALL DATASETS

median_hm <- t(median_hm)
median_hm <- median_hm[order(apply(median_hm,1,median, na.rm = TRUE), decreasing=T),]
median <- apply(median_hm,1, median, na.rm = TRUE)

median_hm <- median_hm[,order(apply(median_hm,2,sum, na.rm = TRUE), decreasing=T)]

median_hm <- cbind(median, median_hm)

pdf("category_all_datasets.pdf", width=20, height=20)
pheatmap(median_hm, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,
         color = brewer.pal(n=9, name="Blues")[3:9], number_color="white",
         angle_col = 45, number_format = "%.2f", gaps_col = 1, cellheight = 40, main="All datasets and categories",
         fontsize_number = 12, cellwidth = 40)
dev.off()