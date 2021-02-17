# install.packages("devtools")
# library(devtools)
# install_github("GIS-SP-Group/RCA")
# library(RCA)
# library(scater)
# BiocManager::install("WGCNA") 

suppressPackageStartupMessages(library(RCA))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(optparse))

pre_clean <- function(sce) {
  
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 1
  sce1 <- sce[keep,]

  return(sce1)
}

main <- function(filename) {
    
  # option_list <- list(
  #   make_option(c("-i", "--input"), default="NA",
  #               help="A RData file containing a SingleCellExperiment object"),
  #   #
  #   make_option(c("-o", "--outdir"), default="NA",
  #               help="A path/name for the results directory")
  # )
  # opt <- parse_args(OptionParser(option_list=option_list))
  # Input<- opt$input
  # Outdir<- opt$outdir
  # 
  
  load(filename)
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/rca.csv")
  
  sce1 <- pre_clean(sce)
  genes <- rowData(sce1)$genes
  genes <- gsub("'", '', genes)
  rowData(sce1) <- data.frame(mat.gene_name = genes)
  #View(rowData(sce1))
  
  start = Sys.time()
  obj <- counts(sce1)
  rownames(obj) <- rowData(sce1)$mat.gene_name
  data_obj <- dataConstruct(obj)
  data_obj <- geneFilt(obj_in = data_obj, method = "default")
  #data_obj <- cellNormalize(data_obj, method = "no_norm")
  
  #IF USING ALREADY NORMALIZED DATASET
  data_obj$fpkm_transformed <- data_obj$fpkm_raw
  data_obj <- dataTransform(data_obj, method = "log10")
  #data_obj <- dataTransform(data_obj, method = "log10")
    
  data_obj <- featureConstruct(data_obj, method="GlobalPanel")
  
  
  #data_obj <- featureConstruct(data_obj, method = "ColonEpitheliumPanel")
  set.seed(20742579)
  data_obj <- cellClust(data_obj, method = "hclust", deepSplit_wgcna = 1,
                        min_group_Size_wgcna = 5)
  
  end = Sys.time()
  
  res <- data.frame(row.names(colData(sce)), data_obj$group_labels_color$groupLabel)
  colnames(res) <- c("cell", "rca")
  
  time = end - start
  write.csv(res, file=folder, quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " rca"), file="time.txt", append=TRUE, sep = "\n")

}

