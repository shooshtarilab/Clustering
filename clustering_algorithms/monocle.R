#devtools::install_github('cole-trapnell-lab/leidenbase')
#BiocManager::install('batchelor')
#devtools::install_github('cole-trapnell-lab/monocle3')
library("monocle3")
library("igraph")

main <- function(filename) {

  load(filename)
  
  #CREATE A CELLDATASET OBJECT
  
  # pheno.data <- colnames(counts(sce)) # Get sample name
  # #pheno.data <- unlist(lapply(pheno.data, function(x) strsplit(x, '_')[[1]][1])) # Get cell type from sample name
  # pheno.data.df <- data.frame(type=pheno.data) # Must be data frame object
  # rownames(pheno.data.df) <- colnames(counts(sce)) # Rownames must match expression data
  # feature.data.df <- data.frame(type=rownames(counts(sce)))
  # rownames(feature.data.df) <- rownames(counts(sce))
  # pd <- new('AnnotatedDataFrame', data = pheno.data.df)
  # fd <- new("AnnotatedDataFrame", data = feature.data.df)
  cds <- new_cell_data_set(counts(sce))
  
  #cds <- estimateSizeFactors(cds)
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  print(dataset)
  folder = paste0("./", dataset, "/monocle.csv")
  start = Sys.time()
  
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 50)
  
  cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "PCA")
  
  cds <- cluster_cells(cds, cluster_method = 'louvain', reduction_method = "tSNE", k=20, verbose=T)
  
  res <- data.frame(cell = row.names(colData(sce)), monocle = cds@clusters$tSNE$clusters)
  
  end = Sys.time()
  time = end - start
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " monocle"), file="time.txt", append=TRUE, sep = "\n")
}

datasets <- c("vangalen_AML_nonTumour_counts.RData")

for (i in datasets) {
  main(i)
}


