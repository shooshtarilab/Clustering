BiocManager::install("SingleCellExperiment")
library(Seurat)
library(SingleCellExperiment)

main <- function(filename) {
  
   print(filename)
  
   dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
   folder = paste0("./", dataset, "/seurat.csv")
     
   load(filename)
   
   start = Sys.time()
   mat <-  CreateSeuratObject(counts = counts(sce),
                                min.cells = 3,  project = "test", min.features = 0,
                                names.field = 1, names.delim = "_", meta.data = NULL)
    
    
    mat <- NormalizeData(object = mat, normalization.method = "LogNormalize", scale.factor = 10000,
                          assay.type = "RNA", display.progress = FALSE)
    
    mat <- FindVariableFeatures(object = mat, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

    mat <- ScaleData(object = mat, 
                      model.use = "linear", use.umi = FALSE, do.scale = TRUE,
                      do.center = TRUE, scale.max = 10, block.size = 1000,
                      min.cells.to.block = 3000, display.progress = FALSE, assay.type = "RNA",
                      do.cpp = TRUE, check.for.norm = TRUE, do.par = FALSE, num.cores = 1)
    
    set.seed(635465465)
    mat <- RunPCA(object = mat, pc.genes = mat@var.genes, do.print = FALSE, reduction.name = "pca",
                   reduction.key = "PC", assay.type = "RNA", seed.use = 42)
    mat <- JackStraw(mat, dims = 20, num.replicate = 100, prop.freq = 0.01,
                      maxit = 1000)
  
    set.seed(44487152)
    

    mat <- FindNeighbors(object = mat, reduction = "pca", dims = 1:10, k.params=20, prune.SNN = 1/15)

    mat <- FindClusters(object = mat,
                         modularity.fxn = 1, resolution = 0.6, algorithm = 1, n.start = 100,
                         n.iter = 10, random.seed = 0)
    
    colData(sce)$Seurat <- mat@meta.data$res.0.8
    res <- colData(sce)
    
    res <- mat@meta.data$RNA_snn_res.0.6
    final <- data.frame(colnames(mat), as.numeric(res))
    colnames(final) <- c("cell", "seurat")
    
    end = Sys.time()
    
    time = end-start

    write.csv(final, file=folder, sep=",", quote=FALSE, row.names = TRUE, col.names=TRUE)
    cat(paste0(time, " ", dataset, " seurat"), file="time.txt", append=TRUE, sep = "\n")
  
}

