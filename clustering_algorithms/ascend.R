library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)
library(ascend)
library(SingleCellExperiment)
library("RSpectra")
library("fastcluster")
library(Matrix)
library(tidyverse)
library(dplyr)


main <- function(filename) {
  
    load(filename)
   dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
    folder = paste0("./", dataset, "/ascend.csv")
      
    expression_matrix <- data.frame(counts(sce))
    expression_matrix$gene_name <- rowData(sce)$genes
    expression_matrix <- expression_matrix %>% distinct(expression_matrix$gene_name, .keep_all = TRUE)
    rownames(expression_matrix) <- expression_matrix$gene_name
    expression_matrix$gene_name=NULL
    expression_matrix$`expression_matrix$gene_name`=NULL

    ave.counts <- rowMeans(expression_matrix)

    keep <- ave.counts >= 5 
    expression_matrix <- expression_matrix[keep,]

    row_info_data_frame <- data.frame(row.names = rownames(expression_matrix), mat.gene_name = rownames(expression_matrix))

    rownames(row_info_data_frame) <- row_info_data_frame$mat.gene_name
    rownames(expression_matrix) <- row_info_data_frame[,1]
  
    col_info_data_frame <- as.data.frame(colData(sce)[, c('cells'), drop=FALSE])
    colnames(col_info_data_frame) <- "cell_barcode"
    rownames(col_info_data_frame) <- col_info_data_frame$cell_barcode
    gene_list <- rownames(expression_matrix)
    
    control_list <- list(Mt = row_info_data_frame$mat.gene_name[grep("^MT", row_info_data_frame$mat.gene_name, ignore.case = TRUE)],
                     Rb = row_info_data_frame$mat.gene_name[grep("^Rps|^Rpl", row_info_data_frame$mat.gene_name, ignore.case = TRUE)])

    start = Sys.time()
    em.set <- EMSet(expression_matrix,
                    colInfo = col_info_data_frame,
                    rowInfo = row_info_data_frame,
                    controls = control_list)

    ## Check gene names for EMSet match those in rowInfo
    # if (!(identical(gene_list, rownames(rowInfo)))){
    #   errors <- c(errors, "rowInfo rownames do not match EMSet rownames.")
    # } else{
    #   if (!(identical(as.vector(rowInfo[,1]), gene_list))){
    #     errors <- c(errors, "First column of rowInfo does not match EMSet rownames and rowInfo rownames.")
    #   }
    # }
    
    #em.set <- filterByControl(em.set, control = 'Mt', pct.threshold = 100)
    #em.set <- filterByControl(em.set, control = 'Rb', pct.threshold = 100)
    
    #em.set <- filterLowAbundanceGenes(em.set, pct.threshold = 20)
    em.set@log$NormalisationMethod="Deconvolution"
    normcounts(em.set) <- counts(em.set)
    logcounts(em.set) <- counts(em.set)
    names(assays(em.set)) <- c("counts", "normcounts", "logcounts")
    #em.set <-scranNormalise(em.set, quickCluster = TRUE, min.mean = 1e-05)
    
    set.seed(737368)
    em.set <- runPCA(em.set)
    
    em.set <- runCORE(em.set, dims = 20, nres = 40, remove.outliers = FALSE)
    
    end = Sys.time()
    missing_cells <- setdiff(row.names(colData(sce)), row.names(em.set@clusterAnalysis$clusteringMatrix))
    if (length(missing_cells)>0){
      df_missing_cells <- data.frame(cell = missing_cells, ascend = "unassigned")
    }
    res <- data.frame(attr(em.set@clusterAnalysis$distanceMatrix, "Labels"), em.set@clusterAnalysis$clusters)
    colnames(res) <- c("cell", "ascend")
    if (length(missing_cells)>0) {
      res <- rbind(res, df_missing_cells)
    }
    #res <- colData(sce)
    #View(em.set@clusterAnalysis$clusters)
    time = end - start
    write.csv(res, file=folder, sep=",", quote=FALSE, col.names=TRUE, row.names = FALSE)
    cat(paste0(time, " ", dataset, " ascend"), file="time.txt", append=TRUE, sep = "\n")
  
}



