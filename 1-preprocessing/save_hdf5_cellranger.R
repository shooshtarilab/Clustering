#PREPARE FOR CELLRANGER

library("DropletUtils")
library(Matrix)
library(mltools)
library(data.table)

load("lambrets_lung_nonTumour_TPM.RData")
cell.ids <- as.character(colData(sce)$cells)
gene.symbol <- as.character(rowData(sce)$genes)
counts(sce) <- log10(counts(sce)+1)
my.counts <- Matrix(as.matrix(counts(sce)), sparse = TRUE)
write10xCounts(path="./filtered_feature_bc_matrix_lambrechts_lung_nonTumour.h5", x=my.counts, 
               gene.symbol=gene.symbol, barcodes=cell.ids, type="HDF5", overwrite = TRUE, version="2")