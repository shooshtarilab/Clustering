#install.packages("/Users/alainamahalanabis/Downloads/cellrangerRkit-master")
devtools::install_github("hb-gitified/cellrangerRkit")

library(cellrangerRkit)
library(SingleCellExperiment)
library(Matrix)
matrix_dir = "./NSCLC_Lambrechts/"
matrix.path <- paste0(matrix_dir, "matrix.mtx")
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
mat <- readMM(file = matrix.path)
gene_info <- read.delim(features.path, stringsAsFactors = FALSE, 
                        sep = "\t", header = FALSE)
gene_info <- data.frame(mat.gene_name = gene_info$V2)


#PENG PANCREATIC MALIGNANT

library(scater)
mat <- readRDS("vangalen_counts.rds")
mat <- read.csv("/Users/alainamahalanabis/Documents/clustering/GEO_datasets/peng_pancreatic_malignant/Data.csv", stringsAsFactors = F)
cells <- read.csv("/Users/alainamahalanabis/Documents/clustering/GEO_datasets/peng_pancreatic_malignant/celllist.tsv", header = F)
mat <- cbind(cells, mat)
colnames(mat)[1] <- "gene"
mat$X=NULL
mat$gene = as.character(mat$gene)
rownames(mat) <- mat$gene
mat$gene=NULL
mat <- mat[rownames(mat) %in% truth$cell,]
mat$genes=NULL
mat <- t(mat)

truth <- read.csv("./peng_pancreatic_allCells/truth.csv")
truth <- truth[truth$category == "Tumour",]
write.csv(truth, "./peng_pancreatic_malignant/truth.csv")

for (i in 1:nrow(mat)) {
  mat[i, "gene"] <- strsplit(mat[i,"gene"], "_")[[1]][2]
}

mat <- mat[!duplicated(mat$gene),]
rownames(mat) <- mat$gene
mat <- mat[,colnames(mat) %in% truth$cell]
#rownames(mat) <- featureData$mat.gene_name
mat <- mat[rownames(mat) %in% my_gene$gene_name,]
featureData <- data.frame(mat.gene_name = rownames(mat))
#mat <- apply(mat, 2, as.numeric)
barcodes <- colnames(mat)
row.names(mat) = NULL
pd = data.frame(id = barcodes, row.names = barcodes)
gene.symb <- featureData$mat.gene_name
rownames(mat) <- featureData$mat.gene_name

gene_bc_matrix <- newGeneBCMatrix(mat = mat, fd = featureData, pd = pd)
sce <- SingleCellExperiment(
  assays = list(counts = as(exprs(gene_bc_matrix), "matrix")), colData = as(phenoData(gene_bc_matrix), 'data.frame'), rowData=as(featureData(gene_bc_matrix), 'data.frame'))
save(sce, file="patel_glioblastoma_malignant_counts.RData")

#JA-melanoma
mat <- read.csv("/Users/alainamahalanabis/Documents/clustering/JA-melanoma/GSE115978_tpm.csv")
truth <- read.csv("melanoma_jerby_cellanno.csv")
truth <- truth[truth$Cell_Type == 0,]
colnames(mat)[1]<-"gene_name"

rownames(mat) <- mat[,1]
featureData <- data.frame(mat.gene_name = rownames(mat))
mat <- mat[,colnames(mat) %in% truth$Cell_ID]
mat <- apply(mat, 2, as.numeric)
#mat <- log((mat+1), 10)
mat <- data.frame(mat)
barcodes <- colnames(mat)
row.names(mat) = NULL
pd = data.frame(id = barcodes, row.names = barcodes)
gene.symb <- featureData$mat.gene_name

gene_bc_matrix <- newGeneBCMatrix(mat = mat, fd = featureData, pd = pd)
sce <- SingleCellExperiment(
  assays = list(counts = as(exprs(gene_bc_matrix), "matrix")), colData = as(phenoData(gene_bc_matrix), 'data.frame'), rowData=as(featureData(gene_bc_matrix), 'data.frame'))
save(sce, file="JA_melanoma_relabelled.RData")

#TIROSH MELANOMA
mat <- read.table("GSE72056_melanoma_single_cell_revised_v2.txt", sep='\t', stringsAsFactors = FALSE)
test <- mat[,mat[3,]!=0 && mat[4,]!=0]
truth <- read.csv("melanoma_tirosh_cellanno.csv")
truth <- truth[truth$Cell_Type == 0,]
labels <- read.csv("./Clustering_Results_allCells/truth.csv")
truth <- merge(truth, labels, by.x = "Cell_ID", by.y = "cell")
truth$Cell_Type = NULLs
colnames(truth)[1] = "cell"
write.csv(truth, file = "Clustering_Results_relabelled/truth.csv")

colnames(mat) <- mat[1,]
View(mat)
#mat <- mat[,mat[3,] == 1]
#mat$gene <- gene
mat <- mat[-c(1:4),]
colnames(mat)[1] = "gene_name"

featureData <- data.frame(mat$gene_name)
mat <- mat[,colnames(mat) %in% truth$Cell_ID]
row.names(mat) <- 1:nrow(mat)
barcodes = colnames(mat)[2:ncol(mat)]
mat$gene_name=NULL
mat <- sapply(mat, as.numeric)
pd = data.frame(id = barcodes, row.names = barcodes)
row.names(mat) = NULL

gene_bc_matrix <- newGeneBCMatrix(mat = mat, fd = featureData, pd = pd)
sce <- SingleCellExperiment(
  assays = list(counts = as(exprs(gene_bc_matrix), "matrix")), colData = as(phenoData(gene_bc_matrix), 'data.frame'), rowData=as(featureData(gene_bc_matrix), 'data.frame'))
save(sce, file="tirosh_melanoma_allCells_counts.RData")

#CHUNG BREAST

mat <- read.table("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", sep='\t', stringsAsFactors = FALSE)
truth <- read.csv("breastcancer_chung_cellanno.csv")
truth <- truth[truth$Cell_Type == 0,]
colnames(mat) <- mat[1,]
featureData <- data.frame(mat$gene_id ,mat$gene_name)
featureData <- featureData[-1,]
row.names(mat) <- 1:nrow(mat)
mat <-mat[-1,]
mat <- mat[, colnames(mat) %in% truth$Cell_ID]
barcodes = colnames(mat)[1:ncol(mat)]
mat <- sapply(mat, as.numeric)
mat <- log10(mat+1)
pd = data.frame(id = barcodes, row.names = barcodes)
row.names(mat) = NULL

#PURAM HEAD AND NECK

mat <- read.csv("GSE103322_HNSCC_all_data.txt", stringsAsFactors = FALSE, sep = "\t")
truth <- read.csv("hnscc_puram_cellanno.csv")
truth <- truth[truth$Cell_Type == 0,]
#truth_file <- t(rbind(colnames(mat), mat["non-cancer cell type", ]))
rownames(truth_file) = NULL
colnames(truth_file) <- c("cell", "truth")
truth_file <- data.frame(truth_file)
#truth_file <- truth_file[truth_file$truth != 0,]
truth_file <- truth_file[order(truth_file$truth),]
truth_file$truth <- ordered(truth_file$truth, labels = c("-Fibroblast", "tumour", "B cell", "Dendritic",
         "Endothelial", "Fibroblast", "Macrophage", "Mast", "myocyte", "T cell"))
write.csv(truth_file, file="./Clustering_Results_allCells/truth.csv")

mat <- mat[,colnames(mat) %in% truth$Cell_ID]
mat <- mat[-c(1:5),]
featureData <- data.frame(mat.gene_name = rownames(mat))
mat <- apply(mat, 2, as.numeric)
mat <- data.frame(mat)
barcodes <- colnames(mat)
row.names(mat) = NULL
pd = data.frame(id = barcodes, row.names = barcodes)

gene_bc_matrix <- newGeneBCMatrix(mat = mat, fd = featureData, pd = pd)
sce <- SingleCellExperiment(
  assays = list(counts = as(exprs(gene_bc_matrix), "matrix")), colData = as(phenoData(gene_bc_matrix), 'data.frame'), rowData=as(featureData(gene_bc_matrix), 'data.frame'))
save(sce, file="puram_headNeck_relabelled.RData")

#TIROSH MELANOMA DATASET

truth <- read.csv("../GEO_datasets/li_crc_allCells/truth.csv", sep = ",")
truth <- data.table(truth)
truth[truth =="B cell", category := "Immune"]
truth[truth =="CAF", category := "Fibroblast"]
truth[truth =="T cell", category := "Immune"]
truth[truth =="Endo", category := "Endothelial"]
truth[truth =="tumour", category := "Tumour"]
truth[truth =="Macro", category := "Immune"]
truth[truth =="NK", category := "Immune"]

write.csv(truth, "../GEO_datasets/li_crc_allCells/truth.csv", row.names = FALSE)

#JA MELANOMA DATASET

truth <- read.csv("../GEO_datasets/JA_melanoma_nonTumour/truth.csv", sep = ",")
truth <- data.table(truth)
truth[truth =="B cell", category := "Immune"]
truth[truth =="NK", category := "Immune"]
truth[truth =="CAF", category := "Fibroblast"]
truth[truth =="T cell", category := "Immune"]
truth[truth =="Endo", category := "Endothelial"]
truth[truth =="Mal", category := "Tumour"]
truth[truth =="T.CD4", category := "Immune"]
truth[truth =="T.CD8", category := "Immune"]
truth[truth =="Macrophage", category := "Immune"]

write.csv(truth, "../GEO_datasets/JA_melanoma_nonTumour/truth.csv", row.names = FALSE)

#CHUNG BREAST DATASET

truth <- read.csv("../GEO_datasets/chung_breast_nonTumour/truth.csv", sep = ",")
truth <- data.table(truth)
truth[truth =="Bcell", category := "Immune"]
truth[truth =="Myeloid", category := "Immune"]
truth[truth =="Stromal", category := "Stromal"]
truth[truth =="Tcell", category := "Immune"]
truth[truth =="Tumour", category := "Tumour"]
write.csv(truth, "../GEO_datasets/chung_breast_nonTumour/truth.csv", row.names = FALSE)

#PENG PANCREATIC DATASET
truth <- read.csv("../GEO_datasets/peng_pancreatic_allCells/truth.csv", sep = ",")
truth <- data.table(truth)

truth[truth =="Stellate cell", category := "Fibroblast"]
truth[truth =="T cell", category := "Immune"]
truth[truth =="Acinar cell", category := "Epithelial"]
truth[truth =="Ductal cell type 1", category := "Epithelial"]
truth[truth =="Ductal cell type 2", category := "Tumour"]
truth[truth =="Endocrine cell", category := "Secretory"]
truth[truth =="Endothelial cell", category := "Endocrine"]
truth[truth =="Fibroblast cell", category := "Fibroblast"]
truth[truth =="Macrophage cell", category := "Immune"]
write.csv(truth, "../GEO_datasets/peng_pancreatic_allCells/truth.csv", row.names = FALSE)

#LAMBRECHTS LUNG DATASET
truth <- read.csv("../GEO_datasets/lambrechts_lung_nonTumour/truth.csv", sep = ",")
truth <- data.table(truth)
truth$category=""

truth[truth =="Flat alveolar type 1 (AT1) cells", category := "Epithelial"]
truth[truth =="Alveolar cell", category := "Epithelial"]
truth[truth =="Cuboidal alveolar type 2 (AT2) cells", category := "Epithelial"]
truth[truth =="Basal cells", category := "Epithelial"]
truth[truth =="B cells", category := "Immune"]
truth[truth =="Cancer cells", category := "Tumour"]
truth[truth =="Cancer cells pt 2", category := "Tumour"]
truth[truth =="Cancer cells pt 3", category := "Tumour"]
truth[truth =="Cancer cells pt 4", category := "Tumour"]
truth[truth =="Cancer cells pt 5", category := "Tumour"]
truth[truth =="Dendritic cells", category := "Immune"]
truth[truth =="Lower quality endothelial cell", category := "Endothelial"]
truth[truth =="Monocyte-derived dendritic cells", category := "Immune"]
truth[truth =="Plasmacytoid dendritic cells", category := "Immune"]
truth[truth =="Epithelial cells", category := "Epithelial"]
truth[truth =="COL12A1-expressing fibroblasts", category := "Fibroblast"]
truth[truth =="GABARAP-expressing fibroblasts", category := "Fibroblast"]
truth[truth =="PLA2G2A-expressing fibroblasts", category := "Fibroblast"]
truth[truth =="Fibroblasts", category := "Fibroblast"]
truth[truth =="Normal  lung fibroblasts ", category := "Fibroblast"]
truth[truth =="Erythroblasts", category := "Immune"]
truth[truth =="Granulocytes", category := "Immune"]
truth[truth =="Langerhans cells", category := "Immune"]
truth[truth =="Lymphatic EC", category := "Endothelial"]
truth[truth =="Macrophages", category := "Immune"]
truth[truth =="Mast cells", category := "Immune"]
truth[truth =="Natural killer cells", category := "Immune"]
truth[truth =="Secretory club cells", category := "Epithelial"]
truth[truth =="Regulatory T cells", category := "Immune"]
truth[truth =="T cells", category := "Immune"]
truth[truth =="COPD-injured alveolar cells", category := "Epithelial"]
truth[truth =="Endothelial cell", category := "Endothelial"]
truth[truth =="Respiratory epithelial cells", category := "Epithelial"]
truth[truth =="TFPI2-expressing fibroblasts", category := "Fibroblast"]
truth[truth =="Tumour endothelial cell", category := "Tumour"]
truth[truth =="CD8+ T cells", category := "Immune"]
truth[truth =="MALT B cells", category := "Immune"]
truth[truth =="Follicular B cells", category := "Immune"]
truth[truth =="COL4A2-expressing fibroblasts", category := "Fibroblast"]

write.csv(truth, "../GEO_datasets/lambrechts_lung_nonTumour/truth.csv", row.names = FALSE)

truth <- data.frame(truth$cell, truth$category)
colnames(truth) <- c("cell", "truth")
truth <- truth[(truth$cell %in% colnames(counts(sce))),]

truth$truth <- as.character(truth$truth)

