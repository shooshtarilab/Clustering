#devtools::install_github("xu-lab/SINCERA")
library(SINCERA)
library(SingleCellExperiment)

main <- function(filename) {
  
  load(filename)
  counts(sce) <- log10(counts(sce)+1)
  
  dataset = paste0(strsplit(filename, "_")[[1]][1], "_", strsplit(filename, "_")[[1]][2], "_", strsplit(filename, "_")[[1]][3])
  folder = paste0("./", dataset, "/sincera.csv")
  start = Sys.time()

  sc <- construct(exprmatrix = data.frame(counts(sce)), samplevector = colnames(counts(sce)))
  
  # Filter out non-expressed genes
  sc <- prefilterGenes(sc, pergroup=FALSE, min.expression=1, min.cells=1, min.samples=1)
  
  # Perform z-score scaling
  sc <- normalization.zscore(sc, pergroup=FALSE)
  
  # do PCA using all genes
  sc <- doPCA(sc, genes=NULL, use.fast = T)
  
  min.samples <- 6
  obj <- prefilterGenes(sc, pergroup=FALSE, min.expression=1, min.cells=10, min.samples=min.samples)
  # genes with at least 0.7 specificity in at least 6 samples
  obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=0.7)
  
  # set the selected genes for clustering
  sc <- setGenesForClustering(sc, value=getGenesForClustering(obj))
  
  sc <- cluster.assignment(sc)
  
  end = Sys.time()
  time = end - start
  res <- data.frame(cell=colnames(sce), sincera = getCellMeta(sc, "CLUSTER"))
  
  write.csv(res, file=folder, sep=",", quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " sincera"), file="time.txt", append=TRUE, sep = "\n")
  
}

datasets <- c("JA_melanoma_nonTumour_TPM.RData",
              "tirosh_melanoma_nonTumour_TPM.RData",
              "darmanis_glioblastoma_nonTumour_TPM.RData")

for (i in datasets) {
  main(i)
}



library(SINCERA)

# load IPF data, which contains two data frames
# - the ipf.exprmatrix contains the TPM values of 540 single cell transcriptomes
# - the ipf.cells contains the cell sample and condition information
data(IPF)

ls()

dim(ipf.exprmatrix)

head(ipf.cells)

table(ipf.cells$Diagnosis) # normal or ipf cells

table(ipf.cells$DataID) # cell sample


# The analysis starts with running the construct function to create an R S4 object, which will hold all the data and analysis results.
# The function takes expression matrix and cell sample information as input
sc <- construct(exprmatrix=ipf.exprmatrix, samplevector=ipf.cells$DataID)

# After contruction, you can use setCellMeta to add more cell information to sincera
# such as adding the condition information (CONTROL or IPF) of cells into sincera object
sc <- setCellMeta(sc, name="CONDITION", value=ipf.cells$Diagnosis)

# use getCellMeta functoin to assess a specific meta data
# In most of the SINCERA functions, cell grouping will be based on the GROUP meta data
# The GROUP meta was initialized to sample information during the construction
table(getCellMeta(sc, name="GROUP"))

# Identify and remove low quality cells.
# The key parameters of running this function include: “min.expression”,
# which specifies the minimum expression value for a gene to be considered as expressed,
# and “min.genes”, which specifies the lower bound of the number of expressed genes in a cell.
sc <- filterLowQualityCells(sc, min.expression=1, min.genes=1000, do.plot = T)

# Set the minimum value of expression to 0.01
sc <- expr.minimum(sc, value=0.01)

# Run batch.analysis function to generate plots that may help identify potential batch differences.
sc <- batch.analysis(sc, analysis=c("distribution"), min.expression=1)

# Filter out non-expressed genes
sc <- prefilterGenes(sc, pergroup=FALSE, min.expression=1, min.cells=1, min.samples=1)

# Perform z-score scaling
sc <- normalization.zscore(sc, pergroup=FALSE)

# do PCA using all genes
sc <- doPCA(sc, genes=NULL, use.fast = T)

# plot the standard deviation of PCA components
plotPCASD(sc, num.pcs = 20)

# do tSNE using PCA components
sc <- doTSNE(sc, genes=NULL, dims = 1:5, use.fast = T)

# plot cells in tSNE spaces
plotRDS(sc, feature.type="tsne")




## first iteration


#In the first iteration, we selected genes for clustering using the following criteria
# genes expressed in at least 10 cells in at least 6 samples
min.samples <- 6
obj <- prefilterGenes(sc, pergroup=TRUE, min.expression=1, min.cells=10, min.samples=min.samples)
# genes with at least 0.7 specificity in at least 6 samples
obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=0.7)

# set the selected genes for clustering
sc <- setGenesForClustering(sc, value=getGenesForClustering(obj))


# use gap statistics to determine the number of clusters
if (FALSE) {
  x <- getExpression(sc, scaled=T, genes=getGenesForClustering(sc))
  
  library(cluster)
  cordist <- function(y) as.dist((1-cor(t(y)))/2)
  hclustForGap <- function(y, k) list(cluster=cutree(hclust(cordist(y), method = "average"),k=k))
  gapstats <- clusGap(t(x),
                      FUN = hclustForGap,
                      K.max = 10,
                      B = 100)
}

library(cluster)
print(ipf.gapstats$iter1) # the gap statistics suggested 5 clusters

# use the default HC algorithm to identify 5 clusters
# the clustering results were saved to the "CLUSTER" metadata
# if update.cellgroup is TRUE, the GROUP meta data will also be updated
sc <- cluster.assignment(sc, k=5)

res <- data.frame(cell=colnames(sce), sincera = getCellMeta(sc, "CLUSTER"))

