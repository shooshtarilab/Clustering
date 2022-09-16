.libPaths("/hpf/largeprojects/ccmbio/amahalanabis/tools/R/library")
library("patchwork")
library("Seurat")
library("ggplot2")
library(SingleCellExperiment)

package.version("Seurat")

best <- list()
best[["li_crc_allCells"]] <- "sc3"
best[["li_crc_nonTumour"]] <- "sc3"
best[["chung_breast_allCells"]] <- "ascend"
best[["chung_breast_nonTumour"]] <- "seurat"
best[["darmanis_glioblastoma_allCells"]] <- "ascend"
best[["darmanis_glioblastoma_nonTumour"]] <- "bigScale"
best[["JA_melanoma_allCells"]] <- "cidr"
best[["JA_melanoma_nonTumour"]] <- "rca"
best[["tirosh_melanoma_allCells"]] <- "seurat"
best[["tirosh_melanoma_nonTumour"]] <- "rca"
best[["vangalen_AML_allCells"]] <- "monocle"
best[["vangalen_AML_nonTumour"]] <- "seurat"
best[["peng_pancreatic_allCells"]] <- "seurat"
best[["lambrechts_lung_allCells"]] <- "cellranger"

all.datasets <- c("tirosh_melanoma_malignant", "chung_breast_malignant",
                  "peng_pancreatic_malignant",
                  "vangalen_AML_malignant",
                  "JA_melanoma_malignant", "lambrechts_lung_malignant",
                  "darmanis_glioblastoma_malignant", "li_crc_malignant")

all.algorithms <- c("altAnalyze",	"ascend", "bigScale",
                    "cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",
                    "phenograph",	"raceid", "rca",	"sc3",	"scran",	"seurat",
                    "simlr", "sincera",	"tscan")

##### TOP ROW OF FIGURES - UMAP COLOURED BY TOP PERFORMING ALGORITHM ACROSS ALL DATASETS

dataset <- "vangalen_AML_malignant"
truth_dir <- "/Users/alainamahalanabis/Documents/clustering/GEO_datasets/vangalen_AML_malignant/vangalen_aml_patients.tsv"

#prefix_dir <- "/hpf/largeprojects/ccmbio/amahalanabis/silhouette/"
prefix_dir <- "/Users/alainamahalanabis/Documents/clustering/GEO_datasets/"
truth_dir <- "/Users/alainamahalanabis/Documents/clustering/analysis_scripts/Clustering_Results/"
count = 0
plt <- list()

for (dataset in all.datasets) {
  print(dataset)
  rdata_dir <- paste0(prefix_dir, dataset, "_counts.RData")
  print (rdata_dir)
  load(rdata_dir)
  #truth_dir <- paste0("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/Clustering_Results/", dataset, ".csv")
  truth <- read.csv(truth_dir, header=TRUE, sep = "\t")
  truth <- truth[truth$cell %in% colnames(counts(sce)),]
  seurat.object <- CreateSeuratObject(counts=counts(sce))
  seurat.object <- NormalizeData(seurat.object, verbose=F)
  seurat.object <- FindVariableFeatures(seurat.object, verbose=F)
  seurat.object <- ScaleData(seurat.object, features = rownames(seurat.object), verbose=F)
  seurat.object <- RunPCA(seurat.object, features = VariableFeatures(seurat.object), verbose=F)
  seurat.object <- RunUMAP(seurat.object, dims = 1:20, verbose=F)
  seurat.object <- RunTSNE(seurat.object, check_duplicates = FALSE, dims = 1:20, verbose=F)
  #truth_dir <- paste0("/Users/alainamahalanabis/Documents/clustering/analysis_scripts/Clustering_Results/", dataset, ".csv")
  #truth <- read.csv(truth_dir, header=TRUE)
  #truth <- truth[truth$cell %in% colnames(counts(sce)),]
  df <- as.data.frame(truth[c("truth", "cell")], row.names=as.character(truth$cell))
  
  df <- as.data.frame(truth[c("patient", "cell")], row.names=as.character(truth$cell))
  count = count + 1  

  seurat.object <- AddMetaData(seurat.object, df)
  title <- paste0(dataset, " ", "truth_patient")
    plt[[count]] <- DimPlot(seurat.object, reduction = "tsne", group.by ="patient") +
    ggtitle(title)
 
  df <- as.data.frame(truth[c("truth", "patient")], row.names=as.character(truth$cell))
  
  count = count + 1 
  seurat.object <- AddMetaData(seurat.object, df)
  title <- paste0(dataset, " ", "Patient id")
    plt[[count]] <- DimPlot(seurat.object, reduction = "tsne", group.by = "patient") +
    ggtitle(title)
  
}

plt[[4]]
ggsave(filename=paste0(all.datasets[1],"_tsne_png.pdf"), width=4.5*2, height=2.8*2)
dev.off()

plt[[1]] + plt[[2]]
ggsave(filename=paste0(all.datasets[6],"_tsne_png.pdf"), width=4.5*4, height=2.8*3)
dev.off()

plt[[3]] + plt[[4]]
ggsave(filename=paste0(all.datasets[2],"_tsne_png.pdf"), width=4.5*2, height=2.8)
dev.off()

plt[[5]] + plt[[6]]
ggsave(filename=paste0(all.datasets[3],"_tsne_png.pdf"), width=4.5*4, height=2.8*3)
dev.off()

plt[[7]] + plt[[8]]
ggsave(filename=paste0(all.datasets[4],"_tsne_png.pdf"), width=4.5*4, height=2.8*3)
dev.off()

plt[[9]] + plt[[10]]
ggsave(filename=paste0(all.datasets[5],"_tsne_png.pdf"), width=4.5*4, height=2.8*3)
dev.off()

plt[[3]]
ggsave(filename=paste0(all.datasets[8],"_tsne_png.pdf"), width=4.5*2, height=2.8*2)
dev.off()

# plt[[1]] + plt[[2]] + plt[[3]] + plt[[4]]
# ggsave(filename=paste0(all.datasets[1],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
# plt[[5]] + plt[[6]] + plt[[7]] + plt[[8]]
# ggsave(filename=paste0(all.datasets[2],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
# plt[[9]] + plt[[10]] + plt[[11]] + plt[[12]]
# ggsave(filename=paste0(all.datasets[3],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
# plt[[13]] + plt[[14]] + plt[[15]] + plt[[16]]
# ggsave(filename=paste0(all.datasets[4],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
# plt[[17]] + plt[[18]] + plt[[19]] + plt[[20]]
# ggsave(filename=paste0(all.datasets[5],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
# plt[[21]] + plt[[22]] + plt[[23]] + plt[[24]]
# ggsave(filename=paste0(all.datasets[6],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
# plt[[25]] + plt[[26]] + plt[[27]] + plt[[28]]
# ggsave(filename=paste0(all.datasets[7],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
# plt[[29]] + plt[[30]] + plt[[31]] + plt[[32]] 
# ggsave(filename=paste0(all.datasets[8],"_tsne_png.pdf"), width=4.5*4, height=2.8, type="pdf") 
# dev.off()
