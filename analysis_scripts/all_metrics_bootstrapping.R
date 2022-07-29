suppressPackageStartupMessages(library(SingleCellExperiment))
source("Functions.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(janeaustenr)
library(forcats)

status <- "Bootstrap-Ensemble" #  "Actual" # 
setwd("/Users/alainamahalanabis/Documents/clustering/analysis_scripts")

##################################################################
#######                 Setting Parameters                ########
##################################################################
ref <- "truth"
all.datasets <-  c("li_crc_malignant","chung_breast_malignant", "JA_melanoma_malignant",
                   "tirosh_melanoma_malignant")

# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
all.algorithms <- c("altAnalyze",	"ascend", "bigScale",
                    "cellranger",	"cidr",	"countClust",	"monocle",	"pcaReduce",
                    "phenograph",	"raceid", "rca",	"sc3",	"scran",	"seurat",
                    "simlr", "sincera",	"tscan")

##################################################################
#######               Calculating 6 metrics               ########
##################################################################
f.measure <- f.measure.low.CI <- f.measure.high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))
ARI <- ARI.low.CI <- ARI.high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))
VI <- VI.low.CI <- VI.high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))
AMI <- AMI.low.CI <- AMI.high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))
majority <- majority.low.CI <- majority.high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))
homogeneity <- homogeneity.low.CI <- homogeneity.high.CI <- matrix(0, nrow=length(all.algorithms), ncol=length(all.datasets), dimnames=list(all.algorithms, all.datasets))

df <- c()
for (dataset in all.datasets) {
  print(dataset)
  data.path <- paste0("./Clustering_Results/", dataset, ".csv")
  data <- read.csv(data.path, header=T)
  
  for (alg in intersect(all.algorithms, colnames(data))) {
    print(alg)
    C <- data[, ref]
    K <- data[, alg]
    
    N <- nrow(data)
    cell <- data[,"cell"]
    bootstrap.func <- function() {
      indx <- sample(1:N, size=N, replace=T)
      C1 <- C[indx]
      K1 <- K[indx]
      cell1 <- cell[indx]
      F.Measure = F.Measure.Func(C1, K1)$F.measure
      Homogeneity = Homogeneity.Func(C1,K1)
      #Homogeneity=1
      Majority = Majority.Func(cell1, C1, K1)$Majority
      ARI = ARI.Func(C1, K1)
      AMI = AMI.Func(C1, K1)
      VI = VI.Func(C1, K1)
      ret <- c(F.Measure, Homogeneity, Majority, ARI, AMI, VI)
      return(ret)
    }
    rslt <- c()
    for (i in 1:10000) {
      #print (i)
      ret <- bootstrap.func()
      if (!(any(is.na(ret)))) {
          rslt<-rbind(rslt, ret)
      }
    }
    #write.csv(rslt, file="rslt.csv")
    # Measure 95% confidence intervals using bootstrapping
    #print("columns in rslt")
    #print(ncol(rslt))   
    f.measure.CI <- quantile(rslt[,1], probs = c(0.025, 0.5, 0.975))
    f.measure.low.CI[alg, dataset] <- round(f.measure.CI[1], 2)
    f.measure.high.CI[alg, dataset] <- round(f.measure.CI[3], 2)
    f.measure[alg, dataset] <- round(f.measure.CI[2], 2)
    
    homogeneity.CI <- quantile(rslt[,2], probs = c(0.025, 0.5, 0.975))
    homogeneity.low.CI[alg, dataset] <- round(homogeneity.CI[1], 2)
    homogeneity.high.CI[alg, dataset] <- round(homogeneity.CI[3], 2)
    homogeneity[alg, dataset] <- round(homogeneity.CI[2], 2)
    
    majority.CI <- quantile(rslt[,3], probs = c(0.025, 0.5, 0.975))
    majority.low.CI[alg, dataset] <- round(majority.CI[1], 2)
    majority.high.CI[alg, dataset] <- round(majority.CI[3], 2)
    majority[alg, dataset] <- round(majority.CI[2], 2)
    
    ARI.CI <- quantile(rslt[,4], probs = c(0.025, 0.5, 0.975))
    ARI.low.CI[alg, dataset] <- round(ARI.CI[1], 2)
    ARI.high.CI[alg, dataset] <- round(ARI.CI[3], 2)
    ARI[alg, dataset] <- round(ARI.CI[2], 2)
    
    AMI.CI <- quantile(rslt[,5], probs = c(0.025, 0.5, 0.975))
    AMI.low.CI[alg, dataset] <- round(AMI.CI[1], 2)
    AMI.high.CI[alg, dataset] <- round(AMI.CI[3], 2)
    AMI[alg, dataset] <- round(AMI.CI[2], 2)
    
    VI.CI <- quantile(rslt[,6], probs = c(0.025, 0.5, 0.975))
    VI.low.CI[alg, dataset] <- round(VI.CI[1], 2)
    VI.high.CI[alg, dataset] <- round(VI.CI[3], 2)
    VI[alg, dataset] <- round(VI.CI[2], 2)

    x <- c(dataset, alg, metric = "FMeasure", f.measure[alg, dataset], f.measure.low.CI[alg, dataset], f.measure.high.CI[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric="Homogeneity", homogeneity[alg, dataset], homogeneity.low.CI[alg, dataset], homogeneity.high.CI[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric="majority", majority[alg, dataset], majority.low.CI[alg, dataset], majority.high.CI[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric = "ARI", ARI[alg, dataset], ARI.low.CI[alg, dataset], ARI.high.CI[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric = "AMI", AMI[alg, dataset], AMI.low.CI[alg, dataset], AMI.high.CI[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric = "VI", VI[alg, dataset], VI.low.CI[alg, dataset], VI.high.CI[alg, dataset])
    df <- rbind(df, x)

   write.csv(df, file="dataframe_results/all_results.csv", row.names = FALSE)
  }
  write.csv(df, file="dataframe_results/all_results.csv", row.names = FALSE)
}
