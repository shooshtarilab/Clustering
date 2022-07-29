suppressPackageStartupMessages(library(SingleCellExperiment))
source("Functions.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(janeaustenr)
library(forcats)

status <- "Bootstrap-Ensemble" #  "Actual" # 

##################################################################
#######                 Setting Parameters                ########
##################################################################
ref <- "truth"
all.datasets <- c(list.files("./ensemble_results"))

# Additional algorithms - To be added later: "bigScale", "raceid", "simlr"
all.algorithms <- c("clue_3",  "clue_4",  "clue_5",  "clue_6",  "clue_7",  "clue_8",  "clue_9",  
                    "clue_10", "clue_11", "clue_12", "clue_13", "clue_14", "clue_15", "clue_16",
                    "SAME_AIC_3",  "SAME_AIC_4",  "SAME_AIC_5",  "SAME_AIC_6",  "SAME_AIC_7",  
                    "SAME_AIC_8", "SAME_AIC_9",  "SAME_AIC_10", "SAME_AIC_11", "SAME_AIC_12",
                    "SAME_AIC_13", "SAME_AIC_14", "SAME_AIC_15", "SAME_AIC_16",
                    "SAME_BIC_3",  "SAME_BIC_4",  "SAME_BIC_5",  "SAME_BIC_6",  "SAME_BIC_7",  
                    "SAME_BIC_8", "SAME_BIC_9",  "SAME_BIC_10", "SAME_BIC_11", "SAME_BIC_12",
                    "SAME_BIC_13", "SAME_BIC_14", "SAME_BIC_15", "SAME_BIC_16",
                    "SAFE_3",  "SAFE_4",  "SAFE_5",  "SAFE_6",  "SAFE_7",  "SAFE_8",  "SAFE_9",
                    "SAFE_10", "SAFE_11", "SAFE_12", "SAFE_13", "SAFE_14", "SAFE_15", "SAFE_16")
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
  data.path <- paste0("./ensemble_results/", dataset)
  data <- read.csv(data.path, header=T)
  
  for (alg in intersect(all.algorithms, colnames(data))) {
    ref <- "truth"
    print(alg)
    C <- data[, ref]
    K <- data[, alg]
    
    N <- nrow(data)
    cell <- data[,"cell"]
    f.measure[alg, dataset] <- F.Measure.Func(C, K)$F.measure

    homogeneity[alg, dataset] <- Homogeneity.Func(C,K)
    
    majority[alg, dataset] <- Majority.Func(cell, C, K)$Majority
    
    ARI[alg, dataset] <- ARI.Func(C, K)
    
    AMI[alg, dataset] <- AMI.Func(C, K)
    
    VI[alg, dataset] <- VI.Func(C, K)

    x <- c(dataset, alg, metric = "FMeasure", f.measure[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric="Homogeneity", homogeneity[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric="majority", majority[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric = "ARI", ARI[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric = "AMI", AMI[alg, dataset])
    df <- rbind(df, x)
    
    x <- c(dataset, alg, metric = "VI", VI[alg, dataset])
    df <- rbind(df, x)

   write.csv(df, file="./ensemble_results.csv", row.names = FALSE)
  }
  write.csv(df, file="./ensemble_results.csv", row.names = FALSE)
}


df <- read.csv("./ensemble_results.csv")
df <- df[df$metric == "Homogeneity",]
df$method <- as.character(df$method)

for (i in 1:nrow(df)) {
  split_name = strsplit(df[i, "method"], "_")
  df[i, "alg"] <- split_name[[1]][1]
  df[i, "num_inputs"] <- split_name[[1]][2]
}

df <- df[df$metric == "Homogeneity",]
df <- df[df$alg == "clue",]


ggplot(df, aes(num_inputs, value))
