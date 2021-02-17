#PACKAGES TO BE INSTALLED
install.packages("clue")
install.packages("rlist")
install.packages("relations")

#LOAD AND SOURCE THESE PACKAGES
library("clue")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/AAA.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/partition.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/ensemble.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/consensus.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/dissimilarity.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/agreement.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/utilities.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/membership.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/registration.R")
source("/Users/alainamahalanabis/Documents/clustering/clue/R/objects.R")
library("relations")
library(rlist)


run_clue <- function(input_dataframe, tools) {
  final <- input_dataframe
  res <- data.frame(cell = final$cell)
  final <- final[tools]
  clue_input <- list()
  for (col in colnames(final)) {
    print (col)
    tmp <-list(cell = final$cell , class_ids=final[,col])
    class(tmp) <- "cluster"
    clue_input <- list.append(clue_input, tmp)
  }
  
  final$cell = NULL
  names(clue_input) <- colnames(final)
  
  ens<-cl_ensemble(list = clue_input)
  
  consensus <- cl_consensus(ens, "HE") 
  
  res <- cbind(res, cluster = c(cl_class_ids(consensus)))
  colnames(res)[ncol(res)] = "clue"
  colnames(res) <- c("cell", "clue")
  write.csv(res, file="clue.csv", row.names=FALSE)
  
}
