Cell clustering analysis 

These are scripts accompanying our paper on automated methods for clustering cell types in scRNAseq data from the tumour microenvironment.

preprocessing: Contains preprocessing scripts. Cellranger requires hdf5 file format as input. SingleCellExperiment Object is created from all data matrixes.
clustering_algorithms: Contains clustering algorithms that we used in our analysis
analysis_scripts: Contains scripts for evaluating clustering performance. After we run clustering algorithms on all datasets, we can run analysis scripts to compute the ARI, F-measure, homogeneity, etc..
data: Contains intermediate data files that are used for generating paper figures. 

