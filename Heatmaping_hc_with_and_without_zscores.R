################################################################################
# 2. Filtering and preprocessing of data for Pathifier
## Author: Miguel Angel Garcia Campos github: https://github.com/AngelCampos
################################################################################
## MICHAL's  #_# METABRIC flavor
# DATA INPUT
# captura argumentos de la linea de comandos
args <- commandArgs(trailingOnly = TRUE)

###########################################
Path_to_your_Matrix<-args[1] # The path to your matrix
# Path_to_your_Matrix<-c("../Results/Splited/Matrices/TCGA/Subexpression_matrix_Basal_from_TCGA_.tsv_Only_ABC_transporters_genes.tsv")
# Path_to_your_Matrix <- choose.files() # choose.dir()
Path_of_Code<-args[2] # The path to your code
# Path_of_Code<-c("./")
Path_of_Results<-args[3] # # where do you want to save your results?
# Path_of_Results<-c("../Results/Clusterig/")
# Path_of_Results <- choose.dir()
# choose.dir(default = "", caption = "Select folder")
Labels <-args[4] # Label for your results
# Labels <- "Basal_Only_ABC_transporters"
###############################################################################
### Installing and/or loading required packages
###############################################################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("som")) {
  install.packages("som", dependencies = TRUE)
  library(som)
}
if (!require("utils")) {
  install.packages("utils", ask =FALSE)
  library("utils")
}

source(paste0(Path_of_Code,"plot_raw_Matrix_png.R"))
source(paste0(Path_of_Code,"plot_raw_Matrix_png_euclidean_wardD2.R"))

########################
### Reading the data and creating the folder for the results
########################
dir.create(Path_of_Results , recursive = TRUE )
M.matrix<-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)

############################
####### Preparing the data
############################
# calculating z-scores
M.matrix_zscores<-som::normalize(M.matrix, byrow = TRUE)
colnames(M.matrix_zscores) <- colnames(M.matrix)

################################
###### Ploting
###############################

color_labels_vector <- rep("red",dim(M.matrix)[2]) # you can replace the color vector
M.matrix <- as.matrix(M.matrix)
class(M.matrix)
plot_raw_Matrix_png(M.matrix,paste0(Labels),paste0(Path_of_Results,Labels,c("_heatmap_clustering_euclidean_complete.png")),color_labels_vector)
plot_raw_Matrix_png(M.matrix_zscores,paste0(Labels),paste0(Path_of_Results,Labels,c("_Z-scores_heatmap_clustering_euclidean_complete.png")),color_labels_vector)

plot_raw_Matrix_png_ecuclidean_wardD2(M.matrix,paste0(Labels),paste0(Path_of_Results,Labels,c("_heatmap_clustering_euclidean_wardD2.png")),color_labels_vector)
plot_raw_Matrix_png_ecuclidean_wardD2(M.matrix_zscores,paste0(Labels),paste0(Path_of_Results,Labels,c("_Z-scores_heatmap_clustering_euclidean_wardD2.png")),color_labels_vector)

