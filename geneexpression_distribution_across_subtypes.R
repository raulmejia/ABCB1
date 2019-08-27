# This program takes a expression matrix and plot the gene expression distribution of a particular
# gene across different predetermined subtypes

################################################################################
# 2. Filtering and preprocessing of data for Pathifier
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/AngelCampos
################################################################################
# Installing and loading the required libraries
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("utils")) {
  install.packages("utils", ask =FALSE)
  library("utils")
}
if (!require("dplyr")) {
  install.packages("dplyr", ask =FALSE)
  library(dplyr)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", ask =FALSE)
  library(ggplot2)
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library(reshape2)
}
if (!require("tidyverse")) {
  install.packages("tidyverse", ask =FALSE)
  library(tidyverse)
}
if (!require("plotly")) {
  install.packages("plotly", ask =FALSE)
  library(plotly)
}


# DATA INPUT from outside
args <- commandArgs(trailingOnly = TRUE)
###########################################
Path_to_your_Matrix<-args[1] # The path to your matrix
# Path_to_your_Matrix<-c("../Data/joined_indicator_METABRIC.txt")
phenoData_path<-args[2] # The path to the information about the samples
# phenoData_path<-c("../Data/joined_indicator_METABRIC.txt")
Path_of_Code<-args[3] # The path to your code
# Path_of_Code<-c("./")
Path_of_Results<-args[4] # # where do you want to save your results?
# Path_of_Results<-c("../Results/RawDistributions/")
Label_for_Results<-args[5] # Label for your results
# Label_for_Results <-"ABCB1_in_METABRIC"
mygene <- args[6]
# mygene <- "ABCB1$"
logtransform <- args[7]
# logtransform <- FALSE
#choose.dir()
#choose.dir(default = "", caption = "Select folder")

######################################
## Loading data and creating folders #
######################################
dir.create(Path_of_Results)
M.matrix<-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)
if (logtransform == TRUE){
  indicator <- M.matrix[1,]
  M.matrix <- M.matrix[-1,]
  M.matrix <- log(M.matrix + 1)
  M.matrix <- rbind( indicator , M.matrix)
  print("Your matrix is being log2 transformed")
} else{
  print("Your matrix wasn't log2 transformed")
}


PhenoData<-read.table(phenoData_path ,header=TRUE , stringsAsFactors = FALSE)

######################################
## Manipulating  the data
######################################

PhenoData_splited <- split(PhenoData, PhenoData$Labels)
names(PhenoData_splited) <- levels(PhenoData$Labels)

#splited_M.matrix <- list()
#for( k in names( PhenoData_splited) ){
#  splited_M.matrix[[k]] <- M.matrix[ ,which(colnames(M.matrix) %in% rownames(PhenoData_splited[[k]]) )]
#}
#names(splited_M.matrix) <- names(PhenoData_splited)
#sum(as.numeric(lapply (splited_M.matrix, function(x) dim(x)[2] )))

MygeneEM <- M.matrix[grep(mygene,rownames( M.matrix)),]

MygeneEMLabels <- data.frame( t(MygeneEM), PhenoData[,1], row.names = rownames(PhenoData) , stringsAsFactors = FALSE) 
colnames( MygeneEMLabels)[2] <- "Labels"
head(MygeneEMLabels) ; colnames(t(MygeneEM))
if (sum(colnames(MygeneEM) == rownames(PhenoData)) / length(colnames(MygeneEM))  == 1 ){
  print("Colnames of the matrix and labels match in the same order")
} else {
  print( "DON'T MATCH colnames of the matrix and Labels in the same order")
}

#########################
## Plotting and saving ##
#########################

pdf(paste0( Path_of_Results,"expresion_of_",mygene,"in",Label_for_Results,".pdf"), height = 7 , width = 10)
ggplot(MygeneEMLabels , aes(t(MygeneEM) , fill = Labels)) + geom_density(alpha = 0.2) + labs(x = "Intesity")  + labs(title = Label_for_Results)
ggplot(data = MygeneEMLabels , aes(t(MygeneEM) , t(MygeneEM))) +  geom_boxplot(aes(colour = Labels)) + labs(y = "Intesity", x = "", title = Label_for_Results)
dev.off()
