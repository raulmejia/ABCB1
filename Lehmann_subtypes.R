################################################################################
# This script receives a expression matrix with HGNC as gene ids
# and returns 
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
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
if (!require("STROMA4")) {
  BiocManager::install("STROMA4", ask =FALSE)
  library("STROMA4")
}
if (!require("breastCancerMAINZ")) {
  BiocManager::install("breastCancerMAINZ", ask =FALSE)
  library("breastCancerMAINZ")
}


# DATA INPUT from outside
args <- commandArgs(trailingOnly = TRUE)
###########################################
Path_to_your_Matrix<-args[1] # The path to your matrix
# Path_to_your_Matrix<-c("../Results/Splited/Matrices/TCGA-log2/Subexpression_matrix_Basal_from_TCGA-log2_.tsv")
#phenoData_path<-args[2] # The path to the information about the samples
# phenoData_path<-c("../Data/joined_indicator_METABRIC.txt")
Path_of_Code<-args[2] # The path to your code
# Path_of_Code<-c("./")
Path_of_Results<-args[3] # # where do you want to save your results?
# Path_of_Results<-c("../Results/Lehmann-STROMA4/")
Label_for_Results<-args[4] # Label for your results
# Label_for_Results <-"Basal-TCGA-log2"
# If you use windows and you have troubles with your paths, please try out this tool: # Path_of_x <-choose.dir()


##### Toy data
library(breastCancerMAINZ)
data(mainz,package = 'breastCancerMAINZ')
head(fData(mainz)[,"Gene.symbol",drop=FALSE])

just.stromal.properties <- assign.properties(ESet=mainz, geneID.column="Gene.symbol",
                                            genelists="Stroma4", n=10, mc.cores=1)

just.lehmann.properties <- assign.properties(ESet=mainz, geneID.column="Gene.symbol",
                                             genelists="TNBCType", n=70, mc.cores=7)
all.properties <- assign.properties(ESet=mainz, geneID.column="Gene.symbol",
                                    genelists=c("Stroma4", "TNBCType"), n=70, mc.cores=7)

property.cols <- grepl( 'properties' , pData(all.properties))
print(apply(pData(all.properties)[, property.cols], 2, table))

detectCores()
head(pData(all.properties))
head(pData(just.lehmann.properties))
fData(mainz)[1:10,1:10]
str(just.stromal.properties)


?[
class(fData(mainz))
attributes(fData(mainz))
colnames(fData(mainz))
?fData
class(mainz)
str(mainz)
attributes(mainz)
names(mainz)
ls()
