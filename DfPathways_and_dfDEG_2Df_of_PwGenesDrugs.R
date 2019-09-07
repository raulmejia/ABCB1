################################################################################
# This program receives a data frame of KEGG pathways, a dataframe of KEGG
# pathways definitions and retrieves a data frame of Pathway, Gene, Drugs
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

# DATA INPUT from outside
args <- commandArgs(trailingOnly = TRUE)
###########################################
Path_to_your_Matrix<-args[1] # The path to your matrix
# Path_to_your_Matrix<-c("../Results/Pathifier/Basal_under_25_stbl_10_threshold_3_75_median_PDSz_ordered_matrix_top_40.txt")
Path_to_Ptwydb <-args[2] # The path to the information about the samples
# Path_to_Ptwydb <-c("../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv")
Path_of_Code <-args[3] # The path to your code
# Path_of_Code<-c("./")
Path_of_Results<-args[4] # # where do you want to save your results?
# Path_of_Results<-c("../Results/Drugs/")
Label_for_Results<-args[5] # Label for your results
# Label_for_Results <-"Df_pw_gn_drug_"
# If you use windows and you have troubles with your paths, please try out this tool: # Path_of_x <-choose.dir()

### Reading data
myDeg_pathwaysDf <- read.table( Path_to_your_Matrix, sep= "\t", header=TRUE, row.names = 1, quote = "" )
PtwDf_defintions <- read.table( Path_to_Ptwydb, sep= "\t", header=TRUE, quote = "" )

