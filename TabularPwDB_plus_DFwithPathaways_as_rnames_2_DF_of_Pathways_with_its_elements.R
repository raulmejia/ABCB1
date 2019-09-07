################################################################################
# This program receives a data frame with the definitions of pathways in a tabular way
# (Pathway names as rownames, and genes as its elements), also as input a data frame 
# of particular interest with Pathways as rownames.  The result will be a new DF
# with Pathways as rownames and genes as elements according to the dataframe with
# definitions of pathways
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
Path_to_Ptwydb <- args[1] # The path to the information about the samples
# Path_to_Ptwydb <-c("../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv"); Path_to_your_Matrix <-c("../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt"); Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Top_Deregulted_Pathways/TCGA/") ; Label_for_Results <-"Df_pw_gn_TCGA"
Path_to_your_Matrix <- args[2] # The path to your matrix
# Path_to_your_Matrix <-c("TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt");
Path_of_Code <- args[3] # The path to your code
# Path_of_Code<-c("./")
Path_of_Results<-args[4] # # where do you want to save your results?
# Path_of_Results<-c("../Results/Top_Deregulted_Pathways/TCGA/")
Label_for_Results<-args[5] # Label for your results
# Label_for_Results <-"Df_pw_gn_TCGA"
## NOTE: for windows users #  If you use windows and you have troubles with your paths, please try out this tool: # Path_of_x <-choose.dir()

##########################################################
### Reading the data #####################################
##########################################################
Df_your_Dereg_Pw <- read.table( Path_to_your_Matrix, sep= "\t", header=TRUE, row.names = 1, quote = "" )
PtwDf_defintions <- read.table( Path_to_Ptwydb, sep= "\t", header=TRUE, quote = "" , row.names = 1)

##########################################################
### Getting the gene inside your Deregulated pathways ####
##########################################################
Coincident_positions <- which( rownames( Df_your_Dereg_Pw ) %in% rownames(PtwDf_defintions) )
Df_yourPw_with_genes <- PtwDf_defintions[Coincident_positions, ]






