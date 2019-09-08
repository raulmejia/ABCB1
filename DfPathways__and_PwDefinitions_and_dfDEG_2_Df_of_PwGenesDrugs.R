################################################################################
##  This program receives a data frame of your preference with pathways as rownames
##  a tabular definition of such pathways and a Data frame with DEG 
##  and retrieves a data frame of Pathway, Gene, Drugs
##  Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
################################################################################
##   Installing and loading the required libraries              ################
################################################################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("utils")) {
  install.packages("utils", ask =FALSE)
  library("utils")
}
source( paste0(Path_of_Results ,"DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R") )

################################################
# Data given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
Path_to_Ptwydb <- args[1] # The path to the information about the samples
# Path_to_Ptwydb <-c("../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv"); Path_to_your_Matrix <-c("../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt"); Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Top_Deregulted_Pathways/TCGA/") ; Label_for_Results <-"Df_pw_gn_TCGA"
Path_to_your_Matrix <- args[2] # The path to your matrix
# Path_to_your_Matrix <-c("TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt");
Path_of_Code <- args[3] # The path to your code
# Path_of_Code<-c("./")
Path_of_Results <-args[4] # # where do you want to save your results?
# Path_of_Results<-c("../Results/Top_Deregulted_Pathways/TCGA/")
Label_for_Results<-args[5] # Label for your results
# Label_for_Results <-"Df_pw_gn_TCGA"
## NOTE: for windows users #  If you use windows and you have troubles with your paths, please try out this tool: # Path_of_x <-choose.dir()
source()
##########################################################
### Reading the data #####################################
##########################################################
Df_your_Dereg_Pw <- read.table( Path_to_your_Matrix, sep= "\t", header=TRUE, row.names = 1, quote = "" )
PtwDf_defintions <- read.table( Path_to_Ptwydb, sep= "\t", header=TRUE, quote = "" , row.names = 1)

##########################################################
### Getting the genes inside your Deregulated pathways ####
##########################################################
Coincident_positions <- which( rownames( Df_your_Dereg_Pw ) %in% rownames(PtwDf_defintions) )
Df_yourPw_with_genes <- PtwDf_defintions[Coincident_positions, ]

##########################################################
###### Drugs 
##########################################################

##########################################################
###### Drugs 
##########################################################


list_yourPw_with_genes <- setNames(split(Df_yourPw_with_genes, seq(nrow(Df_yourPw_with_genes))), rownames(Df_yourPw_with_genes)) ## Making the df as a list

apply( , 1, )

