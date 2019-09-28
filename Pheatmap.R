################################################################################
# 2. Filtering and preprocessing of data for Pathifier
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
################################################################################
###########################################
## Data given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
Path_to_your_Matrix<-args[1] # The path to your matrix
Path_to_myPhenoData <- args[2] # The path to your Feature data
Path_of_Code<-args[3] # The path to your code
Path_of_Results <-args[4] # # where do you want to save your results?
Labels <-args[5] # Label for your results
log2transformation <- args[6] # Do you want log2transformation for your expression matrix ?
# Path_to_your_Matrix<-c("../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_METABRIC_.tsv_only_ABC_transporters_genes.tsv") ; Path_to_myPhenoData <- c("../Results/Lehmann-STROMA4/METABRIC/Lehmann\'s_Subt_and_properties_Numerical_METABRIC-only-Basals-lehmann_s_and_properties.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Clusterig/METABRIC/Heatmaps/Basal_Only_ABCS/") ; Labels <- "Heatmap_METABRIC_Basal_ABC_transporters_Lehmann_clasification" ; log2transformation <- "log2-transformation-no" 
# Path_to_your_Matrix <- choose.files() # choose.dir()

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
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("Hiiragi2013")) {
  BiocManager::install("Hiiragi2013", dependencies = TRUE)
  library(Hiiragi2013)
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap", dependencies = TRUE)
  library(pheatmap)
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("som")) {
  install.packages("som", dependencies = TRUE)
  library(som)
}

########################
### Reading the data and creating the folder for the results
########################
dir.create(Path_of_Results , recursive = TRUE )
M.matrix<-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)

# transforming your matrix to log2
if( log2transformation == "log2-transformation-no"){
  print("Your expression matrix has NOT being log2 Transformed by this program")
  M.matrix[1:5,1:5]
}else{
  if( log2transformation == "log2-transformation-yes" ){
    print("Your expression matrix has being log2 Transformed by this program")
    M.matrix <- log(M.matrix+1)
    M.matrix[1:5,1:5]
  }else{
    print("ERROR: Please type one option 'log2 transformation no' or 'log2 transformation no' ")
    exit("no")
  }  
}

##################################
### Building ExpressionSet Object
##################################
## PhenoData 
myPhenoData <-read.table(Path_to_myPhenoData ,header=TRUE,row.names = 1)
myPhenoData <- as.data.frame(myPhenoData )
My_pData_Anotdf <- AnnotatedDataFrame( data= myPhenoData)
validObject( My_pData_Anotdf )
## ExpresisionSet
myExpressionSet <- ExpressionSet(as.matrix( M.matrix ), phenoData= My_pData_Anotdf  )

########################
### Pheatmap
########################
rowCenter = function(x) { x - rowMeans(x) } # Defining a function that only make the 
PDSmatrix_zscores <- som::normalize( Biobase::exprs( myExpressionSet ) , byrow = TRUE)
colnames(PDSmatrix_zscores) <- colnames( Biobase::exprs( myExpressionSet ))

mymin <- min(exprs( myExpressionSet ))
mymax <- max(exprs( myExpressionSet ))
mydim <- dim(rowCenter(Biobase::exprs( myExpressionSet )))[2]

pdf( file = paste0(Path_of_Results,Labels,"_Z-scores_Euclidean_wardD2.pdf"))
pheatmap(  PDSmatrix_zscores,
           show_rownames = TRUE, show_colnames = FALSE,
           method = c("euclidean"),
           clustering_method = "ward.D2",
           fontsize = 7,
           main = gsub("_"," ",Labels),
           #breaks = seq(mymin, mymax, length = mydim+1 ),
           annotation_col =
             pData(myExpressionSet)[,c("D.stroma.property","B.stroma.property","T.stroma.property","E.stroma.property","MSL.property","M.property","LAR.property","IM.property","BL1.property","BL2.property")],
           annotation_colors = list(
             #D.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #B.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #T.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #E.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #MSL.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #M.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #LAR.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #IM.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL1.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL2.property = c(`1` = "black", `0` = "gray", `-1` = "white")
           ),
           cutree_rows = 4,
           cutree_cols = 4
)
dev.off()


pdf( file = paste0(Path_of_Results,Labels,"_Row_Center_Euclidean_wardD2.pdf") )
#png(file = paste0(Path_of_Results,Labels,"_Row_Center_Euclidean_wardD2.png"),width = 480, height = 480, pointsize = 2)
pheatmap(  rowCenter(Biobase::exprs( myExpressionSet )) ,
           main = gsub("_"," ",Labels),
           show_rownames = TRUE, show_colnames = FALSE,
           method = c("euclidean"),
           clustering_method = "ward.D2",
           fontsize = 7,
           #breaks = seq(mymin, mymax, length = mydim+1 ),
           annotation_col =
             pData(myExpressionSet)[,c("D.stroma.property","B.stroma.property","T.stroma.property","E.stroma.property","MSL.property","M.property","LAR.property","IM.property","BL1.property","BL2.property")],
           annotation_colors = list(
             #D.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #B.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #T.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #E.stroma.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #MSL.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #M.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #LAR.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #IM.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL1.property = c(`1` = "black", `0` = "gray", `-1` = "white"),
             #BL2.property = c(`1` = "black", `0` = "gray", `-1` = "white")
           ),
           cutree_rows = 4,
           cutree_cols = 4
)
dev.off()
