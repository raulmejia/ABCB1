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
Path_of_Results<-args[4] # # where do you want to save your results?
Labels <-args[4] # Label for your results
# Path_to_your_Matrix<-c("../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_METABRIC_.tsv_only_ABC_transporters_genes.tsv") ; Path_to_myPhenoData <- c("../Results/Lehmann-STROMA4/METABRIC/Lehmann\'s_Subt_and_properties_Numerical_METABRIC-only-Basals-lehmann_s_and_properties.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Clusterig/METABRIC/Heatmaps/Basal_Only_ABCS/") ; Labels <- "Heatmap_METABRIC_Basal_ABC_transporters_Lehmann"
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

########################
### Reading the data and creating the folder for the results
########################
dir.create(Path_of_Results , recursive = TRUE )
M.matrix<-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)

## Feature 
myPhenoData <-read.table(Path_to_myPhenoData ,header=TRUE,row.names = 1)
myPhenoData <- as.data.frame(myPhenoData )
My_pData_Anotdf <- AnnotatedDataFrame( data=myFeatureData)
validObject(My_fData_Anotdf)
##################################
### Building ExpressionSet Object
##################################
myExpressionSet <- ExpressionSet(as.matrix( M.matrix ), phenoData= My_pData_Anotdf  )

########################
### Pheatmap
########################

topGenes = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:500]
rowCenter = function(x) { x - rowMeans(x) }
?pheatmap
max(rowCenter(Biobase::exprs(x)[ topGenes, ]))
min(rowCenter(Biobase::exprs(x)[ topGenes, ]))
dim( rowCenter(Biobase::exprs(x)[ topGenes, ]))
mtcars %>%
  group_by(cyl) %>%
  summarise(mean = mean(disp), n = n())
?summarise
groups = group_by(pData(x), sampleGroup) %>%
  summarise(n = n(), color = unique(sampleColour))
groups
groupColor = setNames(groups$color, groups$sampleGroup)
vignette("programming")
?pheatmap
pheatmap( rowCenter(Biobase::exprs(x)[ topGenes, ] ),
          show_rownames = FALSE, show_colnames = FALSE,
          breaks = seq(-5, +5, length = 101),
          annotation_col =
            pData(x)[, c("sampleGroup", "Embryonic.day", "ScanDate") ],
          annotation_colors = list(
            sampleGroup = groupColor,
            genotype = c(`FGF4-KO` = "chocolate1", `WT` = "azure2"),
            Embryonic.day = setNames(brewer.pal(9, "Blues")[c(3, 6, 9)],
                                     c("E3.25", "E3.5", "E4.5")),
            ScanDate = setNames(brewer.pal(nlevels(x$ScanDate), "YlGn"),
                                levels(x$ScanDate))
          ),
          cutree_rows = 4
)

dir.create(Path_of_Results , recursive = TRUE )
M.matrix<-read.table(Path_to_your_Matrix,header=TRUE,row.names = 1)
