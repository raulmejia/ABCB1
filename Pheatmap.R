################################################################################
# 2. Filtering and preprocessing of data for Pathifier
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
################################################################################
###########################################
## Data given by the user
###########################################
args <- commandArgs(trailingOnly = TRUE)
Path_to_your_Matrix<-args[1] # The path to your matrix
# Path_to_your_Matrix<-c("../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_METABRIC_.tsv_only_ABC_transporters_genes.tsv") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/Clusterig/") ; Labels <- "METABRIC_Only_Basal_Only_ABC_transporters"
# Path_to_your_Matrix<-c("Control_with_subexpression_matrix_of_Basal_from_METABRIC_.tsv_only_ABC_transporters_genes.tsv")
# Path_to_your_Matrix <- choose.files() # choose.dir()
Path_of_Code<-args[2] # The path to your code
Path_of_Results<-args[3] # # where do you want to save your results?
# Path_of_Results <- choose.dir()
# choose.dir(default = "", caption = "Select folder")
Labels <-args[4] # Label for your results

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

########################
### Reading the data and creating the folder for the results
########################
library("Hiiragi2013")
data("x")
x
str(x)
dim(Biobase::exprs(x))
class(Biobase::exprs(x))
Biobase::exprs(x)[1:10,1:10]
attributes(x)

library("pheatmap")
topGenes = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:500]
rowCenter = function(x) { x - rowMeans(x) }
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
