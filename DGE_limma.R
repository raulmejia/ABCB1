################################################################################
## This program calculates the Differential expressed genes through limma
## Author: Raúl Alejandro Mejía Pedroza github: https://github.com/raulmejia
################################################################################
### Data given by the user
################################################################################
args <- commandArgs(trailingOnly = TRUE)
exp_mat_path <-args[1]
# exp_mat_path <-c("../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv") ; Labels_path <-c("../Data/Labels_Ctrl_and_NL_Recal_separated_METABRIC.txt") ; Path_of_Code<-c("./") ; Path_of_Results<-c("../Results/DEG/METABRIC/") ; mypvalue = 0.05 ; mylfc = 0.6 ; Label_for_results <-"METABRIC_some_LumAs_vs_Controls_"  
# exp_mat_path <- c("../Data/joined_indicator_METABRIC.txt")  
Labels_path <-args[2]
Path_of_Code <-args[3]
Path_of_Results <-args[4]
mylfc <-args[5]
mypvalue <- args[6]
Label_for_results <- args[7]
#Filter_value <- args[6] # Filter low value genes # METABRIC = 3

###############################################################################
### Installing and/or loading required packages
###############################################################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("illuminaHumanv3.db")) {
  BiocManager::install("illuminaHumanv3.db", ask =FALSE)
  library("illuminaHumanv3.db")
}
if (!require("hgu133plus2.db")) {
  BiocManager::install("hgu133plus2.db", ask =FALSE)
  library("hgu133plus2.db")
}
if (!require("annotate")) {
  BiocManager::install("annotate", ask =FALSE)
  library("annotate")
}
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
if (!require("gProfileR")) {
  BiocManager::install("gProfileR", ask =FALSE)
  library("gProfileR")
}
if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db", ask =FALSE)
  library("org.Hs.eg.db")}

##########################################
###### Reading the data 
##########################################
# Load the data frame (expression matrix)
Exp_Mat <- read.table(exp_mat_path ,header=TRUE , sep="\t")
# New_Exp_Mat <- cbind( Exp_Mat[,1:10] , Exp_Mat[, c(145,147,148, 150, 156, 158, 159, 160, 163, 164, 165, 166)])

# Load the annotation data to extract subtypes
Labels <-read.table(Labels_path)

############################################
### Preprocessing some data given y the user
############################################
# coercing to numeric 
#Filter_value <- as.numeric( Filter_value ) 
mypvalue <- as.numeric(  mypvalue)
mylfc <- as.numeric( mylfc)
dir.create(Path_of_Results, recursive = TRUE)
# Erasing the row of normals
normals_pos <- which(rownames(Exp_Mat) %in% "NORMALS")
Exp_Mat <- Exp_Mat[-normals_pos ,]

###########################################
####### Hadling the labels
###########################################
# convert into factors
samples <-Labels$Labels
# check factors have been assigned
samples
#######################
 positions <- which( rownames(Labels) %in% colnames(New_Exp_Mat) )
 NewLabels <- Labels[positions,]
 NewLabels <- droplevels(NewLabels)
 DfNewLabels <- data.frame(NewLabels)
 rownames( DfNewLabels ) <- colnames(New_Exp_Mat)
 head(DfNewLabels)
 samples <- DfNewLabels$NewLabels

##################### 
 if( length(levels(samples) )==2  ) {
   print("OK: you only have two levels")
 }else{ 
   print("ERROR: you number of labels is different than 2")
   stop("ERROR: you number of labels is different than 2")
 }
 
 if( length(which(levels(samples) %in% "Control"))!=0  ) {
   print("OK: you have a Label named -Control-")
 }else{ 
   print("ERROR: you DON'T have a label named -Control-")
   stop("ERROR: you DON'T have a label named -Control-")
 }
 
#####################   Changing the Label to "Case"
 new_samples <- as.character(samples)
 your_wild_label <- setdiff(levels(samples), "Control")
 new_samples <- gsub( your_wild_label, "Case", new_samples)
 new_samples <- factor(new_samples  , labels = c("Control", "Case") )

######################
 # set up the experimental design
 design <- model.matrix(~0 + samples)
 colnames(design) <- c("Control", "Case") # Como level el segundo 
 
 # fit the linear model to the your expression set "eset"
 eset<-Exp_Mat
 fit <- lmFit(eset, design)
 # set up a contrast matrix to compare tissues v cell line
 contrast.matrix <- makeContrasts( Case_Controls = Case - Control,  levels=design)
 # Now the contrast matrix is combined with the per-probeset linear model fit.
 huvec_fits <- contrasts.fit(fit, contrast.matrix)
 huvec_ebFit <- eBayes(huvec_fits)
 
#####################################
##DGEs with  Affy ids no anotation###
#####################################
 results_only_affys<-list(length(colnames(contrast.matrix)))
 for(i in 1:length(colnames(contrast.matrix))){
   DGEs<-topTable(huvec_ebFit, coef=colnames(contrast.matrix)[i], number=10000, p.value = mypvalue, lfc = mylfc)
   DGEs<-DGEs[order(rownames(DGEs)),]
   # Saving the matrix
   write.table(DGEs, paste(Path_of_Results ,"DGE_",Label_for_results, mypvalue, mylfc, colnames(contrast.matrix)[i], "_ID_with_ILMN.txt", sep=""), sep="\t", quote=FALSE)
   # Saving with no ILMN_ identifiers
   write.table(DGEs[!grepl("ILMN_",rownames(DGEs)),], paste(Path_of_Results ,"DGE_",Label_for_results, mypvalue, mylfc, colnames(contrast.matrix)[i], "_ID_no_ILMN_.txt", sep=""), sep="\t", quote=FALSE)
   results_only_affys[[i]]<-DGEs
 }
 saveRDS(results_only_affys,file= paste(Path_of_Results ,"DGE_",Label_for_results, mypvalue, mylfc, colnames(contrast.matrix)[i], "_ID_with_ILMN.RDS", sep=""))
 
