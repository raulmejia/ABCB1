########################################################################
### This script Get the Differential expresed genes in RNAseq data #####
# Writed by Raúl Alejandro Mejía Pedroza                              ##
########################################################################
#######   Data selected by te user #####################################
########################################################################
exp_mat_path<-c("../Results/Basal/25th_top_high.tsv")
#exp_mat_path <-c("../Results/25th_top_high.tsv")
Labels_path<-c("../Data/Labels_Controls_and_Normal_separated_TCGA.txt")
results_path<-c("../Results/ToyDeSeq/")
Label_for_results <- c("_DGE_TCGA")

mylfctreshold<-1
myp.adj<-0.01
mycoresavaiable<-7
# falta el nombre de las cosas a guardar
# Nota la categoria de referencia se llama "Control" 
################################################
######  Libraries needed    ####################
################################################
if (!require("BiocManager")) {
  biocLite("BiocManager")
  library(BiocManager)}
if (!require("DESeq2")) {
  BiocManager::install("DESeq2")
  library(DESeq2)}
if (!require("BiocParallel")) {
  BiocManager::install("BiocParallel")
  library(BiocParallel)}
if (!require("DEFormats")) {
  BiocManager::install("DEFormats")
  library(DEFormats)}

dir.create(results_path,recursive = TRUE)
########################################################################
#######   Loading the data #############################################
########################################################################
Exp_Mat <- read.table(exp_mat_path ,header=TRUE , sep= "\t")
Exp_Mat <- as.matrix(Exp_Mat)
MyLabels <-read.table(Labels_path)
########################################################################
#######   Extracting the propper labels ################################
########################################################################
Adjusted_Labels <- data.frame(MyLabels[which( rownames(MyLabels) %in% colnames(Exp_Mat) ),])
rownames(Adjusted_Labels) <- rownames(MyLabels)[ which( rownames(MyLabels) %in% colnames(Exp_Mat) )]
sum(as.vector(table(Adjusted_Labels)))- 112

#######################################################################
############   Building the DDS DESeqDataSet      #####################
#######################################################################
colnames(Adjusted_Labels)[1] <- "condition"
MyLabels$condition <- relevel(Adjusted_Labels$condition, ref="Control")
My_dds <- DESeqDataSetFromMatrix(countData = round(Exp_Mat), colData = Adjusted_Labels, design = ~ condition) 

#######################################################################
############   VST of your matrix               #####################
#######################################################################

My_vsd <- varianceStabilizingTransformation(My_dds)
saveRDS(My_vsd,file=paste0(results_path,"vsd_matrix_from",Label_for_results,"_.RDS"))
dists <- dist(t(assay(My_vsd)))
plot(hclust(dists))

#######################################################################
############   Now the DGE function               #####################
#######################################################################
time1<-proc.time()

My_dds_from_vsd_data <- DESeqDataSetFromMatrix(countData = round(assay(My_vsd)), colData = Adjusted_Labels, design = ~ condition) 
  register(MulticoreParam(mycoresavaiable))
  My_DESeq <- DESeq( My_dds_from_vsd_data,parallel = TRUE)
  register(MulticoreParam(mycoresavaiable))
  results_DESeq <-results( My_DESeq,parallel = TRUE)
time2 = proc.time()-time1


# Saving the results
save(My_DESeq,file=paste0(results_path,"DESeq_object_of_",Label_for_results,".RData"))
save(results_DESeq,file=paste0(results_path,"DESeq_results",Label_for_results,".RData"))

register(MulticoreParam(mycoresavaiable))
lfc1_results_DESeq <-results(My_DESeq,parallel = TRUE, lfcThreshold= mylfctreshold)

time3 = proc.time()-time1

# Saving the results and timming
save(lfc1_results_DESeq_list,file=paste0(results_path,"lfc1_results_DESeq",Label_for_results,".RData"))
write.table(c(time2,time3),file=paste0(results_path,Label_for_results,"Timing_runfunctions.txt"))


  pminus10_3<-grepl("TRUE",lfc1_results_DESeq$padj < myp.adj)
  padj10_3_lfc1_results_DESeq <-lfc1_results_DESeq[pminus10_3,]

#saving the results
save(padj10_3_lfc1_results_DESeq,file=paste0(results_path,"padj10_3_lfc1_results_DESeq",Label_for_results,".RData"))
write.table(padj10_3_lfc1_results_DESeq, 
            file=paste0(results_path,"padj10_3_lfc1_results_DESeq",Label_for_results,".tsv"),
            sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE
            )
?write.table






#######################################################################
############   Cutting the matrices               #####################
#######################################################################
#  
List_of_matrices<-list()
List_of_dfs <-list()
categories <-unique(as.character(AllLabels$Labels))
categories <-setdiff(categories,"Control")

#let's cut the initial matrix and dataframe
for(k in (categories)){
  # Making the vector to subsetting the data.
  selecteditems<-(grepl("Control",as.character(AllLabels$Labels)) | grepl(k,as.character(AllLabels$Labels)))
  # Subsetting the dataframe "ALLLabels"
  List_of_dfs[[k]]<-data.frame(as.character(AllLabels[selecteditems,]))
  rownames(List_of_dfs[[k]])<-rownames(AllLabels)[selecteditems]
  colnames(List_of_dfs[[k]])<-"condition"
  # # Subsetting the matrix "Exp_Mat"
  List_of_matrices[[k]]<-Exp_Mat[,rownames(AllLabels)[selecteditems]]
}

lapply(List_of_matrices,dim)
table(AllLabels)+112
#######################################################################
############    Building the dds objects          #####################
#######################################################################
dds_list<-list()
for(w in (categories)){
  List_of_dfs[[w]]$condition <- relevel(List_of_dfs[[w]]$condition, ref="Control")
  dds_list[[w]]<- DESeqDataSetFromMatrix(countData = round(List_of_matrices[[w]]),
                                colData = List_of_dfs[[w]],
                                design = ~ condition)  
}

#######################################################################
############   Now the DGE function               #####################
#######################################################################

library("BiocParallel")

time1<-proc.time()
DESeq_list<-list()
results_DESeq_list<-list()
for(x in categories){
  register(MulticoreParam(mycoresavaiable))
  DESeq_list[[x]]<- DESeq(dds_list[[x]],parallel = TRUE)
  register(MulticoreParam(mycoresavaiable))
  results_DESeq_list[[x]]<-results(DESeq_list[[x]],parallel = TRUE)
}
time2 = proc.time()-time1

# Saving the results
save(DESeq_list,file=paste0(results_path,"DESeq_list_TCGA-Mike.RData"))
save(results_DESeq_list,file=paste0(results_path,"results_DESeq_list_TCGA-Mike.RData"))

lfc1_results_DESeq_list<-list()
for(x in categories){
  register(MulticoreParam(mycoresavaiable))
  lfc1_results_DESeq_list[[x]]<-results(DESeq_list[[x]],parallel = TRUE, lfcThreshold= mylfctreshold)
}
time3 = proc.time()-time1

# Saving the results and timming
save(lfc1_results_DESeq_list,file=paste0(results_path,"lfc1_results_DESeq_list_TCGA-Mike.RData"))
write.table(c(time2,time3),file="Timing_runfunctions.txt")


padj10_3_lfc1_results_DESeq_list<-list()
for(k in categories){
  pminus10_3<-grepl("TRUE",lfc1_results_DESeq_list[[k]]$padj < myp.adj)
  padj10_3_lfc1_results_DESeq_list[[k]]<-lfc1_results_DESeq_list[[k]][pminus10_3,]
}
#saving the results
save(padj10_3_lfc1_results_DESeq_list,file=paste0(results_path,"padj10_3_lfc1_results_DESeq_list_TCGA-Mike.RData"))


