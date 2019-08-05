if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("dplyr")) {
  install.packages("dplyr", ask =FALSE)
  library(dplyr)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", ask =FALSE)
  library(ggplot2)
}
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library(limma)
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
if (!require("gage")) {
  BiocManager::install("gage", ask =FALSE)
  library("gage")
}
if (!require("KEGGREST")) {
  BiocManager::install("KEGGREST", ask =FALSE)
  library("KEGGREST")
}
if (!require("stringr")) {
  install.packages("stringr", ask =FALSE)
  library("stringr")
}

#Labels_Controls_and_Normal_separated_TCGA.txt


setwd("~/Documentos/ABCB1/ABCB1_code/")
dir.create("../Results/Basal")
results_path <-c("../Results/Basal/")
TCGAem_file_path <-c("../Data/Control_and_ Basal_with_indicator.txt")
#args = commandArgs(trailingOnly=TRUE)
#df_path <-args[1]

TCGAem_plus_indicator <-read.table(file=TCGAem_file_path, sep="\t", quote = "", stringsAsFactors = FALSE)
Indicator_Row <-TCGAem_plus_indicator[1,]
TCGAem <- TCGAem_plus_indicator[-1,]
TCGAem <- log(TCGAem + 1)
ABCem <-TCGAem[grep("ABC",rownames(TCGAem)),]
#ABCem <-TCGAem[grep("^TAP2",rownames(TCGAem)),]
#ABCem <-TCGAem[grep("^TAP1",rownames(TCGAem)),]
#ABCem[1:22,1:5]
TCGAem <- rbind(Indicator_Row,TCGAem)

sum(str_detect( colnames(ABCem) , "\\.01")); colnames(ABCem)[566] ; sum(str_detect( colnames(ABCem) , "\\.11"))

meltedABC <- melt(t(ABCem))
colnames(meltedABC)[c(1,2,3)]<-c("Sample","Gene","Expression_Value")

head(meltedABC,n=5)
meltedABC[910:930,]
dim(meltedABC)
meltedABC$health <- as.character(str_detect(meltedABC$Sample, "\\.11"))
meltedABC$id_and_condition <- paste0(meltedABC$Sample , meltedABC$health)
meltedABC$gene_and_condition <- paste0(meltedABC$Gene , meltedABC$health)
#meltedABClog2[,"Expression_Value"]<- log2(meltedABClog2[,"Expression_Value"] +1)

## looking for the treshold
p1 <- ggplot(meltedABC[1:(10*248),], aes(Expression_Value, fill = gene_and_condition)) + geom_density(alpha = 0.2)
p1
p1 <- ggplot(meltedABC[1:(248),], aes(Expression_Value, fill = gene_and_condition)) + geom_density(alpha = 0.2)
p1
p2 <- ggplot(meltedABC[((10*248)+1):(20*248),], aes(Expression_Value, fill=gene_and_condition)) + geom_density(alpha = 0.2)
p2

ggplot(meltedABC[((10*248)+1):(20*248),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())

ggplot(meltedABC[((0*248)+1):(25*248),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())

pdf(paste0(results_path,"ABCB1_boxplot_heallthvsTumours.pdf"), height = 7 , width = 7)
ggplot(meltedABC[((14*248)+1):(15*248),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())
dev.off()

pdf(paste0(results_path,"ABCB1_boxplot_heallthvsTumours.pdf"), height = 7 , width = 14)
ggplot(data = meltedABC[((0*248)+1):(20*248),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
ggplot(data = meltedABC[((21*248)+1):(40*248),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
ggplot(data = meltedABC[((41*248)+1):(50*248),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
dev.off()

pB1 <- ggplot(meltedABC[((14*248)+1):(15*248),], aes(Expression_Value, fill = gene_and_condition)) + geom_density(alpha = 0.2)
pB1
pB1 <- ggplot(meltedABC[((14*248)+1):(15*248),], aes(y=Expression_Value, fill = gene_and_condition)) + geom_boxplot(alpha = 0.2)
pB1
pB1 <- ggplot(meltedABC[((14*248)+1):(15*248),], aes(y=Expression_Value, x = gene_and_condition)) + geom_violin(alpha = 0.2)
pB1
pB1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1, binwidth = 0.1)

p3 <- ggplot(meltedABC[((20*248)+1):(30*248),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p3
p4 <- ggplot(meltedABC[((30*248)+1):(40*248),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p4
p4 <- ggplot(meltedABC[((32*248)+1):(34*248),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p4
# Interesante ABCC6 y ABCC6P1 tienen distribuciones bimodales, que pacientes seran?
p5 <- ggplot(meltedABC[((49*248)+1):(50*248),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p5
# More clear binomial distributions between healthies and tumours
# End Visualization

##### Spliting by one gene expression level

source("splitdf_by_gene_level_in_tumours.R")
source("ConvertKEGGList_to_DataFrame.R")
# let's split our data frame acoording to the levels of ABCB1, 
# we will build subdata frames, that only keep tumours with a ABCB1 expression level above the mean , under the mean , and so on...

list_TCGA_splited_by_ABCB1_levels_in_tumours <- splitdf_by_gene_level_in_tumours(TCGAem, 1:112, 113:248, "ABCB1")

lapply(list_TCGA_splited_by_ABCB1_levels_in_tumours,dim) # let's check

list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators <- list_TCGA_splited_by_ABCB1_levels_in_tumours
list_TCGA_splited_by_ABCB1_levels_in_tumours[[1]][1:5,1:5]
list_TCGA_splited_by_ABCB1_levels_in_tumours[[2]][1:5,1:5]
list_TCGA_splited_by_ABCB1_levels_in_tumours[[3]][1:5,1:5]
list_TCGA_splited_by_ABCB1_levels_in_tumours[[4]][1:5,1:5]

#for(  k in 1:length(list_TCGA_splited_by_ABCB1_levels_in_tumours) ){
  #indicator_row<-c(rep(1,112),rep(0,dim(list_TCGA_splited_by_ABCB1_levels_in_tumours[[k]])[2] - 112))
  #list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators[[k]]<-rbind(indicator_row, list_TCGA_splited_by_ABCB1_levels_in_tumours[[k]])
  #list_TCGA_splited_by_ABCB1_levels_in_tumours[[1]][1:5,1:5]
  #rownames(list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators[[k]])[1]<- "NORMAL"
#}

for(  k in 1: length(list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators)){
  write.table(list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators[[k]]  ,
              file=paste0(results_path,names(list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators)[k],".tsv" ),
              sep="\t", quote=FALSE , row.names= TRUE, col.names= TRUE  )    
}

## Now pathifier
#ext 1152 Jorge alvarado

# First updating KEGG
kegg_gsets <-kegg.gsets(species = "hsa", id.type = "kegg") # Downloading the most recent keggdb
source("ConvertKEGGList_to_DataFrame.R")
KEGG_pathways_in_df <-list_of_chr_to_df(kegg_gsets$kg.sets)

# Translating the kgid to genesymbols
Kegghsa <- keggList("hsa")
GeneSymbol_kid <- str_extract(Kegghsa, "[:alnum:]+")
names(GeneSymbol_kid) <- gsub("hsa:","",names(Kegghsa))
kegg_sets_kid_gs <- kegg_gsets$kg.sets
for( w in 1:length(kegg_gsets$kg.sets) ){
  kegg_sets_kid_gs[[ names(kegg_gsets$kg.sets)[w] ]] <- GeneSymbol_kid[which(  names(GeneSymbol_kid) %in%  kegg_gsets$kg.sets[[ names(kegg_gsets$kg.sets)[w]   ]])]
}
kegg_sets_kid_gs <- lapply(kegg_sets_kid_gs,unique)
KEGG_pathways_in_df_genesymbols<-list_of_chr_to_df( kegg_sets_kid_gs)

write.table(KEGG_pathways_in_df_genesymbols, file=paste0(results_path,"KEGG_pathways_in_df_genesymbol.tsv") ,sep="\t",col.names = FALSE,quote = FALSE)
write.table(KEGG_pathways_in_df, file=paste0(results_path,"KEGG_pathways_in_df_in_keggids.tsv") ,sep="\t",col.names = FALSE,quote = FALSE)
save(KEGG_pathways_in_df_genesymbols,file=paste0(results_path,"KEGG_pathways_in_df_in_genesymbols.RData"))
save(KEGG_pathways_in_df,file=paste0(results_path,"KEGG_pathways_in_df.RData"))

#KEGG_pathways_in_df_genesymbols[1:5,1:5]
#rownames(KEGG_pathways_in_df_genesymbols)[grep("hsa01522", rownames(KEGG_pathways_in_df_genesymbols))]
#KEGG_pathways_in_df_genesymbols[grep("hsa01522", rownames(KEGG_pathways_in_df_genesymbols)),]
#str(KEGG_pathways_in_df_genesymbols)
#intersect(as.character(KEGG_pathways_in_df_genesymbols[grep("hsa01522", rownames(KEGG_pathways_in_df_genesymbols)),1:97]),
#as.character(KEGG_pathways_in_df_genesymbols[grep("hsa01524", rownames(KEGG_pathways_in_df_genesymbols)),1:73]))

#length(as.character(KEGG_pathways_in_df_genesymbols[grep("hsa01524", rownames(KEGG_pathways_in_df_genesymbols)),1:73]))
#length(as.character(KEGG_pathways_in_df_genesymbols[grep("hsa01522", rownames(KEGG_pathways_in_df_genesymbols)),1:97]))
#length(intersect(as.character(KEGG_pathways_in_df_genesymbols[grep("hsa01522", rownames(KEGG_pathways_in_df_genesymbols)),1:97]),
#          as.character(KEGG_pathways_in_df_genesymbols[grep("hsa01524", rownames(KEGG_pathways_in_df_genesymbols)),1:73])))


# go to script
dir.create( paste0(results_path,"Pathifier/"))
# Falta correr estos globales de nuevo porque actualice el script y ahora estoy usando los procesadores con el basal , queda de tare
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_10" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_low_ABCB1_stbl_10" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/high.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_high_ABCB1_stbl_10" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/25th_top_high.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_above_75_stbl_10" )

system("Rscript Pathifier_Args_3stabilizing_Filtervalue3_75.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_3" ) # lo más rápido
# Rscript Pathifier_Args_100stabilizing_Filtervalue3_75.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_100
# Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_10_Filter3
---
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Basal/25th_top_low.tsv ../Results/Basal/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Basal/Pathifier/ Basal_under_25_stbl_10_threshold_3" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Basal/low.tsv ../Results/Basal/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Basal/Pathifier/ Basal_low_ABCB1_stbl_10_threshold_3" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Basal/high.tsv ../Results/Basal/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Basal/Pathifier/ Basal_high__ABCB1_stbl_10_threshold_3" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Basal/25th_top_high.tsv ../Results/Basal/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Basal/Pathifier/ Basal_above_75_stbl_10_threshold_3" )

system("Rscript Pathifier_Args_100stabilizing_Filtervalue3.R ../Results/Basal/25th_top_low.tsv ../Results/Basal/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Basal/Pathifier/ Basal_under_25_stbl_100_threshold_3" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Basal/25th_top_low.tsv ../Results/Basal/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Basal/Pathifier/ Basal_under_25_stbl_10_threshold_3_75")




#### Pharmacological data bases
if (!require("rDGIdb")) {
  BiocManager::install("rDGIdb", ask =FALSE)
  library(rDGIdb)
}
genes <- c("TNF", "AP1", "AP2", "XYZA")
result <- queryDGIdb(genes)
str(result)
resultSummary(result)
detailedResults(result)
byGene(result)
searchTermSummary(result)
geneCategories()
####
PDSzpaths <-read.table( file="../Results/Basal/Pathifier/Basal_under_25_stbl_10_threshold_3_75_median_PDSz_ordered_matrix.txt" , sep="\t" , quote = "" , stringsAsFactors = FALSE )
head(PDSzpaths)
str(KEGG_pathways_in_df_genesymbols)

ordered_by_PDS <- KEGG_pathways_in_df_genesymbols[rownames(PDSzpaths),]
#adf<-ordered_by_PDS
hola<- queryDGIdb(as.character(ordered_by_PDS[1,]))
str(hola)
hola$interaction
resultSummary(hola)

#w <-1
pathwaydf2drugs <- function(adf){
  list_of_drugdfs <- list()
  for(w in 1:dim(adf)[1] ){
    pw_DGI<-  queryDGIdb(as.character(adf[w,]))
    gene_drug_df <- resultSummary(pw_DGI)[,1:2]
    list_of_drugdfs[[w]] <- data.frame(rep(rownames(adf)[w],dim(gene_drug_df)[1]), gene_drug_df, stringsAsFactors = FALSE)
    colnames(list_of_drugdfs[[w]])[1]<- "Pathway"
  }
  names(list_of_drugdfs) <- rownames(adf)
  return(list_of_drugdfs)
}
pw_gen_drug_df <-pathwaydf2drugs(ordered_by_PDS)

pathwaydf2drugs_detailedres <- function(adf){
  list_of_drugdfs <- list()
  for(w in 1:dim(adf)[1] ){
    pw_DGI<-  queryDGIdb(as.character(adf[w,]))
    gene_drug_df <- detailedResults(pw_DGI)[,2:4]
    list_of_drugdfs[[w]] <- data.frame(rep(rownames(adf)[w],dim(gene_drug_df)[1]), gene_drug_df, stringsAsFactors = FALSE)
    colnames(list_of_drugdfs[[w]])[1]<- "Pathway"
  }
  names(list_of_drugdfs) <- rownames(adf)
  return(list_of_drugdfs)
}
pw_gen_drug_df_list_detailed <- pathwaydf2drugs_detailedres(ordered_by_PDS)

#pw_gen_drug_df <-pathwaydf2drugs(ordered_by_PDS[1:5,])
#length(pw_gen_drug_df)
#names(pw_gen_drug_df)
dir.create(paste0(results_path,"Drugs/"))
saveRDS(pw_gen_drug_df,file=paste0(results_path,"Drugs/","pw_gen_drug_df_lists.RDS"))
saveRDS(pw_gen_drug_df_list_detailed ,file=paste0(results_path,"Drugs/","pw_gen_drug_df_list_detailed.RDS"))
dflists_to_abigdf <- function(mylist){
  #  for(){
  BigDf<-mylist[[1]]
  for( z in 2:length(mylist)){
    BigDf <- rbind( BigDf , mylist[[z]])
  }
  return(BigDf)
}

pw_gen_drug_One_BigDf <- dflists_to_abigdf(pw_gen_drug_df)
write.table(pw_gen_drug_One_BigDf, file=paste0(results_path,"Drugs/","pw_gen_drug_One_BigDf.tsv") ,sep="\t",col.names = TRUE,row.names=FALSE,quote = FALSE)
pw_gen_drug_One_BigDf_including_interactions <- dflists_to_abigdf(pw_gen_drug_df_list_detailed)
write.table(pw_gen_drug_One_BigDf_including_interactions, file=paste0(results_path,"Drugs/","pw_gen_drug_One_BigDf_including_interactions.tsv") ,sep="\t",col.names = TRUE,row.names=FALSE,quote = FALSE)
head(pw_gen_drug_One_BigDf_including_interactions, 110)


head(pw_gen_drug_One_BigDf, n=100)
pw_gen_drug_One_BigDf[order(pw_gen_drug_One_BigDf[,2]),]
pw_gen_drug_df[[4]]
dim(pw_gen_drug_One_BigDf)


## differential expresion 
names(list_TCGA_splited_by_ABCB1_levels_in_tumours)
colnames(list_TCGA_splited_by_ABCB1_levels_in_tumours[[4]])
All_labels<-read.table("/media/raulmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/DGE/Data/Labels_Controls_and_Normal_separated_TCGA.txt")
Ctrl_basal_ABCB125thlow_labels <- data.frame(All_labels[which(rownames(All_labels) %in% colnames(list_TCGA_splited_by_ABCB1_levels_in_tumours[[4]])),])
rownames(Ctrl_basal_ABCB125thlow_labels) <-  rownames(All_labels)[which(rownames(All_labels) %in% colnames(list_TCGA_splited_by_ABCB1_levels_in_tumours[[4]]))]
colnames(Ctrl_basal_ABCB125thlow_labels)<- "Labels"
write.table(Ctrl_basal_ABCB125thlow_labels, file=paste0(results_path,"Ctrl_basal_ABCB125thlow_labels.tsv") ,sep="\t",col.names = TRUE)


####
pw_gen_drug_One_BigDf_including_interactions_onlyDEG <- pw_gen_drug_One_BigDf_including_interactions[which(pw_gen_drug_One_BigDf_including_interactions[,2] %in% rownames(padj10_3_lfc1_results_DESeq_list[[1]])),]
write.table(pw_gen_drug_One_BigDf_including_interactions_onlyDEG, file=paste0(results_path,"pw_gen_drug_One_BigDf_including_interactions_onlyDEG.tsv") ,sep="\t",col.names = TRUE)

# Lo mejorcito
# Jaccard 
if (!require("OmicsMarkeR")) {
  BiocManager::install("OmicsMarkeR", ask =FALSE)
  library(OmicsMarkeR)
}

for( k in 1:length(kegg_sets_kid_gs)){
  
}
as.character(kegg_sets_kid_gs[[1]])
jaccard()

("OmicsMarkeR"
#

Rscript /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/ABCB1_code/Pathifier_Args_100stabilizing_Filtervalue3_75.R /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/25th_top_low.tsv /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/KEGG_pathways_in_df_genesymbol.tsv /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/ABCB1_code/ /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/ All_subtypes_low_ABCB1
paste0(results_path,"low",".tsv" )
# 3.75