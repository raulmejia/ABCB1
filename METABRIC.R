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

dir.create(results_path <-c("../Results/METABRIC/"))
setwd("~/Documentos/ABCB1/ABCB1_code/")
results_path <-c("../Results/METABRIC/")
TCGAem_file_path <-c("/media/raulmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/Pathifier/METABRIC/Data/joined_indicator.txt")
#args = commandArgs(trailingOnly=TRUE)
#df_path <-args[1]

TCGAem <-read.table(file=TCGAem_file_path, sep="\t", quote = "", stringsAsFactors = FALSE)
TCGAem <- log(TCGAem + 1)

ABCem <-TCGAem[grep("ABC",rownames(TCGAem)),]

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
p1 <- ggplot(meltedABC[1:(10*1141),], aes(Expression_Value, fill = gene_and_condition)) + geom_density(alpha = 0.2)
p1
p2 <- ggplot(meltedABC[((10*1141)+1):(20*1141),], aes(Expression_Value, fill=gene_and_condition)) + geom_density(alpha = 0.2)
p2

ggplot(meltedABC[((10*1141)+1):(20*1141),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())

ggplot(meltedABC[((0*1141)+1):(25*1141),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())

pdf(paste0(results_path,"ABCB1_boxplot_heallthvsTumours.pdf"), height = 7 , width = 7)
ggplot(meltedABC[((14*1141)+1):(15*1141),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())
dev.off()


pdf(paste0(results_path,"ABCB1_boxplot_heallthvsTumours.pdf"), height = 7 , width = 14)
ggplot(data = meltedABC[((0*1141)+1):(20*1141),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
ggplot(data = meltedABC[((21*1141)+1):(40*1141),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
ggplot(data = meltedABC[((41*1141)+1):(50*1141),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
dev.off()


pB1 <- ggplot(meltedABC[((14*1141)+1):(15*1141),], aes(Expression_Value, fill = gene_and_condition)) + geom_density(alpha = 0.2)
pB1
pB1 <- ggplot(meltedABC[((14*1141)+1):(15*1141),], aes(y=Expression_Value, fill = gene_and_condition)) + geom_boxplot(alpha = 0.2)
pB1
pB1 <- ggplot(meltedABC[((14*1141)+1):(15*1141),], aes(y=Expression_Value, x = gene_and_condition)) + geom_violin(alpha = 0.2)
pB1
pB1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1, binwidth = 0.1)

p3 <- ggplot(meltedABC[((20*1141)+1):(30*1141),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p3
p4 <- ggplot(meltedABC[((30*1141)+1):(40*1141),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p4
p4 <- ggplot(meltedABC[((32*1141)+1):(34*1141),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p4
# Interesante ABCC6 y ABCC6P1 tienen distribuciones bimodales, que pacientes seran?
p5 <- ggplot(meltedABC[((49*1141)+1):(50*1141),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p5
# La distribucion de ABCG8 tambien esta curiosa (el 50)
# p5 nada solo unas de piquitos finales
#p4 interesante dos bimodales y una con un piquito al final
# p3 con una de un piquito final
# p2 contiene un ABCB algo muy interesante
# p1 falta análisis

##### Spliting by one gene expression level

source("splitdf_by_gene_level_in_tumours.R")
source("ConvertKEGGList_to_DataFrame.R")
# let's split our data frame acoording to the levels of ABCB1, 
# we will build subdata frames, that only keep tumours with a ABCB1 expression level above the mean , under the mean , and so on...
list_TCGA_splited_by_ABCB1_levels_in_tumours <- splitdf_by_gene_level_in_tumours(TCGAem, 1:112, 113:1141, "ABCB1")

names(list_TCGA_splited_by_ABCB1_levels_in_tumours)
lapply(list_TCGA_splited_by_ABCB1_levels_in_tumours,dim) # let's check
list_TCGA_splited_by_ABCB1_levels_in_tumours[[1]][1:5,1:5]

list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators <- list_TCGA_splited_by_ABCB1_levels_in_tumours
list_TCGA_splited_by_ABCB1_levels_in_tumours[[1]][1:5,1:5]
list_TCGA_splited_by_ABCB1_levels_in_tumours[[2]][1:5,1:5]
list_TCGA_splited_by_ABCB1_levels_in_tumours[[3]][1:5,1:5]
list_TCGA_splited_by_ABCB1_levels_in_tumours[[4]][1:5,1:5]

list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators[[1]][1:5,1:5]
for(  k in 1:length(list_TCGA_splited_by_ABCB1_levels_in_tumours) ){
  indicator_row<-c(rep(1,112),rep(0,dim(list_TCGA_splited_by_ABCB1_levels_in_tumours[[k]])[2] - 112))
  list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators[[k]]<-rbind(indicator_row, list_TCGA_splited_by_ABCB1_levels_in_tumours[[k]])
  #list_TCGA_splited_by_ABCB1_levels_in_tumours[[1]][1:5,1:5]
  rownames(list_TCGA_splited_by_ABCB1_levels_in_tumours_indicators[[k]])[1]<- "NORMAL"
}

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

paste0(results_path,"low",".tsv" )

# go to script
dir.create( paste0(results_path,"Pathifier/"))
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_10" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_low_ABCB1_stbl_10" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/high.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_high_ABCB1_stbl_10" )
system("Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/25th_top_high.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_above_75_stbl_10" )

system("Rscript Pathifier_Args_3stabilizing_Filtervalue3_75.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_3" ) # lo más rápido
# Rscript Pathifier_Args_100stabilizing_Filtervalue3_75.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_100
# Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/25th_top_low.tsv ../Results/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/ All_subtypes_under_25_stbl_10_Filter3


# Lo mejorcito


Rscript /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/ABCB1_code/Pathifier_Args_100stabilizing_Filtervalue3_75.R /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/25th_top_low.tsv /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/KEGG_pathways_in_df_genesymbol.tsv /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/ABCB1_code/ /media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/ All_subtypes_low_ABCB1
paste0(results_path,"low",".tsv" )


# 3.75
