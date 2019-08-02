if (!require("dplyr")) {
  install.packages("dplyr", ask =FALSE)
  library(dplyr)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", ask =FALSE)
  library(ggplot2)
}
if (!require("limma")) {
  install.packages("limma", ask =FALSE)
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


#args = commandArgs(trailingOnly=TRUE)
#df_path <-args[1]

TCGAem_file <-c("/media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/DGE/Data/expMatrix_TCGA_cBioPortal_no_males.txt")
TCGAem <-read.table(file=TCGAem_file, sep="\t", quote = "", stringsAsFactors = FALSE)
TCGAem[1:5,1:4]
TCGAem <- log(TCGAem + 1)


ABCem <-TCGAem[grep("ABC",rownames(TCGAem)),]
rownames(ABCem)


sum(str_detect( colnames(ABCem) , "\\.01"))
colnames(ABCem)[566]
sum(str_detect( colnames(ABCem) , "\\.11"))


colnames(ABCem)
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
p1 <- ggplot(meltedABC[1:(10*920),], aes(Expression_Value, fill = gene_and_condition)) + geom_density(alpha = 0.2)
p1
p2 <- ggplot(meltedABC[((10*920)+1):(20*920),], aes(Expression_Value, fill=gene_and_condition)) + geom_density(alpha = 0.2)
p2

ggplot(meltedABC[((10*920)+1):(20*920),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())

ggplot(meltedABC[((0*920)+1):(25*920),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())

pdf("/media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/ABCB1_boxplot_heallthvsTumours.pdf", height = 7 , width = 7)
ggplot(meltedABC[((14*920)+1):(15*920),], aes(x = Gene, y = Expression_Value, fill = health)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge())
dev.off()

pdf("/media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/Results/boxplots_ABCs.pdf", height = 7 , width = 14)
ggplot(data = meltedABC[((0*920)+1):(20*920),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
ggplot(data = meltedABC[((21*920)+1):(40*920),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
ggplot(data = meltedABC[((41*920)+1):(50*920),], aes(Gene, Expression_Value)) +
  geom_boxplot(aes(colour = health))
dev.off()


pB1 <- ggplot(meltedABC[((14*920)+1):(15*920),], aes(Expression_Value, fill = gene_and_condition)) + geom_density(alpha = 0.2)
pB1
pB1 <- ggplot(meltedABC[((14*920)+1):(15*920),], aes(y=Expression_Value, fill = gene_and_condition)) + geom_boxplot(alpha = 0.2)
pB1
pB1 <- ggplot(meltedABC[((14*920)+1):(15*920),], aes(y=Expression_Value, x = gene_and_condition)) + geom_violin(alpha = 0.2)
pB1
pB1 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1, binwidth = 0.1)

p3 <- ggplot(meltedABC[((20*920)+1):(30*920),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p3
p4 <- ggplot(meltedABC[((30*920)+1):(40*920),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p4
p4 <- ggplot(meltedABC[((32*920)+1):(34*920),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p4
# Interesante ABCC6 y ABCC6P1 tienen distribuciones bimodales, que pacientes seran?
p5 <- ggplot(meltedABC[((49*920)+1):(50*920),], aes(Expression_Value, fill = Gene)) + geom_density(alpha = 0.2)
p5
# La distribucion de ABCG8 tambien esta curiosa (el 50)
# p5 nada solo unas de piquitos finales
#p4 interesante dos bimodales y una con un piquito al final
# p3 con una de un piquito final
# p2 contiene un ABCB algo muy interesante
# p1 falta análisis

##### Extracting the more 

source("/media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/ABCB1/ABCB1_code/splitdf_by_gene_level_in_tumours.R")
# let's split our data frame acoording to the levels of ABCB1, 
# we will build subdata frames, that only keep tumours with a ABCB1 expression level above the mean , under the mean , and so on...
list_TCGA_splited_by_ABCB1_levels_in_tumours <- splitdf_by_gene_level_in_tumours(TCGAem, 1:112, 113:920, "ABCB1")

lapply(list_TCGA_splited_by_ABCB1_levels_in_tumours,dim) # let's check
list_TCGA_splited_by_ABCB1_levels_in_tumours[[1]][1:5,1:5]

## Now pathifier


