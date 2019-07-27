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

args = commandArgs(trailingOnly=TRUE)
df_path <-args[1]

TCGAem_file <-c("/media/rmejia/ADATA/boba-bk-postsismo/rmejia/Documents/Doctorado/ProyectoDoctorado/Pipe_post_IncosistencyPatways/4_Enrichment/DGE/Data/expMatrix_TCGA_cBioPortal_no_males.txt")
TCGAem <-read.table(file=TCGAem_file, sep="\t", quote = "", stringsAsFactors = FALSE)
TCGAem[1:5,1:5]


ABCem <-TCGAem[grep("ABC",rownames(TCGAem)),]
rownames(ABCem)
ABCcolors<-c(rep("purple",14),rep("blue",9),rep("cyan",14),rep("green",4),rep("yellow",1),rep("pink",3),rep("red",5))

dim(ABCem)
meltedABC<- melt(t(ABCem))
head(meltedABC,n=5)
meltedABC[910:930,]
str(meltedABC)

meltedABClog2<-meltedABC
meltedABClog2[,"value"]<- log2(meltedABClog2[,"value"])+1
p1 <- ggplot(meltedABClog2[1:(10*920),], aes(value, fill = Var2)) + geom_density(alpha = 0.2)
p1
p2 <- ggplot(meltedABClog2[((10*920)+1):(20*920),], aes(value, fill = Var2)) + geom_density(alpha = 0.2)
p2
p3 <- ggplot(meltedABClog2[((20*920)+1):(30*920),], aes(value, fill = Var2)) + geom_density(alpha = 0.2)
p3
p4 <- ggplot(meltedABClog2[((30*920)+1):(40*920),], aes(value, fill = Var2)) + geom_density(alpha = 0.2)
p4
p5 <- ggplot(meltedABClog2[((40*920)+1):(50*920),], aes(value, fill = Var2)) + geom_density(alpha = 0.2)
p5
# p5 nada solo unas de piquitos finales
#p4 interesante dos bimodales y una con un piquito al final
# p3 con una de un piquito final
# p2 contiene un ABCB algo muy interesante
# p1 falta análisis