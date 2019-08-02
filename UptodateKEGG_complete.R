source("http://bioconductor.org/biocLite.R")
if (!require("gage")) {
  biocLite("gage", ask =FALSE)
  library(gage)
}
kegg_gsets<-kegg.gsets(species = "hsa", id.type = "kegg")
dir.create("../Results")
save(kegg_gsets,file="../Results/uptdateKeggPathways.RData")

if (!require("KEGGREST")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("KEGGREST", ask =FALSE)
  library("KEGGREST")
}
if (!require("stringr")) {
  install.packages("stringr", ask =FALSE)
  library("stringr")
}
results_path<-c("./")

list_of_chr_to_df<-function(X){
  # function that converts a list of character vectors to a data frame filled up with NA in blank spaces
  # X is a list of 310 lists, each "sublist" is a a list of 1 element, that element is a character vector
  mymax<-max(unlist(lapply(X,length)))
  Num_pathways<-length(lapply(X,length))
  mylist<-list()
    for(k in 1:Num_pathways){
      N_faltantes<-mymax-length(X[[k]])
      mylist[[k]]<-c(X[[k]],rep("NA",N_faltantes))
    }
      names(mylist)<-names(X)
      dd<-data.frame(matrix(nrow=Num_pathways,ncol=mymax))
      for(i in 1:length(mylist)){
        dd[i,] <- mylist[[i]]
      }
      rownames(dd)<-names(mylist)
  return(dd)
}

kegg_gsets <-kegg.gsets(species = "hsa", id.type = "kegg") # Downloading the most recent keggdb
kegg_gsets$kg.sets

KEGG_pathways_in_df <-list_of_chr_to_df(kegg_gsets$kg.sets)

Kegghsa <- keggList("hsa")
GeneSymbol_kid <- str_extract(Kegghsa, "[:alnum:]+")
names(GeneSymbol_kid) <- gsub("hsa:","",names(Kegghsa))
str(kegg_gsets$kg.sets)
kegg_sets_kid_gs <- kegg_gsets$kg.sets
kegg_gsets$kg.sets

for( w in 1:length(kegg_gsets$kg.sets) ){
  kegg_sets_kid_gs[[ names(kegg_gsets$kg.sets)[w] ]] <- GeneSymbol_kid[which(  names(GeneSymbol_kid) %in%  kegg_gsets$kg.sets[[ names(kegg_gsets$kg.sets)[w]   ]])]
}
kegg_sets_kid_gs <- lapply(kegg_sets_kid_gs,unique)

KEGG_pathways_in_df_genesymbols<-list_of_chr_to_df( kegg_sets_kid_gs)

write.table(KEGG_pathways_in_df_genesymbols, file=paste0(results_path,"KEGG_pathways_in_df_genesymbol.tsv") ,sep="\t",col.names = FALSE,quote = FALSE)
write.table(KEGG_pathways_in_df, file=paste0(results_path,"KEGG_pathways_in_df_in_keggids.tsv") ,sep="\t",col.names = FALSE,quote = FALSE)

save(KEGG_pathways_in_df_genesymbols,file=paste0(results_path,"KEGG_pathways_in_df_in_genesymbols.RData"))
save(KEGG_pathways_in_df,file=paste0(results_path,"KEGG_pathways_in_df.RData"))




