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
# Examples of use
#KEGG_pathways_in_df<-list_of_chr_to_df(kegg_gsets$kg.sets)
#write.table(KEGG_pathways_in_df,file="KEGG_pathways_in_df.tsv",sep="\t",col.names = FALSE,quote = FALSE)
#save(KEGG_pathways_in_df,file="KEGG_pathways_in_df.RData")
