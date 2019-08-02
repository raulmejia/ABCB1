
Split_by_thisgene_level<- function( mymatrix, gene){
  mymatrix_25th_top_high <- mymatrix[ , mymatrix[gene,] > quantile(as.matrix(mymatrix[gene,]) , .75) ]
  mymatrix_high <- mymatrix[, mymatrix[gene,] > mean(as.matrix(mymatrix[gene,]))]
  mymatrix_low <- mymatrix[,mymatrix[gene,] < mean(as.matrix( mymatrix[gene,]))]
  mymatrix_25th_top_low <- mymatrix[ , mymatrix[gene,] < quantile(as.matrix(mymatrix[gene,]) , .25) ]
  result<- list(mymatrix_25th_top_high, mymatrix_high, mymatrix_low, mymatrix_25th_top_low)
  names(result) <- c("25th_top_high","high","low","25th_top_low")
  return(result)
}

merge_with_healthem <- function(list_splited, healthem){
  #list_splited <- list_of_splited_TCGAem_Tumour_by_ABCB1_levels
  #healthem <- TCGAem_Health
  list_splited[["25th_top_high"]] <- cbind(list_splited[["25th_top_high"]] , healthem )
  list_splited[["high"]] <- cbind(list_splited[["high"]] , healthem )
  list_splited[["25th_top_low"]] <- cbind(list_splited[["25th_top_low"]] , healthem )
  list_splited[["low"]] <- cbind(list_splited[["low"]] , healthem )
  return(list_splited)
}

splitdf_by_gene_level_in_tumours <- function(df, healthy_positions, tumour_positions,mygene){
  #This function split a data frame in health and disease, and create sub data frames of the
  # disease for example the 4th df has eliminated all the patients (columns) if 
  # their level of the choosed gene is below than the mean, below than 25 th percentil
  # then that submatrix will be pasted to the Healthy df
  # the function returns a list whose elements are each data frame  according that condition
  
  df_Health <- df[,healthy_positions]
  df_Tumour <- df[,tumour_positions]
  list_of_splited_df_Tumour_by_thisgene_levels <- Split_by_thisgene_level(df_Tumour,mygene)
  list_of_splited_df_by_thisgene_levels <- merge_with_healthem(list_of_splited_df_Tumour_by_thisgene_levels, df_Health)
}
