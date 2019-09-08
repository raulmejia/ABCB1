### Pharmacological data bases
## Loading required libraries
if (!require("rDGIdb")) {
  BiocManager::install("rDGIdb", ask =FALSE)
  library(rDGIdb)
}

################################
## Data given by the user ######
################################
args <- commandArgs(trailingOnly = TRUE)
gene_path <- args[1]
#gene_path <- c("../Data/Datatoy_some_genes.tsv")
results_path <- args[2]
# results_path <-c("../Results/ToyDrugs/")
Label <- args[3]
# Label <- "unosgenes"

mygenes <-read.table( file=gene_path , sep="\t" , quote = "" , stringsAsFactors = FALSE )
mygenes <- mygenes[,1]
myquery <- queryDGIdb(mygenes)
#  resultSummary(myquery)[1:5,]
#  detailedResults(myquery)[,c(1,3,4)]

dir.create(paste0(results_path,"Drugs/"))
saveRDS(detailedResults(myquery)[,c(1,3,4)],file=paste0(results_path,"Drugs/",Label,"_gen_drug_df_lists.RDS"))

#pw_gen_drug_One_BigDf <- dflists_to_abigdf(pw_gen_drug_df)
write.table(detailedResults(myquery)[,c(1,3,4)], file=paste0(results_path,"Drugs/",Label,"_gen_drugDf.tsv") ,sep="\t",col.names = TRUE,row.names=FALSE,quote = FALSE)

# Other examples
#genes <- c("TNF", "AP1", "AP2", "XYZA")
#result <- queryDGIdb(genes)
#resultSummary(result)
#detailedResults(result)
#byGene(result)
#searchTermSummary(result)
#geneCategories()
#### reading the data