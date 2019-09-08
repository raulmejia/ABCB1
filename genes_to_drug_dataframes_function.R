### Pharmacological data bases
## Loading required libraries
if (!require("rDGIdb")) {
  BiocManager::install("rDGIdb", ask =FALSE)
  library(rDGIdb)
}

################################
## Data given by the user ######
################################
AnottateDrugsForthisCharacter <- function(myvector){
  myquery <- queryDGIdb(mygenes)
  return(detailedResults(myquery)[,c(1,3,4)])
}

