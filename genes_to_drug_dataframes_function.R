### Pharmacological data bases
## Loading required libraries
if (!require("rDGIdb")) {
  BiocManager::install("rDGIdb", ask =FALSE)
  library(rDGIdb)
}

################################
## Data given by the user ######
################################
AnottateDrugsForthisCharacter <- function(mygenes){
  myquery <- queryDGIdb(mygenes)
  return(detailedResults(myquery)[,c(1,3,4)])
}

