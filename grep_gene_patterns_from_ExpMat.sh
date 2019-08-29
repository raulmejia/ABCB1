# this program extract the rows of a matrix that match with a certain pattern
# 1) folder_to_search
# 2) patter_of_your_matrices
# 3) patter_of_your genes
# 4) Label_for save you results

find $1 -name $2 -exec sh -c 'awk "NR==1 || /'"$3"'/"  {} > {}_'"$4"'_genes.tsv' \;

#find ../Results/Splited/Matrices/ -name *Subexpression_matrix*_.tsv -exec sh -c 'awk "NR==1 || /ABC/"  {} > {}_ABC_genes.tsv' \;


