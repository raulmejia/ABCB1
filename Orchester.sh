Rscript Split_matrix_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/Submatrices_ohne_controls/ TCGA

Rscript Controls_pasted_with_Submatrices_generated_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/SubMatrices_with_controls/ TCGA

Rscript splitdf_by_gene_level_in_tumours_Command_Args.R ../Results/Splited/SubMatrices_with_controls/Control_with_subexpression_matrix_of_Basal_from_TCGA_.tsv  ABCB1 ../Results/Matrices_splited_by_gene/ABCB1/ TCGA_Basal 

## Getting the updated KEGG database
Rscript UptodateKEGG.R ../Results/KEGGDB/

####################################
## Runing Pathifier
####################################
Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_under_percentile_25_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_under_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_high.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_above_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_high.tsv 
 ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_above_percentile_75_stbl_10


for names in ../Results/Splited/SubMatrices_with_controls/Subexpression_matrix_*from_TCGA_.tsv ; do
new_names=${names%_.tsv}
	echo $new_names;
done;





for names in ../Results/Splited/Matrices/Subexpression_matrix_*from_TCGA_.tsv ; do
new_names=${names%_.tsv}
	echo $new_names;
done;


