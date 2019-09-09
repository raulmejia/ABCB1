
##############################
#######
####### TCGA
#######
##############################

Rscript Split_matrix_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/Submatrices_ohne_controls/ TCGA

Rscript Controls_pasted_with_Submatrices_generated_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/SubMatrices_with_controls/ TCGA

Rscript splitdf_by_gene_level_in_tumours_Command_Args.R ../Results/Splited/SubMatrices_with_controls/Control_with_subexpression_matrix_of_Basal_from_TCGA_.tsv  ABCB1 ../Results/Matrices_splited_by_gene/ABCB1/ TCGA_Basal 

####################################
## Runing Pathifier
####################################
## Getting the updated KEGG database
Rscript UptodateKEGG.R ../Results/KEGGDB/

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_under_percentile_25_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_under_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_high.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_above_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_high.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_above_percentile_75_stbl_10

####################################
## DGE genes
####################################
## only log 2 transformed

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05 0.6 0.05 7

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_low.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_under_per50_only_log2transformed_lfc2_of0_6_padjof0_05 0.6 0.05 7

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_high.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_above_per50_only_log2transformed_lfc2_of0_6_padjof0_05 0.6 0.05 7

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_high.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_above_per75_only_log2transformed_of0_6_padjof0_05 0.6 0.05 7

## Big DF Pathway Target Drug Interaction LogFC 

Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/TCGA/ TCGA_Basal_under_percentile_25_top20_DEG_log2only

--- Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/TCGA/ TCGA_Basal_under_percentile_25_top20_DEG_log2only

####################
####################
#### METABRIC ######
####################
####################
Rscript Split_matrix_by_labels.R ../Data/joined_indicator_METABRIC.txt ../Data/Labels_Ctrl_and_NL_Recal_separated_METABRIC.txt ../Results/Splited/Submatrices_ohne_controls/ METABRIC

Rscript Controls_pasted_with_Submatrices_generated_by_labels.R ../Data/joined_indicator_METABRIC.txt ../Data/Labels_Ctrl_and_NL_Recal_separated_METABRIC.txt ../Results/Splited/SubMatrices_with_controls/ METABRIC

Rscript splitdf_by_gene_level_in_tumours_Command_Args.R ../Results/Splited/SubMatrices_with_controls/Control_with_subexpression_matrix_of_Basal_from_METABRIC_.tsv ABCB1 ../Results/Matrices_splited_by_gene/ABCB1/ METABRIC_Basal 

####################################
## Runing Pathifier
####################################
## Getting the updated KEGG database
Rscript UptodateKEGG.R ../Results/KEGGDB/

Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/METABRIC/ METABRIC_Basal_under_percentile_25_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/METABRIC/ METABRIC_Basal_under_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_high.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/METABRIC/ METABRIC_Basal_above_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_high.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/Basal/METABRIC/ METABRIC_Basal_above_percentile_75_stbl_10


####################################
## DGE genes
####################################
## only log 2 transformed
#Usar limma en lugar de DESeq2

Rscript DGE_limma2.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv ./ ../Results/DEG/METABRIC/ 0.6 0.05 METABRIC_Basal_ABCB1_under_per25_lfc0_6_padj0_05

Rscript DGE_limma2.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_low.tsv ./ ../Results/DEG/METABRIC/ 0.6 0.05 METABRIC_Basal_ABCB1_under_per50_lfc0_6_padj0_05

Rscript DGE_limma2.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_high.tsv ./ ../Results/DEG/METABRIC/ 0.6 0.05 METABRIC_Basal_ABCB1_above_per50_lfc0_6_padj0_05

Rscript DGE_limma2.R ../Results/Matrices_splited_by_gene/ABCB1/METABRIC_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_high.tsv ./ ../Results/DEG/METABRIC/ 0.6 0.05 METABRIC_Basal_ABCB1_above_per75_lfc0_6_padj0_05

## Big DF Pathway Target Drug Interaction LogFC 

Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/METABRIC/METABRIC_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt ../Results/DEG/METABRIC/DGE_METABRIC_Basal_ABCB1_under_per25_lfc0_6_padj0_050.050.6Case_Controls_ID_with_ILMN.txt ./ ../Results/BigDfPTD/METABRIC/ METABRIC_Basal_under_percentile_25_top20_DEG_



## Big DF Pathway Target Drug Interaction LogFC 

Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/METABRIC/ METABRIC_Basal_under_percentile_25_top20_DEG_log2only


Rscript DGE_limma2.R ../Data/toyMETABRIC_Controls_LumA.txt ./ ../Results/DEG/METABRIC/ 0.6 0.05 METABRIC_toy_Normals_LumisA




#grep -nr Breast ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix.txt ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_50_stbl_10_median_PDSz_ordered_matrix.txt ../Results/Pathifier/Basal/TCGA/TCGA_Basal_above_percentile_50_stbl_10_median_PDSz_ordered_matrix.txt ../Results/Pathifier/Basal/TCGA/TCGA_Basal_above_percentile_75_stbl_10_median_PDSz_ordered_matrix.txt 


#for names in ../Results/Splited/SubMatrices_with_controls/Subexpression_matrix_*from_TCGA_.tsv ; do
#new_names=${names%_.tsv}
#	echo $new_names;
#done;





for names in ../Results/Splited/Matrices/Subexpression_matrix_*from_TCGA_.tsv ; do
new_names=${names%_.tsv}
	echo $new_names;
done;


