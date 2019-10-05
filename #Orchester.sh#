##############################
#######
####### MAINZ
#######
##############################
Rscript getting_MAINZ_DB.R ../Results/Splited/Submatrices_ohne_controls/ ../Results/Lehmann-STROMA4/MAINZ/ ./ MAINZ-only-Basals  

##############################
#######
####### TCGA
#######
##############################
Rscript Split_matrix_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/Submatrices_ohne_controls/ TCGA

Rscript Controls_pasted_with_Submatrices_generated_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/SubMatrices_with_controls/ TCGA # Generating submatrix of cases and controls (A)

Rscript splitdf_by_gene_level_in_tumours_Command_Args.R ../Results/Splited/SubMatrices_with_controls/Control_with_subexpression_matrix_of_Basal_from_TCGA_.tsv  ABCB1 ../Results/Matrices_splited_by_gene/ABCB1/ TCGA_Basal 

####################################
## Runing Pathifier (B)
####################################
## Getting the updated KEGG database
Rscript UptodateKEGG.R ../Results/KEGGDB/

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_under_percentile_25_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_low.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_under_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_high.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./  ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_above_percentile_50_stbl_10

Rscript Pathifier_Args_10stabilizing_Filtervalue3_75.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_high.tsv ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ./ ../Results/Pathifier/Basal/TCGA/ TCGA_Basal_above_percentile_75_stbl_10

# spliting pathifier's matrices
for names in ../Results/Pathifier/Basal/TCGA/TCGA_Basal_*median_PDSz_ordered_matrix.txt  ; do
new_names=${names%_.tsv}
	head -n 30 $names > $new_names"_top_30.tsv"
	head -n 40 $names > $new_names"_top_40.tsv"
	head -n 50 $names > $new_names"_top_50.tsv"
	echo $new_names;
done;

####################################
## DGE genes (C)
####################################
## only log 2 transformed

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05 0.6 0.05 7

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_low.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_under_per50_only_log2transformed_lfc2_of0_6_padjof0_05 0.6 0.05 7

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_high.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_above_per50_only_log2transformed_lfc2_of0_6_padjof0_05 0.6 0.05 7

Rscript DESeq2_with_log2_transformed_data_only.R ../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_high.tsv ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ./ ../Results/DEG/TCGA/log2only/ _DGE_TCGA_Basal_ABCB1_above_per75_only_log2transformed_of0_6_padjof0_05 0.6 0.05 7

#########################################################
## Big DF Pathway Target Drug Interaction LogFC   (D) ###
#########################################################

Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/TCGA/ TCGA_Basal_under_percentile_25_top20_DEG_log2only # top 20

--- Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix_Top20.txt ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/TCGA/ TCGA_Basal_under_percentile_25_top20_DEG_log2only

Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix.txt_top_50.tsv ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/TCGA/ TCGA_Basal_under_percentile_25_top50_DEG_log2only # top 50

Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix.txt ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_under_per25_only_log2transformed_lfc2_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/TCGA/ TCGA_Basal_under_percentile_25_ALLpw_DEG_log2only  # ALL pathifier pathways 


##### Above 75
Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/TCGA/TCGA_Basal_above_percentile_75_stbl_10_median_PDSz_ordered_matrix.txt_top_50.tsv ../Results/DEG/TCGA/log2only/padj10_3_lfc1_results_DESeq_DGE_TCGA_Basal_ABCB1_above_per75_only_log2transformed_of0_6_padjof0_05.tsv ./ ../Results/BigDfPTD/TCGA/ TCGA_Basal_above_percentile_75_top50_DEG_log2only # top 50 above 75

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

for names in ../Results/Pathifier/Basal/METABRIC/METABRIC_Basal_*median_PDSz_ordered_matrix.txt  ; do
new_names=${names%_.tsv}
	head -n 30 $names > $new_names"_top_30.tsv"
	head -n 40 $names > $new_names"_top_40.tsv"
	head -n 50 $names > $new_names"_top_50.tsv"
	echo $new_names;
done;


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
Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/METABRIC/METABRIC_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix.txt_top_50.tsv ../Results/DEG/METABRIC/DGE_METABRIC_Basal_ABCB1_under_per25_lfc0_6_padj0_050.050.6Case_Controls_ID_no_ILMN_.txt ./ ../Results/BigDfPTD/METABRIC/ METABRIC_Basal_under_percentile_25_top50_DEG_0_05_lfc0_6

Rscript DfPathways__and_PwDefinitions_and_dfDEG_2_Df_of_PwGenesDrugs.R ../Results/KEGGDB/KEGG_pathways_in_df_genesymbol.tsv ../Results/Pathifier/Basal/METABRIC/METABRIC_Basal_above_percentile_75_stbl_10_median_PDSz_ordered_matrix.txt_top_50.tsv ../Results/DEG/METABRIC/DGE_METABRIC_Basal_ABCB1_above_per75_lfc0_6_padj0_050.050.6Case_Controls_ID_no_ILMN_.txt ./ ../Results/BigDfPTD/METABRIC/ METABRIC_Basal_above_percentile_75_top50_DEG_0_05_lfc0_6

################################
### Lehmann's subtypes
#################################
Rscript Lehmann_subtypes_numerical_results.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_METABRIC_.tsv ./ ../Results/Lehmann-STROMA4/METABRIC/ METABRIC-only-Basals-lehmann_s_and_properties no 7 

Rscript Lehmann_subtypes_numerical_results.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_TCGA_.tsv ./ ../Results/Lehmann-STROMA4/TCGA/ TCGA-only-Basals-lehmann_s_and_properties yes 7

################################
### Clustering of ABC genes
################################
# Extracting only the ABC genes
sh grep_gene_patterns_from_ExpMat.sh ../Results/Splited/Submatrices_ohne_controls/ Subexpression_matrix*Basal*.tsv ^ABC only_ABC_transporters # Extracting only the ABC transporter genes

sh grep_gene_patterns_from_ExpMat.sh ../Results/Splited/SubMatrices_with_controls/ Control*ubexpression_matrix*Basal*.tsv ^ABC only_ABC_transporters # Extracting only the ABC transporter genes  Basal with Controls

# Heatmap and clustering of that samples
Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_METABRIC_.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/METABRIC/Lehmann\'s_Subt_and_properties_Numerical_METABRIC-only-Basals-lehmann_s_and_properties.tsv ./ ../Results/Clusterig/METABRIC/Heatmaps/Basal_Only_ABCS/ Heatmap_METABRIC_Basal_ABC_transporters_Lehmann_clasification log2-transformation-no

Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_TCGA_.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/TCGA/Lehmann\'s_Subt_and_properties_Numerical_TCGA-only-Basals-lehmann_s_and_properties.tsv ./ ../Results/Clusterig/TCGA/Heatmaps/Basal_Only_ABCS/ Heatmap_TCGA_Basal_ABC_transporters_Lehmann_clasification log2-transformation-yes

Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_MAINZ_HCGS_.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/MAINZ/Lehmann\'s_Subt_and_properties_Numerical_MAINZ-only-Basals.tsv ./ ../Results/Clusterig/MAINZ/Heatmaps/Basal_Only_ABCS/ Heatmap_MAINZ_Basal_ABC_transporters_Lehmann_clasification log2-transformation-no

Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_MAINZ_HCGS_collapsed.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/MAINZ/Lehmann\'s_Subt_and_properties_Numerical_MAINZ-only-Basals.tsv ./ ../Results/Clusterig/MAINZ/Heatmaps/Basal_Only_ABCS/ Heatmap_MAINZ_Basal_ABC_transporters_Lehmann_clasification_collapsed_median log2-transformation-no

Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_MAINZ_HCGS_er_neg_.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/MAINZ/Lehmann\'s_Subt_and_properties_Numerical_MAINZ-only-Basalsreneg_.tsv ./ ../Results/Clusterig/MAINZ/Heatmaps/Basal_Only_ABCS/ Heatmap_MAINZ_Basal_ABC_transporters_Lehmann_clasification_erneg log2-transformation-no

Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_MAINZ_HCGS_collapsed_er_neg.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/MAINZ/Lehmann\'s_Subt_and_properties_Numerical_MAINZ-only-Basalsreneg_.tsv ./ ../Results/Clusterig/MAINZ/Heatmaps/Basal_Only_ABCS/ Heatmap_MAINZ_Basal_ABC_transporters_Lehmann_clasification_collapsed_medianerneg log2-transformation-no

Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_MAINZ_HCGS_er_neg_ERBB_normal_.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/MAINZ/Lehmann\'s_Subt_and_properties_Numerical_MAINZ-only-Basalsreneg_ERBB_normal_.tsv ./ ../Results/Clusterig/MAINZ/Heatmaps/Basal_Only_ABCS/ Heatmap_MAINZ_Basal_ABC_transporters_Lehmann_clasification_ERBB_Normal_erneg_ log2-transformation-no

Rscript Pheatmap.R ../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_MAINZ_HCGS_collapsed_er_neg_ERBB_normal_.tsv_only_ABC_transporters_genes.tsv ../Results/Lehmann-STROMA4/MAINZ/Lehmann\'s_Subt_and_properties_Numerical_MAINZ-only-Basalsreneg_ERBB_normal_.tsv ./ ../Results/Clusterig/MAINZ/Heatmaps/Basal_Only_ABCS/ Heatmap_MAINZ_Basal_ABC_transporters_Lehmann_clasification_ERBB_erneg_collapsed_medianerneg log2-transformation-no


Subexpression_matrix_Basal_from_MAINZ_HCGS_er_neg_ERBB_normal_.tsv_only_ABC_transporters_genes.tsv


../Results/Splited/Submatrices_ohne_controls/Subexpression_matrix_Basal_from_MAINZ_HCGS_collapsed_er_neg.tsv_only_ABC_transporters_genes.tsv
#


#grep -nr Breast ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_25_stbl_10_median_PDSz_ordered_matrix.txt ../Results/Pathifier/Basal/TCGA/TCGA_Basal_under_percentile_50_stbl_10_median_PDSz_ordered_matrix.txt ../Results/Pathifier/Basal/TCGA/TCGA_Basal_above_percentile_50_stbl_10_median_PDSz_ordered_matrix.txt ../Results/Pathifier/Basal/TCGA/TCGA_Basal_above_percentile_75_stbl_10_median_PDSz_ordered_matrix.txt 

#for names in ../Results/Splited/SubMatrices_with_controls/Subexpression_matrix_*from_TCGA_.tsv ; do
#new_names=${names%_.tsv}
#	echo $new_names;
#done;

for names in ../Results/Splited/Matrices/Subexpression_matrix_*from_TCGA_.tsv ; do
new_names=${names%_.tsv}
	echo $new_names;
done;

