
# Split TCGA matrix by subtypes
Rscript Split_matrix_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ../Results/Splited/Matrices/TCGA/ TCGA

# from the previous matrices extract the rows that contain ABC genes
sh grep_gene_patterns_from_ExpMat.sh ../Results/Splited/Matrices/TCGA/ *Subexpression_matrix*_.tsv ABC Only_ABC_transporters

# Split TCGA matrix by subtypes
Rscript Split_matrix_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator_log2_transformed.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ../Results/Splited/Matrices/TCGA-log2/ TCGA-log2

# from the previous matrices extract the rows that contain ABC genes
sh grep_gene_patterns_from_ExpMat.sh ../Results/Splited/Matrices/TCGA-log2/ *Subexpression_matrix*_.tsv ABC Only_ABC_transporters

# Ploting trough heatmap
Rscript Heatmaping_hc_with_and_without_zscores.R ../Results/Splited/Matrices/TCGA-log2/Subexpression_matrix_Basal_from_TCGA-log2_.tsv_Only_ABC_transporters_genes.tsv ./ ../Results/Clusterig/TCGA-log2/ Basal_Only_ABC_transporters-log2
