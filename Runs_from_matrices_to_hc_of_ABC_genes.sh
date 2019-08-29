
# Split TCGA matrix by subtypes
Rscript Split_matrix_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt ../Results/Splited/Matrices/TCGA/ TCGA

# from the previous matrices extract the rows that contain ABC genes
sh grep_gene_patterns_from_ExpMat.sh ../Results/Splited/Matrices/TCGA/ *Subexpression_matrix*_.tsv ABC Only_ABC_transporters
