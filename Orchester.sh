Rscript Split_matrix_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/Matrices/ TCGA

Rscript Controls_pasted_with_Submatrices_generated_by_labels.R ../Data/expMatrix_TCGA_cBioPortal_no_males_withindicator.txt ../Data/Labels_Controls_and_Normal_separated_TCGA.txt  ../Results/Splited/Matrices/ TCGA

for names in ../Results/Splited/Matrices/Subexpression_matrix_*from_TCGA_.tsv ; do
new_names=${names%_.tsv}
	echo $new_names;
done;

for names in ../Results/Splited/Matrices/Subexpression_matrix_*from_TCGA_.tsv ; do
new_names=${names%_.tsv}
	echo $new_names;
done;



new_name=${names%".tsv"}
$foo=${foo%"$suffix"}
