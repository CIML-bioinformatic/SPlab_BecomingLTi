###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "04_Dynamics_Monocle"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

# File containing the raw epxression data after filtration of cells from QC and contamination analysis
PATH_RAW_EXPRESSION_MATRIX_FILE = file.path( PATH_PROJECT_OUTPUT, "03_GlobalHeterogeneity_NoContamination", paste0( outputFilesPrefix, "QCfiltered_Contaminationfiltered_raw_expression_matrix.csv"))

# File containing the t-SNE coordinates computed in the previous step
PATH_TSNE_EMBEDDING_FILE = file.path( PATH_PROJECT_OUTPUT, "03_GlobalHeterogeneity_NoContamination", paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_Tsnecoord.tsv"))
PATH_UMAP_EMBEDDING_FILE = file.path( PATH_PROJECT_OUTPUT, "03_GlobalHeterogeneity_NoContamination", paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_Umapcoord.tsv"))

VAR_GENES_MAX_NUMBER = 100


