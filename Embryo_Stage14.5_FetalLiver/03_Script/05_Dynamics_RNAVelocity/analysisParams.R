###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "05_Dynamics_RNAVelocity"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

# Data from analysis of heterogenity with Seurat
PATH_SEURAT_NORMALIZED_EXPRESSION_MATRIX = file.path( PATH_PROJECT_OUTPUT, "03_GlobalHeterogeneity_NoContamination", paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_expression_matrix.csv"))
PATH_SEURAT_CLUSTERS = file.path( PATH_PROJECT_OUTPUT, "03_GlobalHeterogeneity_NoContamination", paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_Seurat_clusters.tsv"))
PATH_SEURAT_PCA_EMBEDDING = file.path( PATH_PROJECT_OUTPUT, "03_GlobalHeterogeneity_NoContamination", paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_PCAcoord.tsv"))

# Data of cell type analysis with Monocle
PATH_SEURAT_TSNE_EMBEDDING_WITHCELLTYPE = file.path( PATH_PROJECT_OUTPUT, "04_Dynamics_Monocle", paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_Tsnecoord_WithCellType.tsv"))
PATH_SEURAT_UMAP_EMBEDDING_WITHCELLTYPE = file.path( PATH_PROJECT_OUTPUT, "04_Dynamics_Monocle", paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_Umapcoord_WithCellType.tsv"))


