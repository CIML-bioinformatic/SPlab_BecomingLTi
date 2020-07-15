###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "03_GlobalHeterogeneity_NoContamination"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

PATH_SEURAT_AFTER_QC_OBJECT = file.path( PATH_PROJECT_OUTPUT, "01_QC", paste0( outputFilesPrefix, "sc10x.rna.seurat.RDS"))

PATH_CONTAMINATION_CELL_FILE = file.path( PATH_PROJECT_OUTPUT, "02_GlobalHeterogeneity", paste0( outputFilesPrefix, "contamination_excluded_all_cells.txt"))



# Maximum number of variable features to keep for PCA analysis
VARIABLE_FEATURES_MAXNB   = 2000;

# Maximum number of variable features to keep for table presentation
VARIABLE_FEATURES_SHOWTOP = 200;

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION = 1.0;

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 20;  # Number of dimensions to use from PCA results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 20;       # Number of marker genes to show in report and tables (NULL for all)




