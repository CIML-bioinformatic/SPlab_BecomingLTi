###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "01_QC"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Quality Control and normalization"

PATH_CELLRANGER_ANALYSIS = file.path( PATH_PROJECT_RAWDATA, "CellRanger_Counts")
PATH_MOUSE_RNA_COUNT_TABLE = file.path( PATH_CELLRANGER_ANALYSIS, "outs/filtered_feature_bc_matrix")

#### Filtering / Normalization

LOAD_MIN_CELLS     = 3;    # Retain cells with at least this many features (annotations)
LOAD_MIN_FEATURES  = 200;  # Retain annotations appearing in at least this many cells

# Cells with number of UMIs outside the range will be excluded
FILTER_UMI_MIN     = 2000;
FILTER_UMI_MAX     = NULL;

# Cells with number of genes outside the range will be excluded
FILTER_FEATURE_MIN = NULL;
FILTER_FEATURE_MAX = 7500;

# Cells with percentage of mitochondrial gene above this value will be excluded
FILTER_MITOPCT_MAX = 10;

# Cells with percentage of ribosomal gene below this value will be excluded
FILTER_RIBOPCT_MIN = 10;

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)

