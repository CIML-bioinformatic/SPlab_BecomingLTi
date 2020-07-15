###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#


#### General

GLOBAL_DESCRIPTION = "Becoming LTi - Embryo Stage13.5 Periphery - CellRanger V3"

SCIENTIFIC_GROUP = "SPlab"
SCIENTIFIC_PROJECT_NAME = "BecomingLTi"
EXPERIMENT_PROJECT_NAME = "Embryo_Stage13.5_Periphery_CellRangerV3"

SAMPLE_NAME = "Embryo_Stage13.5_Periphery"
SAMPLE_ID = "106310444"

#### List of genes to be monitored

SIGNATURE_CLP = c( "Flt3", "Id2", "Ly6a")
SIGNATURE_alphaLP = c( "Id2", "Itga4", "Itgb7", "Rorc", "Icos", "Ly6a")
SIGNATURE_LTiP = c( "Cxcr6", "Rorc", "Zbtb16", "Flt3", "Id2", "Icos", "Tox")
SIGNATURE_LTi = c( "Lta", "Il22", "Il23r", "Ccr6", "Cd4", "Cxcr5", "Id2", "Il7r", "Ltb")
SIGNATURE_ILCP = c( "Zbtb16", "Rorc", "Ccr6", "Icos", "Tox")
SIGNATURE_ILC1 = c( "Il2rb", "Rorc", "Eomes", "Tbx21")
SIGNATURE_ILC2 = c( "Il13", "Rora", "Il1rl1", "Rorc", "Klrg1", "Il18r1", "Il17rb")
SIGNATURE_ILC3 = c( "Rorc", "Ccr6", "Lta", "Cd4", "Ncr1", "Il18r1", "Il22")
SIGNATURE_NK = c( "Eomes")

# Genes related to proliferation
SIGNATURE_PROLIFERATION = c( "Ccnb1", "Ccnb2", "Ccna2", "Cdc20", "Plk1", "Ube2c", "Cenpf", "Cenpa", "Aurka", "Mki67")

# List all the MONITORED_GENES
MONITORED_GENES = list(  CLP = SIGNATURE_CLP,
			 alphaLP = SIGNATURE_alphaLP,
                         LTiP = SIGNATURE_LTiP,
                         LTi = SIGNATURE_LTi,
                         ILCP = SIGNATURE_ILCP,
                         ILC1 = SIGNATURE_ILC1,
                         ILC2 = SIGNATURE_ILC2,
                         ILC3 = SIGNATURE_ILC3,
			 NK = SIGNATURE_NK,
                         PROLIFERATION = SIGNATURE_PROLIFERATION)

#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "/mnt/NAS7", 
                                        SCIENTIFIC_GROUP,
                                        "BIOINFO_PROJECT",
                                        SCIENTIFIC_PROJECT_NAME,
                                        EXPERIMENT_PROJECT_NAME)

PATH_PROJECT_RAWDATA = file.path( PATH_PROJECT, "00_RawData")
PATH_PROJECT_OUTPUT = file.path( PATH_PROJECT, "05_Output")

# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_",     
                            EXPERIMENT_PROJECT_NAME, "_"
                            #startTimeFileName, "_",
                            #paramsHash, "_"
)

#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



