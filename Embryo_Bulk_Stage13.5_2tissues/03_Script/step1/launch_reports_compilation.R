# ##########################################################################
# This script is used to launch the compilation of the RMarkdown report
# independently of Rstudio interface
# ##########################################################################

WORKING_DIR = file.path( Sys.getenv( "WORKING_DIR"), "Embryo_Bulk_Stage13.5_2tissues")

STEP = "step1"

RAW_DATA_DIR = file.path( WORKING_DIR, "00_RawData")
SCRIPT_DIR = file.path( WORKING_DIR, "03_Script", STEP)
OUTPUT_DIR = file.path( WORKING_DIR, "05_Output", STEP)

rmarkdown::render( input = file.path( SCRIPT_DIR, "SPlab_ILCYou_PrimaryAnalysis.Rmd"),
                   output_dir = OUTPUT_DIR,
                   output_file  = "SPlab_ILCYou_PrimaryAnalysis.html",
                   quiet = FALSE)
