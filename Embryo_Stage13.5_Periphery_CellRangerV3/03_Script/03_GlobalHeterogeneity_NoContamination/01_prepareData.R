# #########################################
# This script reads the sc10x data from
# QC and filtering analysis
# #########################################


# READ DATA
# ---------

## @knitr loadData

# Load the data after QC
sc10x.rna.seurat = readRDS( file.path( PATH_SEURAT_AFTER_QC_OBJECT))
cat("<br>Number of cells in loaded data:", length( Cells( sc10x.rna.seurat)))

# Remove contination as identified in previous step
cells_to_exclude = read.table( file = PATH_CONTAMINATION_CELL_FILE, sep="\t", header = TRUE)
sc10x.rna.seurat = subset( sc10x.rna.seurat, cells = setdiff( Cells( sc10x.rna.seurat), cells_to_exclude$cellid))
cat("<br>Number of cells after removing contamination:", length( Cells( sc10x.rna.seurat)))


# Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x.rna.seurat))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));

# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( sc10x.rna.seurat))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)

# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit);

