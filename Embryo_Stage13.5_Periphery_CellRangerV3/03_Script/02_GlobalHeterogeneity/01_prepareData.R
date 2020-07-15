# #########################################
# This script reads the sc10x data from
# QC and filtering analysis
# #########################################


# READ DATA
# ---------

## @knitr loadData

sc10x.rna.seurat = readRDS( file.path( PATH_SEURAT_FILTERED_OBJECT))

cat("<br>Number of cells in loaded data:", length( Cells( sc10x.rna.seurat)))


# Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x.rna.seurat))));
matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x.rna.seurat))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));

# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( sc10x.rna.seurat))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)

# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit);

