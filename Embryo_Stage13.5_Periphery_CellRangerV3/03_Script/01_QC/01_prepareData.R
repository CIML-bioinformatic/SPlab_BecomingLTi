# #########################################
# This script reads and filters sc10x  data
# #########################################


# -----------------------------------------------
# -----------------------------------------------
# READ DATA
# -----------------------------------------------
# -----------------------------------------------

## @knitr loadData

sc10x.rna.seurat = CreateSeuratObject( counts = Read10X( PATH_MOUSE_RNA_COUNT_TABLE), 
                            min.cells = LOAD_MIN_CELLS, 
                            min.features = LOAD_MIN_FEATURES, 
                            project = SAMPLE_NAME);


# Attribute a numeric ID to each cell (easier than barcode)
sc10x.rna.seurat[[ "numID"]] = 1:length( Cells( sc10x.rna.seurat));


# # Remove eventual NULL (empty) list elements from list of genes to be monitored
# monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
# if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
# MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];
# 
# # Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
# #matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x.rna.seurat))));
# matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x.rna.seurat))));
# monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
# if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));
# # Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
# MONITORED_GENES = relist( rownames( GetAssayData( sc10x.rna.seurat))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)
# # Finally remove names that did not match
# MONITORED_GENES = lapply( MONITORED_GENES, na.omit);


# -----------------------------------------------
# -----------------------------------------------
# FILTER DATA
# -----------------------------------------------
# -----------------------------------------------

## @knitr filterData

### Identify mitocondrial genes in matrix
mito.genes = grep( pattern = "^mt-", x = rownames(GetAssayData(object = sc10x.rna.seurat, slot = "counts")), value = TRUE, ignore.case = TRUE)
if(length(mito.genes)==0){
  warning( "No mitochondrial genes could be identified in this dataset.");
} else{
  # Compute percentage of mitochondrial transcripts in each cell and 
  # Add the mitocondrial gene percentage as meta information in the Seurat object
  # Note : Specify features instead of pattern so we can use ignore.case in grep call
  sc10x.rna.seurat[["percent.mito"]] <- PercentageFeatureSet(sc10x.rna.seurat, features=mito.genes)
}

## Identify Ribisomal genes
# Compute percentage of mitochondrial transcripts in each cell and 
# Add the ribosomal gene percentage as meta information in the Seurat object
# Note : Specify features instead of pattern so we can use ignore.case in grep call
ribo.genes = grep( pattern = "^Rp[sl][[:digit:]]", x = rownames(GetAssayData(object = sc10x.rna.seurat, slot = "counts")), value = TRUE, ignore.case = FALSE)
if(length(ribo.genes)==0){
  warning( "No ribosomal genes could be identified in this dataset.");
} else{
  sc10x.rna.seurat[["percent.ribo"]] <- PercentageFeatureSet(sc10x.rna.seurat, features=ribo.genes)
}

### Identify cells that will be rejected based on specified thresholds
# ....................................................................

# Identify cells to reject based on UMI numbers
nUMI.drop = logical( length( Cells(sc10x.rna.seurat)));
if( ! is.null( FILTER_UMI_MIN)){
  nUMI.drop = nUMI.drop | (sc10x.rna.seurat[["nCount_RNA", drop=TRUE]] < FILTER_UMI_MIN);
}
if( ! is.null( FILTER_UMI_MAX)){
  nUMI.drop = nUMI.drop | (sc10x.rna.seurat[["nCount_RNA", drop=TRUE]] > FILTER_UMI_MAX);
}

# Identify cells to reject based on number of expressed genes
nGene.drop = logical( length(Cells(sc10x.rna.seurat)));
if( ! is.null( FILTER_FEATURE_MIN)){
  nGene.drop = nGene.drop | (sc10x.rna.seurat[["nFeature_RNA", drop=TRUE]] < FILTER_FEATURE_MIN);
}
if( ! is.null( FILTER_FEATURE_MAX)){
  nGene.drop = nGene.drop | (sc10x.rna.seurat[["nFeature_RNA", drop=TRUE]] > FILTER_FEATURE_MAX);
}

# Identify cells with high percentage of mitocondrial genes
mito.drop = logical( length(Cells(sc10x.rna.seurat)));
if( length(mito.genes) && (! is.null(FILTER_MITOPCT_MAX))){
  mito.drop = (sc10x.rna.seurat[["percent.mito", drop=TRUE]] > FILTER_MITOPCT_MAX);
}

# Identify cells with low percentage of ribosomal genes
ribo.drop = logical( length(Cells(sc10x.rna.seurat)));
if( length( ribo.genes) && (! is.null( FILTER_RIBOPCT_MIN))){
  ribo.drop = (sc10x.rna.seurat[["percent.ribo", drop=TRUE]] < FILTER_RIBOPCT_MIN);
}

### Plot distributions of #UMIs, #Genes and %Mito among cells
# ....................................................................

# Gather data to be visualized together (cell name + numID + metrics)
cellsData = cbind( "Cell" = colnames( sc10x.rna.seurat), # cell names from rownames conflict with search field in datatable
                   sc10x.rna.seurat[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x.rna.seurat[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x.rna.seurat[["percent.ribo", drop = TRUE]], digits = 5)) else NULL);

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste, c( "",
                                            "Cell ID: ", "# UMIs: ", "# Genes: ",
                                            if(length( mito.genes)) "% Mito: ",
                                            if(length( ribo.genes)) "% Ribo: "), cellsData, sep=""), sep="\n"))

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
panelWidth = 180;
panelHeight = 400;
lypanel_umis  = plotViolinJitter(cbind( cellsData, outliers = nUMI.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~nCount_RNA, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "# UMIs", thresholdHigh = FILTER_UMI_MAX, thresholdLow = FILTER_UMI_MIN, panelWidth = panelWidth, panelHeight = panelHeight);
lypanel_genes = plotViolinJitter(cbind( cellsData, outliers = nGene.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~nFeature_RNA, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "# Genes", thresholdHigh = FILTER_FEATURE_MAX, thresholdLow = FILTER_FEATURE_MIN, panelWidth = panelWidth, panelHeight = panelHeight);
lypanel_mitos = if(length(mito.genes)) plotViolinJitter(cbind( cellsData, outliers = mito.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~percent.mito, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "% Mito", thresholdHigh = FILTER_MITOPCT_MAX, panelWidth = panelWidth, panelHeight = panelHeight) else NULL;
lypanel_ribos = if(length(ribo.genes)) plotViolinJitter(cbind( cellsData, outliers = ribo.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~percent.ribo, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "% Ribo", thresholdHigh = FILTER_RIBOPCT_MIN, panelWidth = panelWidth, panelHeight = panelHeight) else NULL;


# Set panels as a list and define plotly config
distribPanelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
distribPanelsList = lapply( distribPanelsList, config, displaylogo = FALSE, 
                            toImageButtonOptions = list(format='svg'), 
                            modeBarButtons = list(list('toImage'), list('zoom2d', 'pan2d', 'resetScale2d')));

# Control layout using flex because subplot is limited regarding plot sizing and alignment
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          lapply( distribPanelsList, div, style = paste("flex : 0 0 auto; margin: 5px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)

### Plot scatterplots of pairs of #UMIs, #Genes and %Mito among cells
# ....................................................................

lypanel_feature_vs_numi= plotScatterPlotly(cbind( cellsData, outliers = nUMI.drop | nGene.drop), xAxisFormula = ~nCount_RNA, yAxisFormula = ~nFeature_RNA, xAxisTitle = "# UMI", yAxisTitle = "# Features", pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, panelWidth = 380, panelHeight = 380);
lypanel_ribo_vs_mito = if(length(ribo.genes) && length( mito.genes)) plotScatterPlotly(cbind( cellsData, outliers = ribo.drop | mito.drop), xAxisFormula = ~percent.mito, yAxisFormula = ~percent.ribo, xAxisTitle = "% Mito genes", yAxisTitle = "% Ribo genes", pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, panelWidth = 380, panelHeight = 380) else NULL;

# Set panels as a list and define plotly config
scatterPanelsList = list( lypanel_ribo_vs_mito, lypanel_feature_vs_numi );
scatterPanelsList = lapply( scatterPanelsList, config, displaylogo = FALSE, 
                            toImageButtonOptions = list(format='svg'), 
                            modeBarButtons = list(list('toImage'), list('zoom2d', 'pan2d', 'resetScale2d')));

# Control layout using flex because subplot is limited regarding plot sizing and alignment
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          lapply( scatterPanelsList, div, style = paste("flex : 0 0 auto; margin: 5px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)

### Display information on the removed cells
# ....................................................................

cat( "<br>Number of cells removed based on number of UMIs:", sum( nUMI.drop));
cat( "<br>Number of cells removed based on number of genes:", sum( nGene.drop));
if(exists( "mito.drop")) cat( "<br>Number of cells removed based on high percentage of mitochondrial transcripts:", sum( mito.drop));
if(exists( "ribo.drop")) cat( "<br>Number of cells removed based on low percentage of ribosomal transcripts:", sum( ribo.drop));

### Identify the cells to exclude as union of cells with low nb UMI, low nb of expressed genes, 
### high percentage of mitocondrial genes or low percentage of ribosomal genes
# ....................................................................

sc10x.rna.seurat[["outlier"]] = (nUMI.drop | nGene.drop | (if(exists( "mito.drop")) mito.drop else FALSE) | (if(exists( "ribo.drop")) ribo.drop else FALSE));
cat("<br><br>Removed cells after filters:", sum( unlist(sc10x.rna.seurat[["outlier"]] )));
cat("<br>Remaining cells after filters:", sum( ! unlist(sc10x.rna.seurat[["outlier"]] )));


### Record which cells got rejected
# ....................................................................

# Export the excluded cells to file
excluded_cells_nbUMI_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_nbUMI.txt"))
write.table( data.frame( cellid = names( which( nUMI.drop))), file= excluded_cells_nbUMI_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of excluded cells by nbUMI to file: <BR>", excluded_cells_nbUMI_outfile)

excluded_cells_nbGene_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_nbGene.txt"))
write.table( data.frame( cellid = names( which( nGene.drop))), file= excluded_cells_nbGene_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of excluded cells by nbGene to file: <BR>", excluded_cells_nbGene_outfile)

excluded_cells_pctMito_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctMito.txt"))
write.table( data.frame( cellid = names( which( mito.drop))), file= excluded_cells_pctMito_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of excluded cells by pctMito to file: <BR>", excluded_cells_pctMito_outfile)

excluded_cells_pctRibo_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctRibo.txt"))
write.table( data.frame( cellid = names( which( ribo.drop))), file= excluded_cells_pctRibo_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of excluded cells by PctRibo to file: <BR>", excluded_cells_pctRibo_outfile)


### Exclude the outliers cells from the Seurat RNA object
# ....................................................................
cat("<BR><br>Number of cells in Seurat object before filtering:" , length( Cells( sc10x.rna.seurat)))
sc10x.rna.seurat = subset( sc10x.rna.seurat, subset = (outlier == FALSE))
cat("<br>Number of cells in Seurat object after filtering:" , length( Cells( sc10x.rna.seurat)))


# NORMALIZE DATA
# --------------

## @knitr normalizeData

# Normalize RNA data with log normalization
sc10x.rna.seurat = NormalizeData( object = sc10x.rna.seurat,
                                  normalization.method = DATA_NORM_METHOD,
                                  scale.factor = DATA_NORM_SCALEFACTOR,
                                  verbose = .VERBOSE);

sc10x.rna.seurat = ScaleData( object    = sc10x.rna.seurat,
                              do.center = DATA_CENTER,
                              do.scale  = DATA_SCALE,
                              vars.to.regress = DATA_VARS_REGRESS,
                              verbose = .VERBOSE);

# SAVE SEURAT OBJECT TO RDATA
# ---------------------------

## @knitr saveSeuratRObject

saveRDS( sc10x.rna.seurat, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sc10x.rna.seurat.RDS")))

