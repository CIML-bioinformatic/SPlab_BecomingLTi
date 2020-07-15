# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# IDENTIFY CLUSTERS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_identifyClusters

sc10x.rna.seurat <- RunPCA( object   = sc10x.rna.seurat, 
                           features = VariableFeatures( sc10x.rna.seurat),
                           npcs     = PCA_NPC,
                           verbose  = .VERBOSE);

# # Compute the score of the cells according to group of monitored genes
# for( monitoredGroup in names( MONITORED_GENES)){
#   sc10x <- AddModuleScore( object = sc10x, 
#                            features = MONITORED_GENES[[ monitoredGroup]],
#                            ctrl = length(MONITORED_GENES[[ monitoredGroup]]), 
#                            name = monitoredGroup);
# }

# Identify clusters of cells by graph approach
sc10x.rna.seurat <- FindNeighbors(object  = sc10x.rna.seurat, 
                                   dims    = 1:30,
                                   verbose = .VERBOSE);

sc10x.rna.seurat <- FindClusters(object              = sc10x.rna.seurat, 
                                  # reduction.type     = "pca", 
                                  resolution         = FINDCLUSTERS_RESOLUTION, 
                                  random.seed        = 1511,
                                  temp.file.location = "/tmp/",
                                  verbose            = .VERBOSE);

# Show number of cells in each cluster
clusterCount = as.data.frame( table( ClusterNumber = sc10x.rna.seurat[[ "seurat_clusters" ]]), responseName = "CellCount");

# Create datatable
datatable( clusterCount, 
           class = "compact",
           rownames = FALSE,
           colnames = c("Cluster", "Nb. Cells"),
           options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          columnDefs = list( # Center all columns
                                            list( targets = 0:(ncol(clusterCount)-1),
                                            className = 'dt-center')),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          paging = FALSE, # Disable pagination (show all)
                          processing = TRUE, 
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE));


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# IDENTIFY CLUSTERS - SPLIT STATISTICS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_identifyClusters_splitStats

# Show UMIs Gens and Mitochondrial content split by cluster

# Gather data to be visualized together (cell name + numID + metrics + Cluster)
cellsData = cbind( "Cell" = colnames( sc10x.rna.seurat), # cell names from rownames conflict with search field in datatable
                   sc10x.rna.seurat[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = as.numeric(format(sc10x.rna.seurat[["percent.mito", drop = TRUE]], digits = 5)),
                   "Cluster" = Idents( sc10x.rna.seurat));

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste, c( "", "Cell ID: ", "# UMIs: ", "# Genes: ", "% Mito: "), cellsData[-ncol(cellsData)], sep=""), sep="\n"))

# Define size for panels (or assembled figure when using subplot)
panelWidth = 90 * nlevels(cellsData[["Cluster"]]);
panelHeight = 600;

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
lypanel_umis  = plotViolinJitter( cellsData, xAxisFormula = ~as.numeric(Cluster), yAxisFormula = ~nCount_RNA, colorFormula = ~scales::alpha(hue_pal()(nlevels( Cluster)), 0.7)[Cluster], yAxisTitle = "# UMIs", hoverText = hoverText, traceName = ~paste( "Cluster", Cluster), xTicklabels = levels(cellsData[["Cluster"]]), panelWidth = panelWidth, panelHeight = panelHeight);
lypanel_genes = plotViolinJitter( cellsData, xAxisFormula = ~as.numeric(Cluster), yAxisFormula = ~nFeature_RNA, colorFormula = ~scales::alpha(hue_pal()(nlevels( Cluster)), 0.7)[Cluster], yAxisTitle = "# Genes", hoverText = hoverText, traceName = ~paste( "Cluster", Cluster), xTicklabels = levels(cellsData[["Cluster"]]), panelWidth = panelWidth, panelHeight = panelHeight);
lypanel_mitos = plotViolinJitter( cellsData, xAxisFormula = ~as.numeric(Cluster), yAxisFormula = ~percent.mito, colorFormula = ~scales::alpha(hue_pal()(nlevels( Cluster)), 0.7)[Cluster], yAxisTitle = "% Mito", hoverText = hoverText, traceName = ~paste( "Cluster", Cluster), xTicklabels = levels(cellsData[["Cluster"]]), panelWidth = panelWidth, panelHeight = panelHeight);

# Set panels as a list, define plotly config, and remove legend
panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos);
panelsList = lapply(panelsList, config, displaylogo = FALSE, 
                    toImageButtonOptions = list(format='svg'), 
                    modeBarButtons = list(list('toImage'), 
                                          list('zoom2d', 'pan2d', 'resetScale2d')));

# Group plotly violin/jitter panels so we can sycnchronise axes and use highlight on the full set
plotPanels = layout( subplot( panelsList,
                              nrows = 3,
                              shareX = TRUE,
                              titleY = TRUE),
                     xaxis = list(title = "Seurat Cluster",
                                  showgrid = TRUE,
                                  tickvals = seq(nlevels(cellsData[["Cluster"]])),
                                  ticktext = levels(cellsData[["Cluster"]])),
                     showlegend = FALSE, # Remove eventual legends (does not mix well with subplot and highlight)
                     autosize = TRUE);

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          div(plotPanels, style = paste("flex : 0 0 auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PCA - COMPUTATION
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_pca

# Plot PCA, highlighting seurat clusters for combination of dimensions
pcaClustersList = apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims) 
          {
                    ( DimPlot( object=sc10x.rna.seurat, reduction="pca", dims = dims) + 
                        scale_color_discrete( name = "Cluster") +  # Fix legend title
                        theme( legend.position = "none"));         # Remove legend
          });
# Plot figures (ggplot)
resNULL = lapply(pcaClustersList, print);


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PCA - HIGHLIGHT UMI COUNTS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_pca_umisCounts

# Same PCA plots but highlight UMI counts
pcaUmisList = apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims) { ( FeaturePlot( sc10x.rna.seurat, feature = "nCount_RNA", reduction = "pca", dims = dims) +
                                                                            theme( legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"), plot.title = element_blank())); });
# Plot figures (ggplot)
resNULL = lapply(pcaUmisList, print);


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PCA : highlight gene counts
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_pca_genesCounts

# Same PCA plots but highlight feature counts
pcaGenesList = apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims) { ( FeaturePlot( sc10x.rna.seurat, feature = "nFeature_RNA", reduction = "pca", dims = dims) + theme( legend.position = "none", 
                                                                                                                                                              plot.margin = margin(0, 0, 0, 0, "cm"),
                                                                                                                                                              plot.title = element_blank())); });
# Plot figures (ggplot)
resNULL = lapply(pcaGenesList, print);



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PCA to features correlation
# ---------------------------------------------------------------
# ---------------------------------------------------------------
## @knitr heterogeneity_pca_correlations

# Isolate the PCA location of cells in first dimensions, together with UMI and genes counts for correlation analysis (makes use of cbind recycling to repeat values for each stacked PC)
relationToPC = suppressWarnings( cbind( 
  stack( as.data.frame( Embeddings( sc10x.rna.seurat, reduction = "pca")[ rownames( sc10x.rna.seurat[[]]), paste0( "PC_", 1:PCA_PLOTS_NBDIMS) ])),
  sc10x.rna.seurat[[ c( "nCount_RNA", "nFeature_RNA") ]],
  Cluster = Idents( sc10x.rna.seurat)));

# Plot relationship of UMIs and genes counts with PCs(combine plots using '/' from 'patchwork' lib)
print((ggplot( data = relationToPC, aes( x = values, y = nCount_RNA)) 
       + facet_wrap( ~ind) 
       + stat_binhex( bins = 60)
       #+ geom_point(aes(col = Cluster), alpha = 0.5) 
       + geom_smooth( method = 'lm') 
       + stat_cor(method="spearman") 
       + ylab("# UMIs")
       + theme(axis.title.x = element_blank(), 
               axis.text.x = element_blank(), 
               axis.ticks.x= element_blank(), 
               plot.margin= unit(c(1, 1, -0.5, 0.5), "lines")))
       / 
       (ggplot( data = relationToPC, aes( x = values, y = nFeature_RNA)) 
        + facet_wrap( ~ind) 
        + stat_binhex( bins = 60)
        #+ geom_point(aes(col = Cluster), alpha = 0.5) 
        + geom_smooth( method = 'lm') 
        + stat_cor(method="spearman") 
        + xlab("PC values") 
        + ylab("# Genes")));



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PCA : look at pCA loadings
# ---------------------------------------------------------------
# ---------------------------------------------------------------
## @knitr heterogeneity_pca_loadings

# Plot PCA loadings
pcaLoadingsList = apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  namesPC=paste0( "PC_", dims);
  # Get the loading values for concerned PCs
  loadingsMatrix = Loadings( sc10x.rna.seurat, reduction = "pca")[ , namesPC ];
  # Sort features by average absolute value and convert as DF with features names as column
  loadingsMatrix = head( loadingsMatrix[ order( apply( loadingsMatrix, 1, function(x){ mean( abs( x)) }), decreasing = TRUE), ], PCA_PLOTS_NBFEATURES);
  loadingsDF = data.frame( loadingsMatrix, features = rownames( loadingsMatrix));

  # Define symmetric and consistent axes for group of plots
  axesLimit = max( abs( loadingsMatrix));

  # Plot arrows and features name
  return(ggplot( data = loadingsDF, aes_( x = as.name( namesPC[1]), y = as.name( namesPC[2]))) +
    coord_cartesian( xlim = c( -axesLimit, axesLimit), ylim = c( -axesLimit, axesLimit)) +
    geom_text_repel( aes( label = features), max.iter = 10000) +
    geom_segment( x = 0 , y = 0, aes_( xend = as.name(namesPC[1]), yend = as.name(namesPC[2])), col = "#00000044", arrow = arrow( length = unit( 2, "mm"))));
});

# Plot figures (ggplot)
resNULL = lapply(pcaLoadingsList, print)



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# DIMENSIONAL REDUCTION (TSNE/UMAP)
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_dimReduc
sc10x.rna.seurat = RunTSNE( sc10x.rna.seurat, dims = 1:DIMREDUC_USE_PCA_NBDIMS, seed.use = 1511);
sc10x.rna.seurat = RunUMAP( sc10x.rna.seurat, dims = 1:DIMREDUC_USE_PCA_NBDIMS, seed.use = 1511);


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# DIMENSION REDUCTION - INTERACTIVE PLOT - Clusters
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_dimReduc_interactivePlot_clusters

# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

Idents( sc10x.rna.seurat) = "seurat_clusters"

## Do it in pure plotly
dimReducData = cbind(as.data.frame(Embeddings(sc10x.rna.seurat, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     Cluster = as.numeric(as.character(Idents( sc10x.rna.seurat))));

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReduc = plot_ly(dimReducData, 
                       x = as.formula( paste( "~", colnames( dimReducData)[1])), 
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       name = ~paste("Cluster", Cluster),
                       color = ~Cluster,
                       colors = hue_pal()( nlevels( Idents( sc10x.rna.seurat))),
                       height = 600) %>%
  add_markers( marker = list( size = 5,
                              opacity = 0.6),
               text = hoverText,
               hoverinfo = "text+name") %>%
  # Compute median of coordinates for each group and use it as new data for Cluster text annotations (type='scatter', mode='text')
  add_trace( data = as.data.frame( do.call( rbind, by( dimReducData, dimReducData[["Cluster"]], apply, 2, median))), 
             x = as.formula( paste( "~", colnames( dimReducData)[1])), 
             y = as.formula( paste( "~", colnames( dimReducData)[2])), 
             type = "scatter", 
             mode = "text", 
             text = ~Cluster,
             textfont = list( color = '#000000', size = 20),
             hoverinfo = 'skip') %>%
  hide_colorbar() %>%
  hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list(format='svg'),
          modeBarButtons = list(list('toImage'),
                                list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
     div(plotDimReduc, style = paste("flex : 0 0 600px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")))
#)




# ---------------------------------------------------------------
# ---------------------------------------------------------------
# DIMENSION REDUCTION - INTERACTIVE PLOT - Samples
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_dimReduc_interactivePlot_samples

# Interactive plot using plotly
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

Idents( sc10x.rna.seurat) = "MULTI_Sample"

dimReducData = cbind(as.data.frame(Embeddings(sc10x.rna.seurat, reduction = ifelse(exists("useReduction"), useReduction, "umap"))),
                     MULTI_Sample = as.character(Idents( sc10x.rna.seurat)));

# Reusing 'hoverText' summarizing cell info (generated for split violin plots)
plotDimReducBySamples = plot_ly(dimReducData, 
                       x = as.formula( paste( "~", colnames( dimReducData)[1])), 
                       y = as.formula( paste( "~", colnames( dimReducData)[2])),
                       name = ~paste("Sample", MULTI_Sample),
                       color = ~MULTI_Sample,
                       colors = hue_pal()( nlevels( Idents( sc10x.rna.seurat))),
                       width = 800,
                       height = 600) %>%
  add_markers( marker = list( size = 5,
                              opacity = 0.6),
               text = hoverText,
               hoverinfo = "text+name") %>%
  hide_colorbar() %>%
  #hide_legend() %>%
  config( displaylogo = FALSE,
          toImageButtonOptions = list(format='svg'),
          modeBarButtons = list(list('toImage'),
                                list('zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d')));

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
     div(plotDimReducBySamples, style = paste("flex : 0 0 600px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;")))
#)


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# DIMENSIONNALITY REDUCTION - SAMPLES DETAILS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_dimReduc_details_samples

# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

Idents( sc10x.rna.seurat) = "MULTI_Sample"

# Plot expression values of marker genes (TODO: message if list empty)
resNULL=lapply( sort( unique( as.character(Idents( sc10x.rna.seurat)))), function( sampleName)
{
  cat("###### Sample", sampleName, " \n");
  
  # Extract cells from current cluster to highlight them on a dimreduc plot
  sampleCells = which( Idents( sc10x.rna.seurat) == sampleName);
  # Create a dimreduc plot with current cluster highlighted (useful for large number of clusters)
  print( DimPlot( sc10x.rna.seurat, 
                  reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
                  cols="#44444422", 
                  cells.highlight = sampleCells, 
                  cols.highlight = "#FF000088", 
                  sizes.highlight = 1.5, 
                  order = sampleCells,  # Plot highlighted cells last
                  group.by=NULL) + 
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none") +
           ggtitle( paste( "Sample", sampleName)));
  
  cat(" \n \n"); # Required for '.tabset'
});

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# DIMENSIONALITY REDUCTION - THUMBNAIL OF CLUSTERS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_dimReduc_thumbnail
# Non-interactive with large labels for generating thumbnails when analyzing monitored and marker genes expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

Idents( sc10x.rna.seurat) = "seurat_clusters"
# Replot the projection with colored clusters and large labels
print( suppressWarnings( DimPlot( sc10x.rna.seurat, reduction = ifelse(exists("useReduction"), useReduction, "umap"), label=TRUE, label.size = 8)) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none", 
                plot.margin = margin(0, 0, 0, 0, "cm")));
cat(" \n \n"); # Required for '.tabset'


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - IDENTIFICATION
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes

Idents( sc10x.rna.seurat) <- "seurat_clusters"

# Identify marker genes
markers = FindAllMarkers( object          = sc10x.rna.seurat,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          verbose         = .VERBOSE);

#Print the datatable with all markers
print( datatable( markers[ which( markers$p_val_adj < FINDMARKERS_PVAL_THR), ]), caption = paste( "All markers by clusters with adj.pval <", FINDMARKERS_PVAL_THR))

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_logFC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP)) x else head( x, n = FINDMARKERS_SHOWTOP));
});

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkersHalf = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_logFC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP)) x else head( x, n = floor( FINDMARKERS_SHOWTOP /2)));
});

# Merge marker genes in a single data.frame and render it as datatable
topMarkersDF = do.call( rbind, topMarkers);
topMarkersHalfDF = do.call( rbind, topMarkersHalf);
# Select and order columns to be shown in datatable
topMarkersDT = topMarkersDF[c("gene", "cluster", "avg_logFC", "p_val_adj")]
topMarkersHalfDT = topMarkersHalfDF[c("gene", "cluster", "avg_logFC", "p_val_adj")]



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - TABLE
# ---------------------------------------------------------------
# ---------------------------------------------------------------
## @knitr heterogeneity_markerGenes_table

# Create datatable
datatable( topMarkersDT, 
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP), "All", paste("Top", FINDMARKERS_SHOWTOP)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list( 
                            list( # Center all columns except first one
                              targets = 1:(ncol( topMarkersDT)-1),
                              className = 'dt-center'),
                            list( # Set renderer function for 'float' type columns (LogFC)
                              targets = ncol( topMarkersDT)-2,
                              render = htmlwidgets::JS("function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}")),
                            list( # Set renderer function for 'scientific' type columns (PValue)
                              targets = ncol( topMarkersDT)-1,
                              render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toExponential(4);}"))), 
                          #fixedColumns = TRUE, # Does not work well with filter on this column
                          #fixedHeader = TRUE, # Does not work well with 'scrollX'
                          lengthMenu = list(c( 10, 50, 100, -1),
                                            c( 10, 50, 100, "All")),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          processing = TRUE, 
                          #search.regex= TRUE, # Does not work well with 'search.smart'
                          search.smart = TRUE,
                          select = TRUE, # Enable ability to select rows
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE));



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - HEATMAP
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x.rna.seurat, features = topMarkersDF[["gene"]])));

# Get the matrix of expression and associated clusters from Seurat object 
expMat = as.matrix( GetAssayData( sc10x.rna.seurat));
clusterID = Idents( sc10x.rna.seurat);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkersDF[["gene"]];
topMarkersHalfGenes = topMarkersHalfDF[["gene"]];
clusterOrdering = order( clusterID);

topMarkersGenes_expMat = expMat[topMarkersGenes, clusterOrdering];
topMarkersHalfGenes_expMat = expMat[topMarkersHalfGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Create a matrix for text on mouse over events (customize 'text' that is supposed to show value for adding info)
hoverTextMarkers = topMarkersGenes_expMat;
hoverTextMarkers[] = paste( paste( "Cell ID:", sc10x.rna.seurat[["numID", drop = FALSE]][colnames( topMarkersGenes_expMat)[col( topMarkersGenes_expMat)],]),
                            paste( "Value:", topMarkersGenes_expMat),
                            paste( "Cluster:", sort( clusterID)[col( topMarkersGenes_expMat)]),
                            paste( "Markers cl.:", topMarkersDF[["cluster"]][row( topMarkersGenes_expMat)]), 
                            sep = "<br>");


# Prepare rows and columns annotation bars (cluster groups)
rowsAnnot = data.frame(Cluster = topMarkersDF[["cluster"]]);
rowsAnnotHalf = data.frame(Cluster = topMarkersHalfDF[["cluster"]]);
colsAnnot = data.frame(Cluster = clusterID);

# Create heatmap
heatmapPlot = iheatmap(topMarkersGenes_expMat, 
                       row_labels = TRUE, 
                       col_labels = FALSE, 
                       text = hoverTextMarkers,
                       tooltip = setup_tooltip_options(row = TRUE,
                                                       col = TRUE,
                                                       value = TRUE,
                                                       prepend_row = "Gene: ",
                                                       prepend_col = "Cell: ",
                                                       prepend_value = ""),
                       row_annotation = rowsAnnot, 
                       col_annotation = colsAnnot,
                       row_annotation_colors = list(Cluster = hue_pal()(nlevels(Idents( sc10x.rna.seurat)))),
                       col_annotation_colors = list(Cluster = hue_pal()(nlevels(Idents( sc10x.rna.seurat)))));

# Hack to customize modebar buttons as for plotly native objects
# see: https://github.com/ropensci/iheatmapr/issues/38#issuecomment-445428166
heatmapPlot_widget = iheatmapr::to_widget(heatmapPlot);
# Edit the htmlwidget object itself
heatmapPlot_widget[["x"]][["config"]][["displaylogo"]] = FALSE;
heatmapPlot_widget[["x"]][["config"]][["toImageButtonOptions"]] = list( format='svg');  # This one does not seem to work... TODO: Check issue on GitHub
heatmapPlot_widget[["x"]][["config"]][["modeBarButtons"]]       = list( list( 'toImage'),
                                                                        list( 'zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d'));

# Render the heatmap
heatmapPlot_widget;


# Create a heatmap of Half Top Marker Genes to export to PDF

# -- Create colors for expression and annotations
col_expression <- brewer.pal(9, "Reds")
pal = colorRampPalette( col_expression)
col_cluster = hue_pal()(nlevels(Idents( sc10x.rna.seurat)))
names( col_cluster) = levels( Idents( sc10x.rna.seurat))

# -- Create the heatmap in PDF
pdf( file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "heatmap_halfTopMarkersGenes.pdf")),
     height = 9)

  # -- Create the column and rows annotations
  top_anno = HeatmapAnnotation( Cluster = colsAnnot$Cluster, 
                                col = list( Cluster = col_cluster),
                                show_legend = c( "Cluster" = FALSE))
  
  row_anno = rowAnnotation( Cluster = rowsAnnotHalf$Cluster,
                            col = list( Cluster = col_cluster),
                            show_legend = c( "Cluster" = FALSE))

  ComplexHeatmap::Heatmap( topMarkersHalfGenes_expMat,
                           col = pal( 20),
                           top_annotation = top_anno,
                           right_annotation = row_anno,
                           heatmap_legend_param = list( title = "Expression",  title_position = "leftcenter-rot"),
                           row_names_gp = gpar(fontsize = 4),
                           show_column_names = FALSE,
                           cluster_rows = FALSE,
                           cluster_columns = FALSE,
                           )
dev.off()

# -- Create the heatmap in PDF
pdf( file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "heatmap_halfTopMarkersGenes_horizontal.pdf")),
     width = 9)

  top_anno = HeatmapAnnotation( Cluster = rowsAnnotHalf$Cluster, 
                                col = list( Cluster = col_cluster),
                                show_legend = c( "Cluster" = FALSE))
  
  row_anno = rowAnnotation( Cluster = colsAnnot$Cluster,
                            col = list( Cluster = col_cluster),
                            show_legend = c( "Cluster" = FALSE))
  
  ComplexHeatmap::Heatmap( t( topMarkersHalfGenes_expMat),
                           col = pal( 20),
                           top_annotation = top_anno,
                           right_annotation = row_anno,
                           heatmap_legend_param = list( title = "Expression",  title_position = "leftcenter-rot"),
                           column_names_gp = gpar(fontsize = 4),
                           show_row_names = FALSE,
                           cluster_rows = FALSE,
                           cluster_columns = FALSE,
)
dev.off()


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - EXPRESSION on T-SNE and UMAP
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes_expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of marker genes (TODO: message if list empty)
resNULL=lapply(names(topMarkers), function(clusterName)
{
  cat("###### Cl.", clusterName, " \n");
  
  topMarkersGenes = topMarkers[[clusterName]][["gene"]]
  
  # Extract cells from current cluster to highlight them on a dimreduc plot
  clusterCells = which( Idents( sc10x.rna.seurat) == clusterName);
  # Create a dimreduc plot with current cluster highlighted (useful for large number of clusters)
  print( DimPlot( sc10x.rna.seurat, 
                  reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
                  cols="#44444422", 
                  cells.highlight = clusterCells, 
                  cols.highlight = "#FF000088", 
                  sizes.highlight = 1.5, 
                  order = clusterCells,  # Plot highlighted cells last
                  group.by=NULL) + 
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none") +
           ggtitle( paste( "Cluster", clusterName)));
         
  # Prepare panels as list of ggplots and remove axes title
  markersPanels = lapply( topMarkersGenes, function(featuresName) 
    { 
      FeaturePlot( sc10x.rna.seurat, features = featuresName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE) +
        theme( axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none");
    
    });
  
  resNULL = lapply( markersPanels, print);
  cat(" \n \n"); # Required for '.tabset'

  
});

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - EXPRESSION ON VIOLIN PLOTS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes_expression_violinplot

topMarkersGenes= sort( unique( do.call( c, lapply(names(topMarkers), function( clusterName){
  return( topMarkers[[clusterName]][["gene"]])
}))))

# Gather data to be visualized together (cell name + numID + metrics + Cluster)
cellsDataMarkers = cbind( "Cell" = colnames( sc10x.rna.seurat), # cell names from rownames conflict with search field in datatable
                   sc10x.rna.seurat[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = as.numeric(format(sc10x.rna.seurat[["percent.mito", drop = TRUE]], digits = 5)),
                   "Cluster" = Idents( sc10x.rna.seurat));

# # Add gene expression information on cell data
original_names = names( cellsDataMarkers)
for( marker_gene in topMarkersGenes){
  cellsDataMarkers = cbind( cellsDataMarkers, sc10x.rna.seurat@assays$RNA@data[ marker_gene, colnames( sc10x.rna.seurat)])
}
names( cellsDataMarkers) = c( original_names, topMarkersGenes)

# Generate plotly violin/jitter panels for each marker gene
for( marker_gene in topMarkersGenes){
  
  # Compute the maximal expression value of the marker genen over the clusters
  max_expression = max( cellsDataMarkers[ , marker_gene])
  
  # Compute the number of cells with positive expression in the different clusters
  positive_cellsDataMarkers = cellsDataMarkers[ which( cellsDataMarkers[ , marker_gene] > 0), "Cluster"]
  positive_cellsDataMarkers = as.data.frame( table( positive_cellsDataMarkers))
  names( positive_cellsDataMarkers) = c( "Cluster", "NbCells")
  
  # Compute the number of cells with null expression in the different clusters
  null_cellsDataMarkers = cellsDataMarkers[ which( cellsDataMarkers[ , marker_gene] == 0), "Cluster"]
  null_cellsDataMarkers = as.data.frame( table( null_cellsDataMarkers))
  names( null_cellsDataMarkers) = c( "Cluster", "NbCells")
  
  # Plot the expression of the marker gene for each cluster separatly as violin + jitter
  # and adding the number of positive expressed cells and null expressed cells above each cluster plot
   print( 
     ggplot( cellsDataMarkers, aes_string ( x="Cluster", y=as.name( marker_gene))) +
      geom_violin( aes( col = Cluster, fill=Cluster)) +
      geom_jitter( width=0.2, size = 0.2, shape = 1) + 
      geom_text( data = positive_cellsDataMarkers, aes( x=Cluster, y = 1.15*max_expression, label=NbCells), size = 3) +
      geom_text( data = null_cellsDataMarkers, aes( x=Cluster, y = 1.1*max_expression, label=NbCells), size = 3) +
      theme_minimal() +
      theme( legend.position = "None")
   )
}



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MONITORED GENES
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_monitoredGenes

# Just remind the warning for genes names not in object
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));




## @knitr heterogeneity_monitoredGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x.rna.seurat, features = unlist(MONITORED_GENES))));

# Get the matrix of expression and associated clusters from Seurat object 
expMat = as.matrix( GetAssayData( sc10x.rna.seurat));
clusterID = Idents( sc10x.rna.seurat);

# Select monitored genes and reorder cells to group clusters together
monitoredGenes = unlist( MONITORED_GENES);
clusterOrdering = order( clusterID);

expMat = expMat[monitoredGenes, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Create a matrix for text on mouse over events (customize 'text' that is supposed to show value for adding info)
hoverTextMarkers = expMat;
hoverTextMarkers[] = paste( paste( "Cell ID:", sc10x.rna.seurat[["numID", drop = FALSE]][colnames( expMat)[col( expMat)],]),
                            paste( "Value:", expMat),
                            paste( "Cluster:", sort( clusterID)[col( expMat)]),
                            paste( "Monitored Grp.:", rep(names(MONITORED_GENES), sapply(MONITORED_GENES, length))[row( expMat)]), 
                            sep = "<br>");

# Prepare rows and columns annotation bars (monitored group and cluster respectively)
rowsAnnot = data.frame( Monitored = rep( names( MONITORED_GENES), sapply( MONITORED_GENES, length)));
colsAnnot = data.frame( Cluster = clusterID);

# Create heatmap
heatmapPlot = iheatmap( expMat,
                        row_labels = TRUE,
                        col_labels = FALSE,
                        text = hoverTextMarkers,
                        tooltip = setup_tooltip_options(row = TRUE,
                                                        col = TRUE,
                                                        value = TRUE,
                                                        prepend_row = "Gene: ",
                                                        prepend_col = "Cell: ",
                                                        prepend_value = ""),
                        row_annotation = rowsAnnot,
                        col_annotation = colsAnnot,
                        col_annotation_colors = list( Cluster = hue_pal()( nlevels( Idents( sc10x.rna.seurat)))));

# Hack to customize modebar buttons as for plotly native objects
# see: https://github.com/ropensci/iheatmapr/issues/38#issuecomment-445428166
heatmapPlot_widget = iheatmapr::to_widget(heatmapPlot);
# Edit the htmlwidget object itself
heatmapPlot_widget[["x"]][["config"]][["displaylogo"]] = FALSE;
heatmapPlot_widget[["x"]][["config"]][["toImageButtonOptions"]] = list( format='svg');  # This one does not seem to work... TODO: Check issue on GitHub
heatmapPlot_widget[["x"]][["config"]][["modeBarButtons"]]       = list( list( 'toImage'),
                                                                        list( 'zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d'));

# Render the heatmap
heatmapPlot_widget;


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MONITORED GENES - EXPRESSION - on T-SNE AND UMAP
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_monitoredGenes_expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of monitored genes (TODO: message if list empty)
resNULL = lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("######", monitoredGroup, " \n");

  # List of plots (or error message if feature not found)
  expressionPlots = lapply( MONITORED_GENES[[monitoredGroup]], function(featuresName)
  {
    tryCatch( FeaturePlot( sc10x.rna.seurat, features = featuresName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
                theme( axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       legend.position = "none"),
              error=function(e){return(conditionMessage(e))}); # In case gene name could not be found in object (should have been removed earlier anyway)
  });
  
  resNULL = lapply( expressionPlots, print);
  cat(" \n \n"); # Required for '.tabset'
});


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MONITORED GENES - EXPRESSION ON VIOLIN PLOTS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_monitoredGenes_expression_violinplot

# Gather data to be visualized together (cell name + numID + metrics + Cluster)
cellsDataMonitored = cbind( "Cell" = colnames( sc10x.rna.seurat), # cell names from rownames conflict with search field in datatable
                          sc10x.rna.seurat[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                          "percent.mito" = as.numeric(format(sc10x.rna.seurat[["percent.mito", drop = TRUE]], digits = 5)),
                          "Cluster" = Idents( sc10x.rna.seurat));

all_monitored_genes = sort( unique( unlist( MONITORED_GENES, use.names = FALSE)))

# # Add gene expression information on cell data
original_names = names( cellsDataMonitored)
for( monitored_gene in all_monitored_genes){
  cellsDataMonitored = cbind( cellsDataMonitored, sc10x.rna.seurat@assays$RNA@data[ monitored_gene, colnames( sc10x.rna.seurat)])
}
names( cellsDataMonitored) = c( original_names, all_monitored_genes)

# Generate plotly violin/jitter panels for each monitored gene
for( monitored_gene in all_monitored_genes){
  
  # Compute the maximal expression value of the monitored genen over the clusters
  max_expression = max( cellsDataMonitored[ , monitored_gene])
  
  # Compute the number of cells with positive expression in the different clusters
  positive_cellsDataMonitored = cellsDataMonitored[ which( cellsDataMonitored[ , monitored_gene] > 0), "Cluster"]
  positive_cellsDataMonitored = as.data.frame( table( positive_cellsDataMonitored))
  names( positive_cellsDataMonitored) = c( "Cluster", "NbCells")
  
  # Compute the number of cells with null expression in the different clusters
  null_cellsDataMonitored = cellsDataMonitored[ which( cellsDataMonitored[ , monitored_gene] == 0), "Cluster"]
  null_cellsDataMonitored = as.data.frame( table( null_cellsDataMonitored))
  names( null_cellsDataMonitored) = c( "Cluster", "NbCells")
  
  # Plot the expression of the monitored gene for each cluster separatly as violin + jitter
  # and adding the number of positive expressed cells and null expressed cells above each cluster plot
  print( 
    ggplot( cellsDataMonitored, aes_string ( x="Cluster", y=as.name( monitored_gene))) +
      geom_violin( aes( col = Cluster, fill=Cluster)) +
      geom_jitter( width=0.2, size = 0.2, shape = 1) + 
      geom_text( data = positive_cellsDataMonitored, aes( x=Cluster, y = 1.15*max_expression, label=NbCells), size = 3) +
      geom_text( data = null_cellsDataMonitored, aes( x=Cluster, y = 1.1*max_expression, label=NbCells), size = 3) +
      theme_minimal() +
      theme( legend.position = "None")
  )
}

  
  