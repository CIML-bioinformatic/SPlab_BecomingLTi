# ###############################################################
# This script aims to use Monocle method to classify cells
# with or without using predefined markers.
# ###############################################################

# This script is inpired by the Monocle tutorial:
# http://cole-trapnell-lab.github.io/monocle-release/docs/


# ###############################
# Classify cells by type
# ###############################

## @knitr classify_cells_by_type
set.seed( 123)

#cat("<BR><BR> Using cth_DIFFERENTIAL as cell Type hierarchy")
#SELECTED_HIERARCHY = cth_DIFFERENTIAL

cat("<BR><BR> Using cth_ORIGINAL as cell Type hierarchy")
SELECTED_HIERARCHY = cth_ORIGINAL

# Apply the cell type hierarchy to classify cells
filtered_cds_cell_type <- classifyCells( filtered_cds, SELECTED_HIERARCHY)

# Look at the individual signature to classify cells
filtered_cds_cell_type_list = list()
cell_name_type_list = list()
for( type_name in names( cth_type_list)){
  filtered_cds_cell_type_list[[ type_name]] = classifyCells( filtered_cds, cth_type_list[[ type_name]])
  cell_name_type_list[[ type_name]] = row.names( filtered_cds_cell_type_list[[ type_name]]@phenoData@data)[ which( filtered_cds_cell_type_list[[ type_name]]@phenoData@data$CellType != "Unknown")]
}

# Look at the intersection of individual signatures to locate the ambiguous classification
done_names = vector()
count_matrix = matrix( rep( -1, length( names( cth_type_list))* length( names( cth_type_list))), ncol= length( names( cth_type_list)))
colnames( count_matrix) = names( cth_type_list)
rownames( count_matrix) = names( cth_type_list)
for( type_name_1 in names( cth_type_list)){
  done_names = append( done_names, type_name_1)
  count_matrix[ type_name_1, type_name_1] = length( cell_name_type_list[[ type_name_1]])
  for( type_name_2 in names( cth_type_list)){
    if( !type_name_2 %in% done_names){
      count_matrix[ type_name_1, type_name_2] = length( intersect( cell_name_type_list[[ type_name_1]], cell_name_type_list[[ type_name_2]]))
      count_matrix[ type_name_2, type_name_1] = count_matrix[ type_name_1, type_name_2]
    }
  }
}
print( datatable( count_matrix, caption = "Number of cells identified by individual signatures or pair of signatures"))

# Plot the proportion of classified cells
cat("<HR>")
pie <- ggplot(pData( filtered_cds_cell_type), aes(x = factor(1), fill = factor( CellType))) + geom_bar( width = 1) +
      coord_polar( theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
print( pie)

# Plot the classified cells in the t-SNE map computed in the previous analysis step
TSNE_COORD_DF$CellType = pData( filtered_cds_cell_type)[ row.names( TSNE_COORD_DF), "CellType"]
ggplot( TSNE_COORD_DF) +
  geom_point( aes( x = tSNE_1, y = tSNE_2, col=CellType))
ggplot( TSNE_COORD_DF) +
  geom_point( aes( x = tSNE_1, y = tSNE_2, col=CellType)) +
  facet_wrap( ~ CellType)

# Plot the classified cells in the t-UMAP map computed in the previous analysis step
UMAP_COORD_DF$CellType = pData( filtered_cds_cell_type)[ row.names( UMAP_COORD_DF), "CellType"]
ggplot( UMAP_COORD_DF) +
  geom_point( aes( x = UMAP_1, y = UMAP_2, col=CellType))
ggplot( UMAP_COORD_DF) +
  geom_point( aes( x = UMAP_1, y = UMAP_2, col=CellType)) +
  facet_wrap( ~ CellType)

# Export t-SNE coordinates of cells with cell type
tsne_map_with_cell_type_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_Tsnecoord_WithCellType.tsv"))
cat("<BR>Exporting t-SNE map with cell types:<BR>", tsne_map_with_cell_type_file_path)
write.table( TSNE_COORD_DF, file = tsne_map_with_cell_type_file_path, quote = FALSE, col.names = TRUE, row.names = TRUE, sep=",")

# Export UMAP coordinates of cells with cell type
umap_map_with_cell_type_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "QCfiltered_Normalized_Contaminationfiltered_Umapcoord_WithCellType.tsv"))
cat("<BR>Exporting UMAP map with cell types:<BR>", umap_map_with_cell_type_file_path)
write.table( UMAP_COORD_DF, file = umap_map_with_cell_type_file_path, quote = FALSE, col.names = TRUE, row.names = TRUE, sep=",")

# ##########################################  
# Classify cells using markers
# ##########################################

## @knitr identify_markers_from_celltype
set.seed( 123)

# Select the genes that co-vary with the markers genes used to classify cells
marker_diff <- markerDiffTable( filtered_cds_cell_type[ EXPRESSED_GENES,],
                                SELECTED_HIERARCHY,
                                # residualModelFormulaStr = "~Media + num_genes_expressed",
                                cores = 1)

# Look at the most specific gene for each cell type
cell_type_candidate_clustering_genes <- row.names( subset( marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity( filtered_cds_cell_type[ cell_type_candidate_clustering_genes,], SELECTED_HIERARCHY)
print( datatable( selectTopMarkers( marker_spec, 10)))

# Keep the 500 most specific genes to cluster the cells
semisup_clustering_genes <- unique( selectTopMarkers( marker_spec, 200)$gene_id)
filtered_cds_cell_type <- setOrderingFilter( filtered_cds_cell_type, semisup_clustering_genes)
plot_ordering_genes( filtered_cds_cell_type)
plot_pc_variance_explained( filtered_cds_cell_type, return_all = F)

# Build the T-SNE map
set.seed( 123)
filtered_cds_cell_type <- reduceDimension( filtered_cds_cell_type, max_components = 2, num_dim = 6,
                        norm_method = 'log',
                        reduction_method = 'tSNE',
                        verbose = T)

# Find the cell clusters
#filtered_cds_cell_type <- clusterCells( filtered_cds_cell_type, method = "densityPeak", num_clusters = 9)
#filtered_cds_cell_type <- clusterCells( filtered_cds_cell_type, method = "louvain", k=200)
filtered_cds_cell_type <- clusterCells( filtered_cds_cell_type, method = "DDRTree", num_clusters = 5)

# Plot the data by cell types with one color per cell type
plot_cell_clusters( filtered_cds_cell_type, 1, 2, color = "CellType") +
  # scale_colour_manual(values = c("violet", "blue", "green", "red", "grey")) +
facet_wrap(~CellType)

# Plot the data by cell types with one color per cluster
plot_cell_clusters( filtered_cds_cell_type, 1, 2, color = "Cluster") +
  facet_wrap(~CellType)

# Plot the expression of markers in data per cell type
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "Id2")) + 
    scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
    ggtitle( paste( "Dispersion of", "Id2", "in identified populations"))
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "Rorc")) + 
  scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
  ggtitle( paste( "Dispersion of", "Rorc", "in identified populations"))
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "Zbtb16")) + 
  scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
  ggtitle( paste( "Dispersion of", "Zbtb16", "in identified populations"))
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "Flt3")) + 
  scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
  ggtitle( paste( "Dispersion of", "Flt3", "in identified populations"))


# ###############################
# Classify cells without markers
# ###############################

## @knitr classify_cells_without_marker

disp_table <- dispersionTable( filtered_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
FILTERED_CDS_PSEUDOTIME_UNSUPERVISED <- setOrderingFilter( filtered_cds, unsup_clustering_genes$gene_id)
plot_ordering_genes( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED)

plot_pc_variance_explained( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, return_all = F)

FILTERED_CDS_PSEUDOTIME_UNSUPERVISED <- reduceDimension( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, max_components = 2, num_dim = 6,
                                                         reduction_method = 'tSNE', verbose = T)

FILTERED_CDS_PSEUDOTIME_UNSUPERVISED <- clusterCells( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, num_clusters = 9)

cat("<H5>t-SNE map with clusters</H5>")
plot_cell_clusters( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, 1, 2, color = "Cluster")

for( signature_name in names( MONITORED_GENES)){
  cat("<H5>t-SNE map with genes expression of", signature_name, "</H5>")
  print( plot_cell_clusters( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, 1, 2, color = "Cluster", markers = MONITORED_GENES[[ signature_name]]))
}
