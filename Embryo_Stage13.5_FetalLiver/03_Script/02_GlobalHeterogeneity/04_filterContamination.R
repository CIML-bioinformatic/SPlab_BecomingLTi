# ##############################################################################################
# In this script, we identify the clusters expressing the most some genes expressed by cells
# considered as contamination. Cells in the corresponding cluster are removed
# ##############################################################################################


## @knitr filtering_contamination

# Remove cluster of cells with high expression of genes characteristic of contaminating cells
# ............................................................................................

# Build a data.frame to analyze the distribution of contamination genes along the clusters
SC10X_QCFILTERED_NORMALIZED_DATA_DF = as.data.frame(t(as.matrix( GetAssayData( sc10x.rna.seurat , assay = "RNA"))))
SC10X_QCFILTERED_NORMALIZED_DATA_DF = SC10X_QCFILTERED_NORMALIZED_DATA_DF[ , CONTAMINATION_GENES]
SC10X_QCFILTERED_NORMALIZED_DATA_DF$UniqueCellID = rownames(SC10X_QCFILTERED_NORMALIZED_DATA_DF)
Idents( sc10x.rna.seurat) = "seurat_clusters"
SC10X_QCFILTERED_NORMALIZED_DATA_DF$SeuratClusters = Idents( sc10x.rna.seurat)[ SC10X_QCFILTERED_NORMALIZED_DATA_DF$UniqueCellID]

all_cells_to_filter = vector()
all_cells_to_filter_inside_cluster = vector()
for( characteristic_gene in CONTAMINATION_GENES){
  
  cat("<H4>Contamination indicated by gene", characteristic_gene, "</H4>")
  
  # Show the expression of the characteristic_gene on a UMAP embedding
  print( FeaturePlot( sc10x.rna.seurat, features = characteristic_gene, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE) +
               theme( axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "none")
  )
  
  # -- Look at the cells expressing characteristic_gene
  contamination_df = SC10X_QCFILTERED_NORMALIZED_DATA_DF[ ,c( characteristic_gene, "SeuratClusters")]
  contamination_df[, characteristic_gene] = SC10X_QCFILTERED_NORMALIZED_DATA_DF[ , characteristic_gene] > 0
  contamination_barcodes = row.names( contamination_df)[ contamination_df[ , characteristic_gene]]
  
  # -- Compute the percentage of cells expressing characteristic_gene in each cluster
  cluster_distribution_table = table( contamination_df[ c( "SeuratClusters", characteristic_gene)])
  cluster_contamination_percentage_df = data.frame( ClusterID = rownames( cluster_distribution_table), 
                                                    ClusterSize = cluster_distribution_table[ , "TRUE"] + cluster_distribution_table[ , "FALSE"],
                                                    ExpressingCells = cluster_distribution_table[ , "TRUE"],
                                                    Percentage = signif( 100*cluster_distribution_table[ , "TRUE"] / (cluster_distribution_table[ , "TRUE"] + cluster_distribution_table[ , "FALSE"]), 3)
  )
  print( htmltools::tagList( list( 
    table = datatable( cluster_contamination_percentage_df,
                       caption = paste( "Percentage of cells expressing", characteristic_gene, "in clusters"),
                       rownames = FALSE)
  )))
  
  cat("<BR>Number of cells expressing", characteristic_gene, ":", length( contamination_barcodes))
  
  # -- Identify cluster to filter due to contamination
    contamination_cluster_set = cluster_contamination_percentage_df[ which( cluster_contamination_percentage_df$Percentage > 60), "ClusterID"]
  cat("<BR>ClusterID to remove due to high content of cells expressing", characteristic_gene, ":", paste( contamination_cluster_set, collapse = ","))
  
  # -- Identify all the cells to remove due to association to contaminations
  cells_to_filter_in_cluster = row.names( contamination_df)[ contamination_df$SeuratClusters %in% contamination_cluster_set]
  cells_to_filter = unique( c( contamination_barcodes, cells_to_filter_in_cluster))
  
  cat("<BR>Number of cells identified as contamination in filtered clusters:", length( cells_to_filter_in_cluster))                          
  cat("<BR>Total number of cells identified as contamination:", length( cells_to_filter))
  
  # -- Accumulate the cells to filter out
  all_cells_to_filter_inside_cluster = unique( c( all_cells_to_filter_inside_cluster, cells_to_filter_in_cluster))
  all_cells_to_filter = unique( append( all_cells_to_filter, cells_to_filter))
  
  # -- Export the identified cells to file
  excluded_cells_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "comtamination_excluded_cells_", characteristic_gene, ".txt"))
  write.table( data.frame( cellid = cells_to_filter), file = excluded_cells_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
  cat("<BR><BR>Writing barcodes of excluded cells by ", characteristic_gene, "expression to file: <BR>", excluded_cells_outfile)
}

# Filter the excluded cells in the Seurat object
cat("<H4>Result of contamination analysis</H4>")
cat("<BR><b>Number of cells to filter out:", length( all_cells_to_filter),"</b>")
cat("<BR><b>Number of cells to filter out that are never in filtered clusters:", length( setdiff( all_cells_to_filter, all_cells_to_filter_inside_cluster)),"</b>")
cat("<BR><b>Number of cells before filtering of contamination:", length( Cells( sc10x.rna.seurat)),"</b>")
sc10x.rna.seurat.nocontamination <- subset( x = sc10x.rna.seurat, 
                                         cells = setdiff ( Cells( sc10x.rna.seurat), all_cells_to_filter)
                                        )
cat("<BR><b>Number of cells after filtering of contamination:", length( Cells( sc10x.rna.seurat.nocontamination)), "</b>")

filtered_cells_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "contamination_excluded_all_cells.txt"))
write.table( data.frame( cellid = all_cells_to_filter), file = filtered_cells_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of all excluded cells to file: <BR>", filtered_cells_outfile)

cat("<BR>")

# SAVE FILTERED SEURAT OBJECT TO RDATA
# ---------------------------

## @knitr saveSeuratRObject

saveRDS( sc10x.rna.seurat.nocontamination, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sc10x.rna.seurat.nocontamination.RDS")))
