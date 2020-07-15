# ##############################################################################
# This script aims to produce analysis to caracterise clusters fourn in
# previous step
# ##############################################################################


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PAIRWISE MARKER GENES - IDENTIFICATION
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr characterizeCluster_pairwiseMarkers

Idents( sc10x.rna.seurat) <- "seurat_clusters"

# Create the pairwise comparison by combinations of clusters
seurat_clusters <- sort( unique( Idents( sc10x.rna.seurat)))
pairwise <- combn( seurat_clusters, 2)

# Accumulate the comparisons results in a single dataframe
pairwise_marker_genes_df = data.frame()
for(i in 1:ncol(pairwise)) {
  # find the markers genes if first direction
  pairwise_marker_genes = FindMarkers( sc10x.rna.seurat, 
                                       ident.1 = pairwise[1, i], 
                                       ident.2 = pairwise[2, i],
                                       test.use        = FINDMARKERS_METHOD,
                                       only.pos        = FINDMARKERS_ONLYPOS,
                                       min.pct         = FINDMARKERS_MINPCT,
                                       logfc.threshold = FINDMARKERS_LOGFC_THR,
                                       verbose         = .VERBOSE);
  
  # find the markers genes if second direction
  pairwise_marker_genes_switch = FindMarkers( sc10x.rna.seurat, 
                                              ident.1 = pairwise[2, i], 
                                              ident.2 = pairwise[1, i],
                                              test.use        = FINDMARKERS_METHOD,
                                              only.pos        = FINDMARKERS_ONLYPOS,
                                              min.pct         = FINDMARKERS_MINPCT,
                                              logfc.threshold = FINDMARKERS_LOGFC_THR,
                                              verbose         = .VERBOSE);
  
  # Eliminate the non-significant genes and keep only 100 best
  pairwise_marker_genes = pairwise_marker_genes[ which( pairwise_marker_genes$p_val_adj < FINDMARKERS_PVAL_THR), ]
  pairwise_marker_genes = pairwise_marker_genes[ order( pairwise_marker_genes[ ,"p_val_adj"], decreasing = FALSE), ]
  #pairwise_marker_genes = head( pairwise_marker_genes, 100)
  
  pairwise_marker_genes_switch = pairwise_marker_genes_switch[ which( pairwise_marker_genes_switch$p_val_adj < FINDMARKERS_PVAL_THR), ]
  pairwise_marker_genes_switch = pairwise_marker_genes_switch[ order( pairwise_marker_genes_switch[ ,"p_val_adj"], decreasing = FALSE), ]
  #pairwise_marker_genes_switch = head( pairwise_marker_genes_switch, 100)
  
  # Eliminate row names to allow accumulations in a single dataframe
  pairwise_marker_genes$gene.name = row.names( pairwise_marker_genes)
  row.names( pairwise_marker_genes) = seq( 1, nrow( pairwise_marker_genes), 1)
  pairwise_marker_genes_switch$gene.name = row.names( pairwise_marker_genes_switch)
  row.names( pairwise_marker_genes_switch) = seq( 1, nrow( pairwise_marker_genes_switch), 1)
  
  # Add the information about compared clusters
  pairwise_marker_genes$cluster1 = pairwise[1, i]
  pairwise_marker_genes$cluster2 = pairwise[2, i]
  pairwise_marker_genes_switch$cluster1 = pairwise[2, i]
  pairwise_marker_genes_switch$cluster2 = pairwise[1, i]
  
  # Acculumulate the results
  pairwise_marker_genes_df = rbind( pairwise_marker_genes_df, pairwise_marker_genes)
  pairwise_marker_genes_df = rbind( pairwise_marker_genes_df, pairwise_marker_genes_switch)
}

pairwise_marker_genes_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "pairwise_marker_genes.tsv"))
write.table( pairwise_marker_genes_df, 
             file = pairwise_marker_genes_file_path,
             quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

# Summarize information per cluster and per gene
pairwise_topMarkers = data.frame()
pairwise_markers_summary_df = data.frame()
for( current_cluster in seurat_clusters){
  current_df = pairwise_marker_genes_df[ which( pairwise_marker_genes_df$cluster1 == current_cluster), ]
  new_df = data.frame()
  for( current_gene in unique( current_df$gene.name)){
    marked_cluster_set = as.character( sort( current_df[ which( current_df$gene.name == current_gene), "cluster2"]), decreasing = FALSE)
    new_df = rbind( new_df, data.frame( gene.name = current_gene,
                                        cluster1 = current_cluster,
                                        marked.clusters = paste( marked_cluster_set, collapse = ", "),
                                        marked.clusters.nb = length( marked_cluster_set)))
  }
  new_df = new_df[ order( new_df[ ,"marked.clusters.nb"], decreasing = TRUE), ]
  pairwise_markers_summary_df = rbind( pairwise_markers_summary_df, new_df)
  pairwise_topMarkers = rbind( pairwise_topMarkers, head( new_df, n = FINDMARKERS_SHOWTOP))
}

# Export the summary table to file
pairwise_markers_summary_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "pairwise_marker_genes_summary.tsv"))
write.table( pairwise_markers_summary_df, 
             file = pairwise_markers_summary_file_path,
             quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

# Export the topMarkers summary table to file
pairwise_markers_summary_topMarkers_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "pairwise_marker_genes_summary_topMarkers.tsv"))
write.table( pairwise_topMarkers, 
             file = pairwise_markers_summary_topMarkers_file_path,
             quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PAIRWISE MARKER GENES - CLUSTERS SIMILARITY
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr characterizeCluster_similarity

# Get all the genes that are markers for pairwise comparison
all_deg_genes = unique( pairwise_marker_genes_df[ , "gene.name"])

# Build a matrix of similatiry, showing the complement of the ratio of the markers genes between each pair
# of cluster divided by the total number of pairwise marker genes
comparison_df = data.frame()
for( cluster1 in seurat_clusters){
  
  similarities_set = vector()
  for( cluster2 in seurat_clusters){
    up_genes_vs_cluster2 = pairwise_marker_genes_df[ which( pairwise_marker_genes_df$cluster1 == cluster1
                                                            & pairwise_marker_genes_df$cluster2 == cluster2), ]
    
    similarity = 100*( 1 - (length( unique( up_genes_vs_cluster2$gene.name)) / length( all_deg_genes)))
    similarities_set = append( similarities_set, similarity)
  }
  comparison_df = rbind( comparison_df, similarities_set)
}
names( comparison_df) = seurat_clusters
row.names( comparison_df) = seurat_clusters

# Plot the similarity matrix as heatmap
print( pheatmap( comparison_df, cluster_cols = FALSE, cluster_rows = FALSE))


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PAIRWISE MARKER GENES - FUNCTIONAL ENRICHMENT
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr characterizeCluster_functionalEnrichment

# Define the background
# -- get the complete list of genes from the seurat object
#background = row.names( sc10x.rna.seurat)
background = unique(pairwise_marker_genes_df$gene.name)

# -- associate ENTREZ IDs to gene SYMBOL
eg_all = bitr(background, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db",  drop = FALSE)

# -- remove potential ENTREZ IDs duplicated
eg_all = eg_all[!duplicated(eg_all["SYMBOL"]),]
eg_all = na.omit( eg_all)

# Parse the clusters to look for enrichment analysis on each of them
external_clusters = c( 6,8)
for( cluster in sort( unique( unlist( sc10x.rna.seurat[[ "seurat_clusters"]])))){

  # print( htmltools::tagList( tags$h5( paste( "Enriched Biological Processes for cluster", cluster))))
  # Define the list of genes to test
  # -- get the genes associated as marker for the cluster
  if( cluster %in% external_clusters){
    genes = unique( pairwise_marker_genes_df[ pairwise_marker_genes_df$cluster1 == cluster , "gene.name"])
  }else{
    genes = unique( pairwise_marker_genes_df[ pairwise_marker_genes_df$cluster1 == cluster &
                                              !pairwise_marker_genes_df$cluster2 %in% external_clusters , "gene.name"])
  }
  
  # -- associate ENTREZ IDs to gene SYMBOL
  eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db",  drop = FALSE)
  
  # remove potential ENTREZ IDs duplicated
  eg = eg[ !duplicated(eg[ "SYMBOL"]),]
  eg = na.omit( eg)
  
  # perform GO enrichment (Biological Process)
  ego_BP <- enrichGO(gene          = eg$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     universe = eg_all$ENTREZID,
                     readable      = TRUE)
  
  # Convert result to dataframe
  ego_BP <- data.frame(ego_BP)
  
  if( nrow( ego_BP) > 0){
    # Add a column with the ratio of genes found in the geneset
    ego_BP$Fold <- sapply(ego_BP$GeneRatio, function(x) eval(parse(text=x)))/sapply(ego_BP$BgRatio, function(x) eval(parse(text=x)))
    
    # Shorten the precision of numerical scores
    ego_BP[,c(5,6,7,10)] <- signif(ego_BP[,c(5,6,7,10)],digits = 3)
    
    # Removes row names of dataframe
    rownames(ego_BP) <- NULL
    
    # Print the datatable in report
    print( htmltools::tagList( list( tags$h5( paste( "Enriched Biological Processes for cluster", cluster)),
                                  datatable( ego_BP[ , which( names( ego_BP) != "geneID")], caption = paste( "Enriched BP for cluster", cluster)))))
    
    # Export the result to file
    write.table(ego_BP, file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "_EnrichBP_Cluster", cluster, ".csv")),
                                   sep=",", col.names = T, row.names = FALSE, quote=FALSE)
  }else{
    cat("<BR>No enrichment found for this cluster")
  }
  
  cat("<HR>")

}



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PAIRWISE MARKER GENES - FUNCTIONAL ENRICHMENT FROM NEIGHBORS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr characterizeCluster_functionalEnrichment_fromNeighbors
neighbor_clusters = list( "0" = c( 1,2,7,10,11),
                          "1" = c( 0,2,4,7,9),
                          "2" = c( 0, 1, 7), 
                          "3" = c( 5, 7, 11, 12), 
                          "4" = c( 1,6,7,9,12), 
                          "5" = c( 3,10,11), 
                          "6" = c( 3,4,8), 
                          "7" = c( 0,1,3,4,10,11,12), 
                          "8" = c( 3,6), 
                          "9" = c( 1,4), 
                          "10" = c( 0,11,5), 
                          "11" = c( 0,3,5,7,10), 
                          "12" = c( 3,4,7))

# Parse the clusters to look for enrichment analysis on each of them
external_clusters = c( 6,8)
for( cluster in sort( unique( unlist( sc10x.rna.seurat[[ "seurat_clusters"]])))){
  
  background_clusters = c( cluster, neighbor_clusters[[ paste0( "", cluster)]])
  
  # Define the background
  # -- get the complete list of genes from the seurat object
  #background = row.names( sc10x.rna.seurat)
  background = unique(pairwise_marker_genes_df[ which( pairwise_marker_genes_df$cluster1 %in% background_clusters), "gene.name"])
  
  # -- associate ENTREZ IDs to gene SYMBOL
  eg_all = bitr(background, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db",  drop = FALSE)
  
  # -- remove potential ENTREZ IDs duplicated
  eg_all = eg_all[!duplicated(eg_all["SYMBOL"]),]
  eg_all = na.omit( eg_all)
  
  # print( htmltools::tagList( tags$h5( paste( "Enriched Biological Processes for cluster", cluster, "against neighbor clusters"))))
  # Define the list of genes to test
  # -- get the genes associated as marker for the cluster
  if( cluster %in% external_clusters){
    genes = unique( pairwise_marker_genes_df[ pairwise_marker_genes_df$cluster1 == cluster , "gene.name"])
  }else{
    genes = unique( pairwise_marker_genes_df[ pairwise_marker_genes_df$cluster1 == cluster &
                                                !pairwise_marker_genes_df$cluster2 %in% external_clusters , "gene.name"])
  }
  
  # -- associate ENTREZ IDs to gene SYMBOL
  eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db",  drop = FALSE)
  
  # remove potential ENTREZ IDs duplicated
  eg = eg[ !duplicated(eg[ "SYMBOL"]),]
  eg = na.omit( eg)
  
  # perform GO enrichment (Biological Process)
  ego_BP <- enrichGO(gene          = eg$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     universe = eg_all$ENTREZID,
                     readable      = TRUE)
  
  # Convert result to dataframe
  ego_BP <- data.frame(ego_BP)
  
  if( nrow( ego_BP) > 0){
    # Add a column with the ratio of genes found in the geneset
    ego_BP$Fold <- sapply(ego_BP$GeneRatio, function(x) eval(parse(text=x)))/sapply(ego_BP$BgRatio, function(x) eval(parse(text=x)))
    
    # Shorten the precision of numerical scores
    ego_BP[,c(5,6,7,10)] <- signif(ego_BP[,c(5,6,7,10)],digits = 3)
    
    # Removes row names of dataframe
    rownames(ego_BP) <- NULL
    
    # Print the datatable in report
    print( htmltools::tagList( list( tags$h5( paste( "Enriched Biological Processes for cluster", cluster, "against neighbor clusters")),
                                     datatable( ego_BP[ , which( names( ego_BP) != "geneID")], caption = paste( "Enriched BP for cluster", cluster)))))
    
    # Export the result to file
    write.table(ego_BP, file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "_EnrichBP_Cluster", cluster, "_fromNeighbors.csv")),
                sep=",", col.names = T, row.names = FALSE, quote=FALSE)
  }else{
    cat("<BR>No enrichment found for this cluster")
  }
  
  cat("<HR>")
}
