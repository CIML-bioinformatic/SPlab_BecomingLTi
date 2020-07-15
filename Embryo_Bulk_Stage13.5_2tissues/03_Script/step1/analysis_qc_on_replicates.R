# #####################################################
# This script aims to analyse the replicates of the
# various experiment and compare the replicates
# between each other in more global analysis
# #####################################################

# -------------------------------------------------------------
# Execute a pairs comparison on each condition to analyse 
# -------------------------------------------------------------

## @knitr qc_on_replicates_pairs

for( condition in EXPERIMENT_ID_SET){
  
  comparison_pairs = ggpairs( data = as.data.frame( log1pNORMALIZED_COUNTS)[, which( CONDITION_NAME_COLUMN_SET == condition)],
                              lower = list(continuous = wrap("smooth", alpha = 0.5, size=0.3)),
                              diag = list(continuous = wrap("densityDiag")),
                              title = paste( "Comparison of replicates of ", condition))
  print( comparison_pairs)
}

# -------------------------------------------------------------
# For each subpopulation, execute a PCA analysis
# -------------------------------------------------------------

## @knitr qc_on_replicates_heterogeneity_subpop

all_high_loadings = vector()
widget_list = list()
count_widget = 1
for( subtype in list( EP=c("EP1", "EP2", "EP3", "EP4"), 
                      FL=c("FLP1", "FLP2", "FLP3", "FLP4"), 
                      SUB1=c("EP3", "EP4", "FLP3", "FLP4"),
                      SUB2=c("EP2", "EP3", "FLP2", "FLP3"))){
  cat("<H3>Analyzing heterogeneity in populations ", subtype, "</H3>")
  # acp_df = as.data.frame( t(NORMALIZED_COUNTS))[ which( CONDITION_NAME_COLUMN_SET %in% subtype), ]
  
  # Get the lits of existing subtype from the listed ones
  sample_name_list = CONDITION_NAME_COLUMN_SET[ which( CONDITION_NAME_COLUMN_SET %in% subtype)]
  
  # Build a color palette from the list of subtypes 
  color_palette = brewer.pal( length( unique( sample_name_list)),"Set1")
  sample_colors = unique( unlist( lapply( sample_name_list, function( x){
    return( color_palette[ which( unique( sample_name_list) == x)])
  })))
  
  # Select the sample and subtypes to be analyzed (non WT)
  non_wt_condition_selection =  which( CONDITION_NAME_COLUMN_SET %in% subtype)
  non_wt_sample_selection = colData( dsHTSeq)[ which( colData(dsHTSeq)$population %in% subtype), ]
  se <- SummarizedExperiment( log1pNORMALIZED_COUNTS[ ,non_wt_condition_selection], colData= non_wt_sample_selection)
  pca_and_plot = plotPCA_ade4( DESeqTransform( se ), 
                               intgroup = "population", 
                               individual_label_mapping = REPLICATE_NAME_MAPPING, 
                               legend_label_mapping = POPULATION_NAME_MAPPING,
                               scale_colour = POPULATION_COLOR_MAPPING[ POPULATION_NAME_MAPPING[ subtype]])
  acp = pca_and_plot[[ "pca"]]
  print( pca_and_plot[[ "plot"]])

  cs1_loadings_heatmap = row.names( acp$c1)[ order( abs( acp$c1$CS1), decreasing = TRUE)[1:BEST_NUMBER_HEATMAP]]
  cs2_loadings_heatmap = row.names( acp$c1)[ order( abs( acp$c1$CS2), decreasing = TRUE)[1:BEST_NUMBER_HEATMAP]]
  high_loadings_heatmap = unique( c( cs1_loadings_heatmap, cs2_loadings_heatmap))
  
  # create color palet
  col.pal <- brewer.pal(9,"Blues")

  # define metrics for clustering
  drows1 <- "euclidean"
  dcols1 <- "euclidean"

  # create heatmap of genes on PC1 and PC2
  cat("<H4>Heatmap of samples using genes on PC1 and PC2</H4>")
  heatmap.CS1.CS2_df = as.data.frame( (log1pNORMALIZED_COUNTS[ high_loadings_heatmap, which( CONDITION_NAME_COLUMN_SET %in% subtype)]))
  annotation.CS1.CS2.df = data.frame( Pop = REPLICATE_POPULATION_MAPPING[ names( heatmap.CS1.CS2_df)])
  row.names( annotation.CS1.CS2.df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS1.CS2_df)]
  names( heatmap.CS1.CS2_df) =  REPLICATE_NAME_MAPPING[ names( heatmap.CS1.CS2_df)]
  hm.parameters.CS1.CS2 <- list( heatmap.CS1.CS2_df,
                         color = col.pal,
                         cellwidth =5, cellheight =5,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         fontsize = 6,
                         kmeans_k = NA,
                         show_rownames = T, show_colnames = T,
                         main = "Heatmap (avg, eucl, unsc)",
                         clustering_method = "average",
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         clustering_distance_rows = drows1,
                         clustering_distance_cols = dcols1,
                         annotation_col = annotation.CS1.CS2.df,
                         annotation_colors = list( Pop = POPULATION_COLOR_MAPPING[ POPULATION_NAME_MAPPING[ subtype]]))

  # To draw the heat map on screen
  do.call("pheatmap", hm.parameters.CS1.CS2)
  
  # create heatmap with genes on PC1 only
  cat("<H4>Heatmap of samples using genes on PC1 only</H4>")
  heatmap.CS1_df = as.data.frame( (log1pNORMALIZED_COUNTS[ cs1_loadings_heatmap, which( CONDITION_NAME_COLUMN_SET %in% subtype)]))
  annotation.CS1.df = data.frame( Pop = REPLICATE_POPULATION_MAPPING[ names( heatmap.CS1_df)])
  row.names( annotation.CS1.df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS1_df)]
  names( heatmap.CS1_df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS1_df)]
  hm.parameters.CS1 <- list( heatmap.CS1_df,
                         color = col.pal,
                         cellwidth =5, cellheight =5,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         fontsize = 6,
                         kmeans_k = NA,
                         show_rownames = T, show_colnames = T,
                         main = "Heatmap (avg, eucl, unsc)",
                         clustering_method = "average",
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         clustering_distance_rows = drows1,
                         clustering_distance_cols = dcols1,
                         annotation_col = annotation.CS1.df,
                         annotation_colors = list( Pop = POPULATION_COLOR_MAPPING[ POPULATION_NAME_MAPPING[ subtype]]))
  
  # To draw the heat map on screen
  do.call("pheatmap", hm.parameters.CS1)
  
  # Print the loadings in tables
  print( kable_styling(
     kable( t( acp$c1[ cs1_loadings_heatmap, 1:2]),
            caption="Best loadings in PCA CS1",
            format="html"),
     bootstrap_options = "striped", full_width = F
   )
  )


  print( kable_styling(
     kable( t( acp$c1[ cs2_loadings_heatmap, 1:2]),
            caption="Best loadings in PCA CS2",
            format="html"),
     bootstrap_options = "striped", full_width = F
   )
  )
}

# -------------------------------------------------------------
# For all subpopulation, execute a PCA analysis
# -------------------------------------------------------------

## @knitr qc_on_replicates_heterogeneity_all

sample_name_list = CONDITION_NAME_COLUMN_SET


# Build a color palette for the list of subtypes
color_palette = brewer.pal( length( unique( sample_name_list)),"Paired")
sample_colors = unique( unlist( lapply( sample_name_list, function( x){
  return( color_palette[ which( unique( sample_name_list) == x)])
})))

# Execute the ACP
se <- SummarizedExperiment( log1pNORMALIZED_COUNTS, colData= colData( dsHTSeq))
pca_and_plot = plotPCA_ade4( DESeqTransform( se ), 
                             intgroup = "population",  
                             individual_label_mapping = REPLICATE_NAME_MAPPING,
                             legend_label_mapping = POPULATION_NAME_MAPPING,
                             scale_colour = POPULATION_COLOR_MAPPING)
acp = pca_and_plot[[ "pca"]]
print( pca_and_plot[[ "plot"]])

cs1_loadings = row.names( acp$c1)[ order( abs( acp$c1$CS1), decreasing = TRUE)[1:BEST_NUMBER_HEATMAP]]
cs2_loadings = row.names( acp$c1)[ order( abs( acp$c1$CS2), decreasing = TRUE)[1:BEST_NUMBER_HEATMAP]]
high_loadings = unique( c( cs1_loadings, cs2_loadings))

all_loadings = vector()
for( pcindex in 1: ncol( acp$c1)){
  cs_loadings = row.names( acp$c1)[ order( abs( acp$c1[, pcindex]), decreasing = TRUE)[1:BEST_NUMBER_HEATMAP]]
  all_loadings = append( all_loadings, cs_loadings)
}
all_loadings = unique( all_loadings)


# create color palet
col.pal <- brewer.pal(9,"Blues")

# define metrics for clustering
drows1 <- "euclidean"
dcols1 <- "euclidean"

# create heatmap
cat("<H3>Heatmap of samples with mains loadings on PC1 and PC2</H3>")
heatmap.CS1.CS2_df = as.data.frame( (log1pNORMALIZED_COUNTS[ high_loadings, ]))
annotation.CS1.CS2.df = data.frame( Pop = REPLICATE_POPULATION_MAPPING[ names( heatmap.CS1.CS2_df)])
row.names( annotation.CS1.CS2.df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS1.CS2_df)]
names( heatmap.CS1.CS2_df) =  REPLICATE_NAME_MAPPING[ names( heatmap.CS1.CS2_df)]
hm.parameters.CS1.CS2 <- list( heatmap.CS1.CS2_df,
                               color = col.pal,
                               cellwidth =5, cellheight =5,
                               treeheight_row = 20,
                               treeheight_col = 20,
                               fontsize = 6,
                               kmeans_k = NA,
                               show_rownames = T, show_colnames = T,
                               main = "Heatmap (avg, eucl, unsc)",
                               clustering_method = "average",
                               cluster_rows = TRUE, cluster_cols = TRUE,
                               clustering_distance_rows = drows1,
                               clustering_distance_cols = dcols1,
                               annotation_col = annotation.CS1.CS2.df,
                               annotation_colors = list( Pop = POPULATION_COLOR_MAPPING))

# To draw the heat map on screen
do.call("pheatmap", hm.parameters.CS1.CS2)

cat("<H3>Heatmap of samples with mains loadings on PC1 only</H3>")
heatmap.CS1_df = as.data.frame( (log1pNORMALIZED_COUNTS[ cs1_loadings, ]))
annotation.CS1.df = data.frame( Pop = REPLICATE_POPULATION_MAPPING[ names( heatmap.CS1_df)])
row.names( annotation.CS1.df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS1_df)]
names( heatmap.CS1_df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS1_df)]
hm.parameters.CS1 <- list( heatmap.CS1_df,
                           color = col.pal,
                           cellwidth =5, cellheight =5,
                           treeheight_row = 20,
                           treeheight_col = 20,
                           fontsize = 6,
                           kmeans_k = NA,
                           show_rownames = T, show_colnames = T,
                           main = "Heatmap (avg, eucl, unsc)",
                           clustering_method = "average",
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           clustering_distance_rows = drows1,
                           clustering_distance_cols = dcols1,
                           annotation_col = annotation.CS1.df,
                           annotation_colors = list( Pop = POPULATION_COLOR_MAPPING))

# To draw the heat map on screen
do.call("pheatmap", hm.parameters.CS1)

cat("<H3>Heatmap of samples with mains loadings on PC2 only</H3>")
heatmap.CS2_df = as.data.frame( (log1pNORMALIZED_COUNTS[ cs2_loadings, ]))
annotation.CS2.df = data.frame( Pop = REPLICATE_POPULATION_MAPPING[ names( heatmap.CS2_df)])
row.names( annotation.CS2.df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS2_df)]
names( heatmap.CS2_df) = REPLICATE_NAME_MAPPING[ names( heatmap.CS2_df)]

hm.parameters.CS2 <- list( heatmap.CS2_df,
                           color = col.pal,
                           cellwidth =5, cellheight =5,
                           treeheight_row = 20,
                           treeheight_col = 20,
                           fontsize = 6,
                           kmeans_k = NA,
                           show_rownames = T, show_colnames = T,
                           main = "Heatmap (avg, eucl, unsc)",
                           clustering_method = "average",
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           clustering_distance_rows = drows1,
                           clustering_distance_cols = dcols1,
                           annotation_col = annotation.CS2.df,
                           annotation_colors = list( Pop = POPULATION_COLOR_MAPPING))

# To draw the heat map on screen
do.call("pheatmap", hm.parameters.CS2)

cat("<H3>Heatmap of Spearman correlation between samples using best loadings genes from each PC</H3>")

sample.cor_df <-as.data.frame( cor( log1pNORMALIZED_COUNTS[ all_loadings, ], method=c("spearman")))
annotation.cor.df = data.frame( Pop = REPLICATE_POPULATION_MAPPING[ names( sample.cor_df)])
row.names( annotation.cor.df) = REPLICATE_NAME_MAPPING[ names( sample.cor_df)]
names( sample.cor_df) = REPLICATE_NAME_MAPPING[ names( sample.cor_df)]
row.names( sample.cor_df) = REPLICATE_NAME_MAPPING[ row.names( sample.cor_df)]
cor.parameters <- list( sample.cor_df,
                           #color = c( rev( brewer.pal(3,"Reds")), brewer.pal(9,"Blues")),
                           color = c( "#FF9999", "#FFBBBB", "#FFDDDD", brewer.pal(9,"Blues")),
                           cellwidth =10, cellheight =10,
                           treeheight_row = 20,
                           treeheight_col = 20,
                           fontsize = 10,
                           kmeans_k = NA,
                           show_rownames = T, show_colnames = T,
                           main = "Heatmap (avg, eucl, unsc)",
                           clustering_method = "average",
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           clustering_distance_rows = drows1,
                           clustering_distance_cols = dcols1,
                           annotation_col = annotation.cor.df,
                           annotation_colors = list( Pop = POPULATION_COLOR_MAPPING))
do.call("pheatmap", cor.parameters)

# Plot in table the best loadings ot first two PC
datatable(acp$c1[ cs1_loadings, 1:2], options = list(pageLength = 5),
          caption="Best loadings in PCA CS1")
datatable(acp$c1[ cs2_loadings, 1:2], options = list(pageLength = 5),
          caption="Best loadings in PCA CS2")
