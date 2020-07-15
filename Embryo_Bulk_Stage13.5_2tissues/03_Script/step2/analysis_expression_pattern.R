# ################################################################
# This script aims to classify the gene by their expression
# pattern along population evolution
# ################################################################

## @knitr expression_pattern

# Determine the evolution pattern of each genes along the WT and KO-JUN conditions
# for instance "11" = increase between point1 and point2 anf between point2 and point3
# for instance "1-1" = increase between point1 and point2 and decrease between point2 and point3
pattern_df = do.call( rbind, lapply( row.names( NORMALIZED_MEAN_DF), function( current_gene_name){
  
  pattern_EP = ""
  pattern_EP = paste0( pattern_EP, sign( NORMALIZED_MEAN_DF[ current_gene_name, "EP2"] - NORMALIZED_MEAN_DF[ current_gene_name, "EP1"]))
  pattern_EP = paste0( pattern_EP, sign( NORMALIZED_MEAN_DF[ current_gene_name, "EP3"] - NORMALIZED_MEAN_DF[ current_gene_name, "EP2"]))
  pattern_EP = paste0( pattern_EP, sign( NORMALIZED_MEAN_DF[ current_gene_name, "EP4"] - NORMALIZED_MEAN_DF[ current_gene_name, "EP3"]))
  mean_EP = mean( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "EP1", "EP2", "EP3", "EP4")]))
  sd_EP = sd( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "EP1", "EP2", "EP3", "EP4")]))
  cv_EP = sd_EP/mean_EP
  max_diff_EP = max( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "EP1", "EP2", "EP3", "EP4")])) - min( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "EP1", "EP2", "EP3", "EP4")]))
  
  
  Patter_FLP = ""
  Patter_FLP = paste0( Patter_FLP, sign( NORMALIZED_MEAN_DF[ current_gene_name, "FLP2"] - NORMALIZED_MEAN_DF[ current_gene_name, "FLP1"]))
  Patter_FLP = paste0( Patter_FLP, sign( NORMALIZED_MEAN_DF[ current_gene_name, "FLP3"] - NORMALIZED_MEAN_DF[ current_gene_name, "FLP2"]))
  Patter_FLP = paste0( Patter_FLP, sign( NORMALIZED_MEAN_DF[ current_gene_name, "FLP4"] - NORMALIZED_MEAN_DF[ current_gene_name, "FLP3"]))
  mean_FLP = mean( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "FLP1", "FLP2", "FLP3", "FLP4")]))
  sd_FLP = sd( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "FLP1", "FLP2", "FLP3", "FLP4")]))
  cv_FLP = sd_FLP/mean_FLP
  max_diff_FLP = max( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "FLP1", "FLP2", "FLP3", "FLP4")])) - min( as.numeric( NORMALIZED_MEAN_DF[ current_gene_name, c( "FLP1", "FLP2", "FLP3", "FLP4")]))
  
  max_diff = max( max_diff_EP, max_diff_FLP)
  return( data.frame( patternEP = pattern_EP, maxDiffEP = max_diff_EP, meanEP = mean_EP, sdEP = sd_EP, cvEP = cv_EP,
                      patternFLP = Patter_FLP, maxDiffFLP = max_diff_FLP, meanFLP = mean_FLP, sdFLP = sd_FLP, cvFLP = cv_FLP,
                      maxDiff = max_diff, stringsAsFactors = FALSE))
}))
row.names( pattern_df) = row.names( NORMALIZED_MEAN_DF)


# Identify the different patterns and plot the genes with the same patterns  on EP populations
pattern_set = sort( unique( pattern_df$patternEP))

# Look at the evolution of monitored genes
cat("<H3>Evolution of some selected genes grouped by pattern type along EP2 to EP4</H3>")
for( pattern in pattern_set){
  
  # -- Select the monitored genes with selected pattern
  monitored_genes_with_pattern_df = pattern_df[ which( pattern_df$patternEP == pattern & row.names( pattern_df) %in% MONITORED_GENES), ]
  # -- If there is some monitored genes with the selected pattern, show them
  if( nrow( monitored_genes_with_pattern_df) > 0){
    # -- Get the list of monitored genes with the right pattern
    patterned_gene_set = row.names( monitored_genes_with_pattern_df)
    # -- Select only values of EP population in dataframe
    selected_df = NORMALIZED_MEAN_DF[ patterned_gene_set, c( "EP2", "EP3", "EP4") ]
    selected_df$gene_name = row.names( selected_df)
    # -- Build the long version of the dataframe to use it with ggplot
    selected_long_df = melt( selected_df, id=c( "gene_name"))
    selected_long_df = selected_long_df[ order( selected_long_df$gene_name),]
    # -- Plot the graph of selected genes evolutions
    print( ggplot( selected_long_df ) +
             geom_point( aes( x=variable, y=value, group = gene_name, color = gene_name)) +
             geom_line( aes( x=variable, y=value, group = gene_name, color = gene_name)) +
             ggtitle(  paste( "Pattern:", pattern)) +
             theme( legend.position="bottom") +
             xlab( "Population") + ylab( "Normalized expression") + 
             scale_x_discrete( labels = POPULATION_NAME_MAPPING) + 
             theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black")) +
             theme( axis.text.x = element_text( size = AXIS_TICK_SIZE),
                    axis.text.y = element_text( size = AXIS_TICK_SIZE)) +
             theme( axis.title.x = element_text( size = AXIS_TITLE_SIZE),
                    axis.title.y = element_text( size = AXIS_TITLE_SIZE))+ 
             theme( legend.text = element_text( size=LEGEND_TEXT_SIZE))
    )
  } 
}

# Look at the evolution of genes by functionnal categories
cat("<H3>Evolution of functionnal groups of genes</H3>")
for( category in names( GENE_FUNCTIONAL_CATEGORIES)){
  # -- Get the list of genes to plot
  raw_category_gene_set = GENE_FUNCTIONAL_CATEGORIES[[ category]]
  category_gene_set = intersect( raw_category_gene_set, row.names( NORMALIZED_MEAN_DF))
  
  # -- Select only values of EP or FLP populations in dataframe
  for( population_type_name in names( POPULATION_TYPE)){
    selected_df = NORMALIZED_MEAN_DF[ category_gene_set, POPULATION_TYPE[[ population_type_name]] ]
    selected_df$gene_name = row.names( selected_df)
    # -- Build the long version of the dataframe to use it with ggplot
    selected_long_df = melt( selected_df, id=c( "gene_name"))
    selected_long_df = selected_long_df[ order( selected_long_df$gene_name),]
    # -- Plot the graph of selected genes evolutions
    print( ggplot( selected_long_df ) +
             geom_point( aes( x=variable, y=value, group = gene_name, color = gene_name)) +
             geom_line( aes( x=variable, y=value, group = gene_name, color = gene_name)) +
             ggtitle(  paste( gsub( "_", " ", category), "on", population_type_name)) +
             theme( legend.position="bottom") +
             xlab( "Population") + ylab( "Normalized expression") + 
             scale_x_discrete( labels = POPULATION_NAME_MAPPING) +
             theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black")) +
             theme( axis.text.x = element_text( size = AXIS_TICK_SIZE),
                    axis.text.y = element_text( size = AXIS_TICK_SIZE)) +
             theme( axis.title.x = element_text( size = AXIS_TITLE_SIZE),
                    axis.title.y = element_text( size = AXIS_TITLE_SIZE))+ 
             theme( legend.text = element_text( size=LEGEND_TEXT_SIZE))
    )
  }
}

# Parse the list of patterns to display the corresponding genes
cat("<H3>Evolution of most variable genes grouped by pattern type along EP1 to EP4</H3>")
for( pattern in pattern_set){
  
  # -- Select the genes with the current pattern
  selected_pattern_df = pattern_df[ which( pattern_df$patternEP == pattern),]
  # -- Get the genes with the largest evolutions (top 20)
  selected_pattern_df = selected_pattern_df[ order( selected_pattern_df$maxDiffEP, decreasing = TRUE), ]
  selected_pattern_df = selected_pattern_df[ 1:20,]
  quantiles = quantile( selected_pattern_df$maxDiffEP)
  patterned_gene_set = row.names( selected_pattern_df)
  # -- Select only values of EP population in dataframe
  selected_df = NORMALIZED_MEAN_DF[ patterned_gene_set, c( "EP1", "EP2", "EP3", "EP4")]
  selected_df$gene_name = row.names( selected_df)
  # -- Build the long version of the dataframe to use it with ggplot
  selected_long_df = melt( selected_df, id=c( "gene_name"))
  selected_long_df = selected_long_df[ order( selected_long_df$gene_name),]
  # -- Plot the graph of selected genes evolutions
  print( ggplot( selected_long_df ) +
           geom_point( aes( x=variable, y=value, group = gene_name, color = gene_name)) +
           geom_line( aes( x=variable, y=value, group = gene_name, color = gene_name)) +
           ggtitle(  paste( "Pattern:", pattern), paste( "Quartiles of maximal déviation:", paste( signif( quantiles, 3), collapse =" ; "))) +
           theme( legend.position="bottom") +
           xlab( "Population") + ylab( "Normalized expression") + 
           scale_x_discrete( labels = POPULATION_NAME_MAPPING) +
           theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
           theme( axis.text.x = element_text( size = AXIS_TICK_SIZE),
                  axis.text.y = element_text( size = AXIS_TICK_SIZE)) +
           theme( axis.title.x = element_text( size = AXIS_TITLE_SIZE),
                  axis.title.y = element_text( size = AXIS_TITLE_SIZE))+ 
           theme( legend.text = element_text( size=LEGEND_TEXT_SIZE))
  )

}

