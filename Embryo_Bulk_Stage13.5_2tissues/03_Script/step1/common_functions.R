

## @knitr load_common_functions

#
# Function that convert genes ids from one world 'ENSEMBL, ENTREZID, genesymobol) to an other
#
# ids = list of IDS
# fromKey = key type; toKey = key type we want to convert to
# db = the AnnotationDb object to use.
# ifMultiple = the argument specifies what to do if one source ID maps to several target IDs:
#              should the function return an NA or simply the first of the multiple IDs?
#
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] 
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


# 
# Redefine the plotPCA function of DeSeq2 to use ade4
# 
plotPCA_ade4 = function (object, 
                         intgroup = "condition", 
                         ntop = 500, 
                         individual_label_mapping = NULL, 
                         legend_label_mapping = NULL,
                         scale_colour = NULL) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca_ade4 = dudi.pca( t(assay(object)[select, ]), scannf = FALSE, nf = nrow( object), scale=FALSE)
  percentVar_ade4 <- pca_ade4$eig/sum(pca_ade4$eig)
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  if (length(intgroup) > 1) {
    group = factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    group = colData(object)[[intgroup]]
  }
  d_ade4 <- data.frame(PC1 = pca_ade4$l1$RS1, PC2 = pca_ade4$l1$RS2, group = group, 
                       intgroup.df, name = colnames(object))
  
  plot = ggplot(data = d_ade4, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + 
    geom_text_repel( aes( label= individual_label_mapping[ as.character( name)]),hjust=0, vjust=0) +
    xlab(paste0("PC1: ", round(percentVar_ade4[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar_ade4[2] * 100), "% variance")) +
    coord_fixed() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    theme( axis.text.x = element_text( size = AXIS_TICK_SIZE),
           axis.text.y = element_text( size = AXIS_TICK_SIZE)) +
    theme( axis.title.x = element_text( size = AXIS_TITLE_SIZE),
           axis.title.y = element_text( size = AXIS_TITLE_SIZE))+ 
    theme( legend.text = element_text( size=LEGEND_TEXT_SIZE))
  
  # If a mapping is provided, use it to change the legend labels
   if( !is.null( legend_label_mapping)){
     mapped_labels= unlist( legend_label_mapping[ intersect( levels( d_ade4$group), unique( d_ade4$group))], use.names = FALSE)
     if( ! is.null( scale_colour)){
       plot = plot + scale_colour_manual( values = unname( scale_colour), labels= mapped_labels)
     }else{
       plot = plot + scale_color_manual( values = brewer.pal( name = "Dark2", n = length( mapped_labels)), labels= mapped_labels)
     }
   }
  

  
  return( list( pca = pca_ade4, plot = plot))
}

# 
# Convert the names of the name_set from names in origin_set to names in target_set
# If keep is TRUE, the original name is put in the result next to the converted name bwteen parethesis
# Names that are not in the origin_set are kept as is
# 
convert_names <- function( name_set, origin_set, target_set, keep=TRUE){
  
  result_name_set = vector()
  for( current_name in name_set){
    current_index = which( origin_set == current_name)
    if( length( current_index) >0 && current_index > 0){
      if( keep){
        result_name_set = append( result_name_set, paste0( target_set[ current_index], " (", current_name, ")"))
      }else{
        result_name_set = append( result_name_set, target_set[ current_index])
      }
    }else{
      result_name_set = append( result_name_set, current_name)
    }
  }
  return( result_name_set)
}