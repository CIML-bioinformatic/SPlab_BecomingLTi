# ###############################################################
# This script aims to use Monocle method to order cells
# with pseudotime
# ###############################################################

# This script is inpired by the Monocle tutorial:
# http://cole-trapnell-lab.github.io/monocle-release/docs/

## 
# Define function to get a good color palette for cell counts
##
build_color_palette <- function( colors, local_emb, gradient.range.quantile=0.95){
  
  if(all(sign(colors)>=0)) {
    gradientPalette <- colorRampPalette(c('gray80','red'), space = "Lab")(1024)
  } else {
    gradientPalette <- colorRampPalette(c("blue", "grey70", "red"), space = "Lab")(1024)
  }
  if(all(sign(colors)>=0)) {
    zlim <- as.numeric(quantile(colors,p=c(1-gradient.range.quantile,gradient.range.quantile)))
    if(diff(zlim)==0) {
      zlim <- as.numeric(range(colors))
    }
  } else {
    zlim <- c(-1,1)*as.numeric(quantile(abs(colors),p=gradient.range.quantile))
    if(diff(zlim)==0) {
      zlim <- c(-1,1)*as.numeric(max(abs(colors)))
    }
  }
  # restrict the values
  colors[colors<zlim[1]] <- zlim[1]; colors[colors>zlim[2]] <- zlim[2];
  
  colors <- (colors-zlim[1])/(zlim[2]-zlim[1])
  cols <- gradientPalette[colors[match(rownames(local_emb),names(colors))]*(length(gradientPalette)-1)+1]
  
  return( cols)
}


## @knitr pseudotime_analysis
set.seed( 123)

cat("<H4>Computing pseudo-time representation</H4>")
# Order cells by pseudotime using markers genes from previous analysis
FILTERED_CDS_PSEUDOTIME_SUPERVISED <- setOrderingFilter( filtered_cds_for_pseudotime, signature_genes_for_pseudotime)
plot_ordering_genes( FILTERED_CDS_PSEUDOTIME_SUPERVISED) 
FILTERED_CDS_PSEUDOTIME_SUPERVISED <- reduceDimension( FILTERED_CDS_PSEUDOTIME_SUPERVISED, max_components = 2,
                            method = 'DDRTree')
FILTERED_CDS_PSEUDOTIME_SUPERVISED <- orderCells( FILTERED_CDS_PSEUDOTIME_SUPERVISED)


#
# Plot pseudotime tree by Pseudotime
#
cat("<H5>Plotting pseudotime tree by Pseudotime</H5>")
plot_cell_trajectory( FILTERED_CDS_PSEUDOTIME_SUPERVISED, color = "Pseudotime")

#
# Plot pseudotime tree by CellType
#
cat("<H5>Plotting pseudotime tree by Cell Type</H5>")
plot_cell_trajectory( FILTERED_CDS_PSEUDOTIME_SUPERVISED, color = "CellType") + facet_wrap(~CellType) + theme(legend.position="none")

#
# Plot pseudotime tree by State
#
cat("<H5>Plotting pseudotime tree by State</H5>")
state_colors = c( "firebrick2", "firebrick", "dodgerblue1", "dodgerblue4", "chartreuse3")
plot_cell_trajectory( FILTERED_CDS_PSEUDOTIME_SUPERVISED, color = "State") +
  scale_color_manual(breaks = c("1", "2", "3", "4", "5"), values=state_colors) +
  theme(legend.position = "none")

#
# Show table of distribution of cell type among pseudotime states
#
cat("<H5>Table of distribution of cell type among pseudotime states</H5>")
supervised_cell_state_df = dcast( data.frame( table( pData( FILTERED_CDS_PSEUDOTIME_SUPERVISED)[ , c( "CellType", "State")]))
                       , CellType ~ State, value.var="Freq")
row.names( supervised_cell_state_df) = supervised_cell_state_df$CellType
supervised_cell_state_df = supervised_cell_state_df[ , -which( names( supervised_cell_state_df) == "CellType")]
supervised_cell_state_df = supervised_cell_state_df[ -which( row.names( supervised_cell_state_df) == "Unknown"),]
datatable(  supervised_cell_state_df , caption = "Distribution of cell type among pseudotime states")
d3heatmap( supervised_cell_state_df, colors = "Blues", Rowv = FALSE, Colv = FALSE)

