# ###########################################################################
# This script aims to build the pseudo-time spanning tree from Monocle 
# package
# ###########################################################################

## @knitr build_peudotime_tree_without_marker

diff_test_res <- differentialGeneTest( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED[ EXPRESSED_GENES,],
                                       fullModelFormulaStr = "~Cluster")
ordering_genes <- row.names ( subset(diff_test_res, qval < 0.01))

FILTERED_CDS_PSEUDOTIME_UNSUPERVISED <- setOrderingFilter( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, ordering_genes)
plot_ordering_genes( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED)


FILTERED_CDS_PSEUDOTIME_UNSUPERVISED <- reduceDimension( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, max_components = 2,
                            method = 'DDRTree')
FILTERED_CDS_PSEUDOTIME_UNSUPERVISED <- orderCells( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED)

plot_cell_trajectory( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, color_by = "State")
plot_cell_trajectory( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, color_by = "Cluster")
plot_cell_trajectory( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, color_by = "Cluster") + facet_wrap( ~Cluster)

#
# Plot pseudotime tree by CellType
#
cat("<H5>Plotting pseudotime tree by Cell Type</H5>")
pData( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED)$SupervisedCellType = rep( "NA", length( pData( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED)))
common_cells = intersect( row.names( pData( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED)), row.names( pData( FILTERED_CDS_PSEUDOTIME_SUPERVISED)))
pData( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED)[ common_cells, "SupervisedCellType"] = as.character( pData( FILTERED_CDS_PSEUDOTIME_SUPERVISED)[ common_cells, "CellType"])

plot_cell_trajectory( FILTERED_CDS_PSEUDOTIME_UNSUPERVISED, color = "SupervisedCellType") +
  facet_wrap(~SupervisedCellType) +
  theme(legend.position="none")

