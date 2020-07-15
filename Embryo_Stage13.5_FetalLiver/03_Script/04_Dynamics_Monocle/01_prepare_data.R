# ################################################
# This script aims to read and filter data
# ################################################

# This script is inpired by the Monocle tutorial:
# http://cole-trapnell-lab.github.io/monocle-release/docs/

# READ THE DATA
# -----------------------

## @knitr load_data
raw_expression_df = read.csv( file = PATH_RAW_EXPRESSION_MATRIX_FILE, header = TRUE, sep=",")
cat("<BR>Analyzed data are located on:", PATH_RAW_EXPRESSION_MATRIX_FILE)
cat("<BR>Number of cells:", ncol( raw_expression_df))
cat("<BR>Number of genes:", nrow( raw_expression_df))

TSNE_COORD_DF = read.csv( file = PATH_TSNE_EMBEDDING_FILE, header = TRUE, sep=",")
UMAP_COORD_DF = read.csv( file = PATH_UMAP_EMBEDDING_FILE, header = TRUE, sep=",")

# BUILD THE DATA
# ................

## @knitr build_data
cat("<H5>Building the required data for Monocle</H5>")

feature_df = AnnotatedDataFrame( data.frame( gene_short_name = row.names( raw_expression_df)))
row.names( feature_df) = row.names( raw_expression_df)

# Build the CDS (Cell data set) required by monocle
filtered_cds <- newCellDataSet( as.matrix( raw_expression_df), 
                                featureData = feature_df,
                                expressionFamily=negbinomial.size())

# Compute size factors and dispersion
filtered_cds = estimateSizeFactors( filtered_cds)
filtered_cds = estimateDispersions( filtered_cds)

# FILTER CELLS AND GENES
# ......................

## @knitr filter_data
cat("<H5>Filtering the low quality data</H5>")

# Select the genes with at least 10 cells expressing it at a minimum of 0.1
filtered_cds <- detectGenes( filtered_cds, min_expr = 0.1)
EXPRESSED_GENES <- row.names( subset( fData( filtered_cds), num_cells_expressed >= 10))
print( datatable( fData( filtered_cds)[ EXPRESSED_GENES, ]))

# Look at the total number of mRNA that display the limits of mean +/- 2sd
pData( filtered_cds)$Total_UMI <- Matrix::colSums( Biobase::exprs( filtered_cds))
upper_bound <- 10^(mean(log10(pData( filtered_cds)$Total_UMI)) +
                     2*sd(log10(pData( filtered_cds)$Total_UMI)))
lower_bound <- 10^(mean(log10(pData( filtered_cds)$Total_UMI)) -
                     2*sd(log10(pData( filtered_cds)$Total_UMI)))

print( qplot( Total_UMI, data = pData( filtered_cds), geom = "density") +
         geom_vline(xintercept = lower_bound) +
         geom_vline(xintercept = upper_bound)
)

# Look if the distribution of the expression of expressed genes is lognormal
# -- Log-transform each value in the expression matrix.
log_exprs <- log( Biobase::exprs( filtered_cds)[ EXPRESSED_GENES,])

# -- Standardize each gene, so that they are all on the same scale,
# -- Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt( Matrix::t( scale( Matrix::t( log_exprs))))

# -- Plot the distribution of the standardized gene expression values.
print( qplot(value, geom = "density", data = melted_dens_df) +
         stat_function(fun = dnorm, size = 0.5, color = 'red') +
         xlab("Standardized log( expr)") +
         ylab("Density")
)





