# ################################################
# This script aims to read and filter data
# ################################################

# READ THE DATA
# -----------------------

## @knitr load_data

# Load the data from the Velocyto analysis
input_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( SAMPLE_ID, ".loom"))
ldat <- read.loom.matrices( input_file)
cat("<BR>Analyzed data are located on:", input_file)

# Get the raw data for spliced and unspliced RNA
emat <- ldat$spliced
nmat <- ldat$unspliced

# Look at the prefix and suffix Velocyto added to the cellular barcode
# Habitually, the sampleID is the prefix and a "x" is the suffix e.g. "11230811:AAACCCACATCGGATTx"
loom_prefix = str_match( colnames( emat)[1], "(.*):[A|C|G|T]*(.*)" )[ 1, 2]
loom_suffix = str_match( colnames( emat)[1], "(.*):[A|C|G|T]*(.*)" )[ 1, 3]

if( loom_prefix != SAMPLE_ID){
  cat("<BR><b>WARNING: prefix of cellular barcodes is loom file is not the sample ID: Prefix=", loom_prefix, " while SAMPLE_ID=", SAMPLE_ID, "</b>")
}

### Load the data from the Seurat analysis

# Load normalized expression table and change the cell barcode to adapt to velocyto output
SEURAT_NORMALIZED_EXPRESSION_MATRIX = read.table( PATH_SEURAT_NORMALIZED_EXPRESSION_MATRIX, header=TRUE, row.names = 1, sep=",")
names( SEURAT_NORMALIZED_EXPRESSION_MATRIX) = paste0( loom_prefix, ":", names( SEURAT_NORMALIZED_EXPRESSION_MATRIX), loom_suffix)

# Load the PCA cell coordinates and change the cell barcode to adapt to velocyto output
SEURAT_PCA_EMBEDDING = read.table( PATH_SEURAT_PCA_EMBEDDING, header=TRUE, row.names = 1, sep=",")
row.names( SEURAT_PCA_EMBEDDING) = paste0( loom_prefix, ":", row.names( SEURAT_PCA_EMBEDDING), loom_suffix)

# Load the t-SNE embedding (with cell type infered by Monocle) and change the cell barcode to adapt to velocyto output
SEURAT_TSNE_EMBEDDING_WITHCELLTYPE = read.table( PATH_SEURAT_TSNE_EMBEDDING_WITHCELLTYPE, header=TRUE, row.names = 1, sep=",")
row.names( SEURAT_TSNE_EMBEDDING_WITHCELLTYPE) = paste0( loom_prefix, ":", row.names( SEURAT_TSNE_EMBEDDING_WITHCELLTYPE), loom_suffix)

# Load the UMAP embedding (with cell type infered by Monocle) and change the cell barcode to adapt to velocyto output
SEURAT_UMAP_EMBEDDING_WITHCELLTYPE = read.table( PATH_SEURAT_UMAP_EMBEDDING_WITHCELLTYPE, header=TRUE, row.names = 1, sep=",")
row.names( SEURAT_UMAP_EMBEDDING_WITHCELLTYPE) = paste0( loom_prefix, ":", row.names( SEURAT_UMAP_EMBEDDING_WITHCELLTYPE), loom_suffix)

# Load the Seurat cell clusters mapping and change the cell barcode to adapt to velocyto output
SEURAT_CLUSTERS = read.table( PATH_SEURAT_CLUSTERS, header=TRUE, row.names = 1, sep=",")
row.names( SEURAT_CLUSTERS) = paste0( loom_prefix, ":", row.names( SEURAT_CLUSTERS), loom_suffix)


# FILTER THE DATA
# -----------------------

## @knitr filter_data

# Look at the cells kept in the data from previous analysis (after QC and contamination filtering)
cells_name_to_keep = row.names( SEURAT_TSNE_EMBEDDING_WITHCELLTYPE)

# Restrict to cells that passed pagoda2 filter
emat <- emat[, cells_name_to_keep]
nmat <- nmat[, cells_name_to_keep]

cat("<BR>Number of valid cells for analysis in emat =",length(colnames(emat)), "<BR>")
cat("<BR>Number of valid cells for analysis in nmat =",length(colnames(nmat)), "<BR>")
cat("<BR>Number of valid cells for analysis =",length(intersect(colnames(emat),colnames(nmat))), "<BR>")



