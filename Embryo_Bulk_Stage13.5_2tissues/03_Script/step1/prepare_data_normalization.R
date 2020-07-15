# ###################################################################
# This script aims to normalize the dtata accross all conditions
# using the DEseq noramlization technic
# ###################################################################

## @knitr normalize_data

cat("<H5>Normalize data</H5>")

# Build the dataframe with the columns corresponding to the files associated to the experiments
comparison_df = count_table_all_df[ , EXPERIMENT_FILE_ID_SET]
row.names( comparison_df) = count_table_all_df[ , "external_gene_name"]

# remove lines with no data (only zeros)
comparison_df$sum = apply( comparison_df, 1, sum)
comparison_df = comparison_df[ which( comparison_df$sum > 0),]
comparison_df = comparison_df[ , EXPERIMENT_FILE_ID_SET]

# Change the column name to indicate their association to an experiment name
# and build the dataframe with the list of experiment association
names( comparison_df) = REPLICATE_ID_COLUMN_SET
condition_data = data.frame( population = CONDITION_NAME_COLUMN_SET)

# Create the DESeq2 object
dsHTSeq <- DESeqDataSetFromMatrix(countData=comparison_df ,colData=condition_data, design=~population)
dsHTSeq <- dsHTSeq[ rowSums(counts(dsHTSeq)) > 10, ]
dsHTSeq <- DESeq(dsHTSeq)

dsHTSeq <- estimateSizeFactors( dsHTSeq)

# Get the normalized counts for the whole experiments
NORMALIZED_COUNTS = counts( dsHTSeq, normalized = TRUE)
log1pNORMALIZED_COUNTS = log2(counts(dsHTSeq, normalized=TRUE) + 1)

write.table( NORMALIZED_COUNTS, file= file.path( OUTPUT_DIR, "normalized_counts.csv"), 
             sep=";", quote = FALSE, row.names = TRUE, col.names = TRUE)
