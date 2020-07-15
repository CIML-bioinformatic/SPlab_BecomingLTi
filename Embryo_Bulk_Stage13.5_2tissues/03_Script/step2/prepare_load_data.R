# #############################################################################
# This script aims to load the data from normalized data file produced in
# the previous step.
# #############################################################################


## @knitr load_data

# Prepare datframe with replicate ID to condition association
EXPERIMENT_TABLE_DF = data.frame( 
		condition.id =c(  "EP1", "EP1", "EP1", "EP1", "EP2", "EP2", "EP2", "EP2",
				          "EP3", "EP3", "EP3", "EP3", "EP4", "EP4", "EP4", "EP4",
						  "FLP1", "FLP1", "FLP1", "FLP2", "FLP2", "FLP2",
						  "FLP3", "FLP3", "FLP3", "FLP4", "FLP4"),
		replicate.id = c( "EP1.1", "EP1.2", "EP1.3", "EP1.4", "EP2.1", "EP2.2", "EP2.3", "EP2.4",
				          "EP3.1", "EP3.2", "EP3.3", "EP3.4", "EP4.1", "EP4.2", "EP4.3", "EP4.4",
				          "FLP1.1", "FLP1.2", "FLP1.3", "FLP2.1", "FLP2.2", "FLP2.3",
						  "FLP3.1", "FLP3.2", "FLP3.3", "FLP4.1", "FLP4.3"),
		stringsAsFactors = FALSE
		)

CONDITION_ID_SET = unique( EXPERIMENT_TABLE_DF$condition.id)
CONDITION_NAME_COLUMN_SET = EXPERIMENT_TABLE_DF$condition.id
		
# Load the data from normalized data count produced in previous step
input_file = file.path( INPUT_DIR, "normalized_counts.csv")
cat("<BR>Reading normalized data from file:<BR>", input_file)
raw_normalized_data_df = read.table( input_file, sep=";", header=TRUE, row.names=1)
NORMALIZED_DATA_DF = raw_normalized_data_df[ , EXPERIMENT_TABLE_DF$replicate.id]

# Build a dataframe with the mean of each population over their replicates and the global
# mean and standard deviation over populations
NORMALIZED_MEAN_DF = do.call( rbind, apply( NORMALIZED_DATA_DF, 1, function( row){
					
					values = lapply( CONDITION_ID_SET, function( condition_name){
								return( mean( row[ which( CONDITION_NAME_COLUMN_SET == condition_name & row > 0)]))
							})
					
					values[ which( is.nan( unlist( values)))] = 0
					names( values) = CONDITION_ID_SET
					values = data.frame( values, stringsAsFactors = FALSE)
					return( values)
				}))
NORMALIZED_MEAN_DF$global.mean = apply( NORMALIZED_MEAN_DF, 1, mean)
NORMALIZED_MEAN_DF$global.sd = apply( NORMALIZED_MEAN_DF, 1, sd)

