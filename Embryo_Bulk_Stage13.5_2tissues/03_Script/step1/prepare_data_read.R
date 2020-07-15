
## @knitr read_data

DATA_FILE_1 = file.path( RAW_DATA_DIR, "0025_AHVYCNBGX2_FeatureCounts", "170831_NB501865_0025_AHVYCNBGX2_clean.count")
DATA_FILE_2 = file.path( RAW_DATA_DIR, "0026_AHVW33BGX2_FeatureCounts", "170905_NB501865_0026_AHVW33BGX2.count")

cat("<H5>Load data</H5>")
cat("<BR>The data will be loaded from:")
cat("<BR> *", DATA_FILE_1)
cat("<BR> *", DATA_FILE_2)

INTERESTING_PATHWAYS = c( "mmu04350", "mmu04330", "mmu00100", "mmu04062", "mmu04310", "mmu04110")
INTERESTING_GENES = c("Il17re", "Il17f", "Trance", "TranceR", "CD4", "RORc", "Cxcr6", "Cxcr5", "Cxcr4", "Ccr7", "Gata3", "Flt3", "Csf2", "Tcf7", "Id2", "Nr1d1", "Ets1", "LTa")

# Create the mapping between experiment name, experiment id and expriment files
experiment_table_df = data.frame( 
  experiment.id =c( "EP1", "EP2", "EP3", "EP4",
                    "EP1", "EP2", "FLP1", "FLP2",
                    "EP1", "EP2", "EP3", "EP4",
                    "EP1", "EP2", "EP3", "EP4", 
                    "FLP1", "FLP2", "FLP3", "FLP4", 
                    "EP3", "EP4", "FLP3", "FLP4", 
                    "FLP1", "FLP2", "FLP3", "FLP4"),
  sample.id = c(25,26,27,28,
                29,30,41,46,
                33,34,35,36,
                37,38,39,40,
                13,14,15,48,
                31,32,43,44,
                49,50,51,52),
  sample.file.id=as.character(c(10333525,10333526,10333527,10333528,
                                10333529,10333530,10333545,10333546,
                                10333533,10333534,10333535,10333536,
                                10333537,10333538,10333539,10333540,
                                10333541,10333542,10333543,10333544,
                                10333531,10333532,10333547,10333548,
                                10333549,10333550,10333551,10333552)), stringsAsFactors = FALSE)

EXPERIMENT_ID_SET = as.character( unique( experiment_table_df$experiment.id))


# Read the first data file of counts and change the column names to match with experiment files
count_table_1_df = read.table( file = DATA_FILE_1, header = TRUE, sep="\t")

original_names = names( count_table_1_df)
new_names =vector()
for( name in original_names){
  match = str_match( name, ".*\\.([0-9]+)_S[0-9]+\\.rmdup.bam")
  if( !is.na(match[ ,2])){
    new_names = append( new_names, match[ ,2])
  }else{
    new_names = append( new_names, name)
  }
}
names( count_table_1_df) = new_names

# Read the second data file of counts and change the column names to match with experiment files
count_table_2_df = read.table( file = DATA_FILE_2, header = TRUE, sep="\t")

original_names = names( count_table_2_df)
new_names =vector()
for( name in original_names){
  match = str_match( name, ".*\\.([0-9]+)_S[0-9]+\\.rmdup.bam")
  if( !is.na(match[ ,2])){
    new_names = append( new_names, match[ ,2])
  }else{
    new_names = append( new_names, name)
  }
}
names( count_table_2_df) = new_names

# Merge the data from both data files in a single data frame
count_table_all_df = merge( count_table_1_df, count_table_2_df, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"), all.x=TRUE, all.y=TRUE)

# Associate a gene name to the uniprot ID
# -- interrogate the BioMart database to get gene symbol and description for these genes 
#my_mart = useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org")
my_mart = useMart("ENSEMBL_MART_ENSEMBL",host="aug2017.archive.ensembl.org")
mart<- useDataset( "mmusculus_gene_ensembl", my_mart)
genes.IDs <- getBM(
  attributes= c("ensembl_gene_id","external_gene_name"),
  filters= "ensembl_gene_id",
  values= count_table_all_df$Geneid,
  mart= mart
)
idmap = merge(x = count_table_all_df , y = genes.IDs,by.x="Geneid", by.y="ensembl_gene_id", all.x=TRUE)

# -- identify genes with no gene_name or with duplicated gene names and replaces their name by their original ID
NoID=which(is.na(idmap$external_gene_name))
idmap$external_gene_name[NoID]=as.character(idmap$Geneid[NoID])
DuplicatedGenes=which(duplicated(idmap$external_gene_name))
idmap$external_gene_name[DuplicatedGenes]=as.character(idmap$Geneid[DuplicatedGenes])
count_table_all_df = idmap

# Get the list of replicate files in order of the experiments
raw_experiment_file_id_set = vector()
RAW_REPLICATE_ID_COLUMN_SET = vector()
RAW_CONDITION_NAME_COLUMN_SET = vector()
for( exp_id in EXPERIMENT_ID_SET){
  replicate_set = experiment_table_df[ which( experiment_table_df$experiment.id == exp_id), "sample.file.id"]
  raw_experiment_file_id_set = append( raw_experiment_file_id_set, replicate_set)
  RAW_REPLICATE_ID_COLUMN_SET = append( RAW_REPLICATE_ID_COLUMN_SET, paste( rep( exp_id, length( replicate_set)), seq( 1, length( replicate_set),1), sep="."))
  RAW_CONDITION_NAME_COLUMN_SET = append( RAW_CONDITION_NAME_COLUMN_SET, rep( exp_id, length( replicate_set)))
}

# Filter out the undesired replicates
cat("<H4>Elimination of unwanted replicates</H4>")
cat("<BR>The following replicates will be removed from analysis:", paste( UNWANTED_REPLICATES, collapse=";"))
unwanted_replicates_index = which( RAW_REPLICATE_ID_COLUMN_SET %in% UNWANTED_REPLICATES)
EXPERIMENT_FILE_ID_SET = raw_experiment_file_id_set[ -unwanted_replicates_index]
REPLICATE_ID_COLUMN_SET = RAW_REPLICATE_ID_COLUMN_SET[ -unwanted_replicates_index]
CONDITION_NAME_COLUMN_SET = RAW_CONDITION_NAME_COLUMN_SET[ -unwanted_replicates_index]




