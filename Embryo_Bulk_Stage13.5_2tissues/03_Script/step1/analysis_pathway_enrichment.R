# ############################################################
# This script aims to look at pathway enrichment
# ############################################################

## @knitr pathway_enrichment

## Bimap interface:
x <- org.Mm.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list and build a mapping dataframe
mapped_genes_list <- as.list(x[mapped_genes])
mapping_df =  do.call( rbind, lapply( names( mapped_genes_list), function( entry_name){
  entry_name = as.character( entry_name)
  entry = mapped_genes_list[[ entry_name]] 
  current_df = data.frame( ensid = unlist( entry, use.names = FALSE), entrezid=rep( entry_name, length( entry)))
  return( current_df)
}))

if(length( mapped_genes_list) > 0) {
  count_table_all_df$entrez_id = unlist( sapply( count_table_all_df$Geneid, function( ensembl_id){
    ensembl_id = as.character( ensembl_id)
    id_index = which( mapping_df$ensid == ensembl_id)
    if( length( id_index) == 1 ){
      return( as.character( mapping_df$entrezid[ id_index]))
    }else{
      return( NA)
    }
  }))
}

# Load mouse KEGG pathways
kg.mmu <- kegg.gsets( "mmu" )
kegg.gs2 <- kg.mmu$kg.sets[ kg.mmu$sigmet.idx ]

# Prepare dataframe
# Build the dataframe with the columns corresponding to the files associated to the experiments
mapped_rows = which( !is.na( count_table_all_df$entrez_id) & !duplicated( count_table_all_df$entrez_id))
E = count_table_all_df[ mapped_rows, EXPERIMENT_FILE_ID_SET]
row.names( E) = count_table_all_df[ mapped_rows, "entrez_id"]


# remove lines with no data (only zeros)
E$sum = apply( E, 1, sum)
E = E[ which( E$sum > 0),]
E = E[ , EXPERIMENT_FILE_ID_SET]

# Change the column name to indicate their association to an experiment name
# and build the dataframe with the list of experiment association
names( E) = REPLICATE_ID_COLUMN_SET

# For each pair of condition, get the most enriched KEGGS pathways 
# in genes differentially expressed between conditions
for( i in 1:(length( EXPERIMENT_ID_SET)-1)){
  
  for( j in (i+1): length( EXPERIMENT_ID_SET)){

    # Get the condition names
    ref_name = EXPERIMENT_ID_SET[ i]
    samp_name = EXPERIMENT_ID_SET[ j]
    
    # analyze the enriched matways
    ref_col_indexes = which( CONDITION_NAME_COLUMN_SET == ref_name)
    samp_col_indexes = which( CONDITION_NAME_COLUMN_SET == samp_name)
    res <- gage( E, gsets= kegg.gs2, 
                 ref= ref_col_indexes, 
                 samp= samp_col_indexes,
                 compare= "unpaired", same.dir= FALSE )
    
    # Print a table with the top most enriched patways
    print( kable_styling( 
      kable( data.frame( qval=res$greater[ 1:10, "q.val"]),
             caption=paste( "Best enriched pathways for", samp_name, "vs", ref_name),
             format="html", format.args=list( digits=4, scientific=TRUE)),
          bootstrap_options = "striped", full_width = F
      )
    )
    
#    print( xtable( data.frame( qval=res$greater[ 1:10, "q.val"]), display=c("s", "e"),
#           caption=paste( "Best enriched pathways for", samp_name, "vs", ref_name)))
    
    interesting_pathways_df = do.call( rbind, lapply( INTERESTING_PATHWAYS, function( pathname){
      path_index = which( grepl( pathname, row.names( res$greater)))
      if( length( path_index) == 1){
        name = as.character( row.names( res$greater)[ path_index])
        qval = res$greater[ path_index, "q.val"]
        signif = (qval < 0.05)
        return( data.frame( name = name, qval = qval, signif = signif, stringsAsFactors = FALSE))
      }
    }))
  
    print( kable_styling( 
        kable( interesting_pathways_df,
               caption=paste( "Interesting pathways enrichment status for", samp_name, "vs", ref_name),
               format="html", format.args=list( digits=4, scientific=TRUE)),
        bootstrap_options = "striped", full_width = F
      )
    )

  }
}