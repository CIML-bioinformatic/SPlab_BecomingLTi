# ############################################################
# This script aims to look at pathways and GO terms enrichment
# ############################################################
# Adapted from http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2#gene-set-enrichment-analysis

# To use this script you need to define
# DESEQ_RESULT : the result of e DEseq2 analysis with log2 shrinkage as for instance 
#    DESEQ_RESULT = lfcShrink( results( DESeq( dds), contrast = <some constrast> )), ......)
# The first 3 lines of the script must be checked to see if they match with the used DEseq results columns/rown names

## @knitr gsea_analysis_with_deseq2

#
# Adapt the DESeq result to have columns with entrez-ids and gene symbols
#
# -- add entrez_id of genes to the results
DESEQ_RESULT$hgnc_symbol <- row.names( DESEQ_RESULT)
# -- add entrez_id of genes to the results
DESEQ_RESULT$entrezid <- convertIDs( row.names( DESEQ_RESULT), "SYMBOL", "ENTREZID", org.Mm.eg.db )
# -- filter out the genes without entrez_id of without adjusted p-values
res2 <- DESEQ_RESULT[ DESEQ_RESULT$entrezid %in% keys( reactome.db, "ENTREZID" ) & !is.na( DESEQ_RESULT$padj ) , ]

# Test if the used database is one that is available for analysis
USED_DATABASE =  match.arg( USED_DATABASE, c( "REACTOME", "GO-CC", "GO-MF", "GO-BP", "KEGG"))

#
# Prepare geneset table for REACTOME DB
#
if( USED_DATABASE == "REACTOME"){
  # Build a table mapping the gene entrez_id to the Reactome ID it is part of
  reactomeTable <- AnnotationDbi::select( reactome.db,
                                          keys=as.character( res2$entrezid), keytype="ENTREZID",
                                          columns=c("ENTREZID","REACTOMEID") )
  # Build a table indicating if e gene entrez_id is part of the reactome
  incm <- do.call( rbind, with(reactomeTable, tapply( 
    ENTREZID, factor(REACTOMEID), function(x) res2$entrezid %in% x ) ))
  colnames(incm) <- res2$entrez
  # Build a map of reactome ID to name
  map_id_to_name = do.call( rbind, lapply( unique( reactomeTable$REACTOMEID), function( reactomeID){
    reactome_name = reactomePATHID2NAME[[ reactomeID]]
    reactome_name = ifelse( is.null( reactome_name), "NA", reactome_name)
    data.frame( ID=reactomeID, NAME=reactome_name, stringsAsFactors = FALSE)
  }))
  row.names( map_id_to_name) = map_id_to_name$ID
}

#
# Prepare geneset table for GO DB
#
if( USED_DATABASE == "GO-CC" || USED_DATABASE == "GO-BP" || USED_DATABASE == "GO-MF"){

  # Identify the ontology to use from the database name
  ontology = gsub( "GO-", "", USED_DATABASE, fixed=TRUE)

  # Build a table mapping the gene entrez_id to the Reactome ID it is part of
  goTable <- AnnotationDbi::select( org.Mm.eg.db, 
                                    keys=as.character( res2$entrezid), keytype="ENTREZID", 
                                    columns=c("ENTREZID","GO"))
  
  goTable = goTable[ which( goTable$ONTOLOGY == ontology), ]
  
  # Build a table indicating if e gene entrez_id is part of the reactome
  incm <- do.call( rbind, with( goTable, tapply( 
                                ENTREZID, factor( GO), function(x) res2$entrezid %in% x ) ))
  colnames(incm) <- res2$entrez
  # Build a map of GO ID to GO term
  map_id_to_name = select(GO.db, keys=goTable$GO, columns=c("TERM", "ONTOLOGY"), keytype="GOID")
  map_id_to_name = map_id_to_name[ which( map_id_to_name$ONTOLOGY == ontology), c( "GOID", "TERM")]
  map_id_to_name = map_id_to_name[ !duplicated( map_id_to_name$GOID), ]
  names( map_id_to_name) = c( "ID", "NAME")
  row.names( map_id_to_name) = map_id_to_name$ID
}

#
# Test the various genesets of the chosen database against the DEseq results
#

# Remove all genesets with less than 20 genes and more than 200
within <- function(x, lower, upper) (x>=lower & x<=upper)
incm <- incm[ within(rowSums(incm), lower=20, upper=200), ]

# Compute the significance of difference of expression of genes in the geneset respec to 0 (in log2FoldChange)
# using a t-test
gseaResult <- do.call( rbind, lapply( rownames(incm), function( genesetID ) {
                            isMember <- incm[ genesetID, ]
                            geneset_name = map_id_to_name[ genesetID, "NAME"]
                            geneset_name = ifelse( is.null( geneset_name), "NA", geneset_name)
                            data.frame( 
                              genesetID   = genesetID,
                              genesetName = geneset_name,
                              numGenes    = sum( isMember),
                              avgLFC      = mean( res2$log2FoldChange[isMember] ),
                              sdLFC       = sd( res2$log2FoldChange[isMember] ),
                              zValue      = mean( res2$log2FoldChange[isMember] ) /sd( res2$log2FoldChange[isMember] ),
                              strength    = sum( res2$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
                              pvalue      = t.test( res2$log2FoldChange[ isMember ] )$p.value,
                              geneNames   = paste( res2$hgnc_symbol[isMember], collapse="/"),
                              stringsAsFactors = FALSE 
                            ) 
                          }))

# Adjust with Benjamini-Hochberg procedure the p-values and reorder column for final result
gseaResult$padjust <- p.adjust( gseaResult$pvalue, "BH" )
gseaResult = gseaResult[ , c( "genesetID", "genesetName", "numGenes", "avgLFC", "sdLFC", "zValue", "strength", "pvalue", "padjust", "geneNames")]

# Extract the significant information at 5% FDR
gseaResultSignif <- gseaResult[ gseaResult$padjust < 0.05, ]

# Display the result
datatable( gseaResultSignif[ , c( "genesetID", "genesetName", "numGenes", "avgLFC", "sdLFC", "zValue", "strength", "pvalue", "padjust")], caption = paste( "Significant genesets in", USED_DATABASE))

# Export result to file
write.table( gseaResultSignif, 
             file= file.path( OUTPUT_DIR, paste0( control_experiment_names, "_vs_", treatment_experiment_names, "_GSEA_DEseq2_", USED_DATABASE, ".txt")), 
                                        sep=",", quote = TRUE, row.names = FALSE, col.names = TRUE)

