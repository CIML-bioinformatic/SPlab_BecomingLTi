# ###########################################################
# This script aims to produce DEG analysis between pair
# of subpopulations using the global normalization
# previously executed
# ###########################################################

## @knitr pairwise_deseq

BEST_DE_GENES_SET = vector()
widget_list = list()
count_widget = 1
for( i in 1:(length( EXPERIMENT_ID_SET)-1)){

  for( j in (i+1): length( EXPERIMENT_ID_SET)){

    # Get the name of the control and treatment experiments
    control_experiment = EXPERIMENT_ID_SET[ i]
    treatment_experiment = EXPERIMENT_ID_SET[ j]
    # widget_list[[ count_widget]] = paste( "<H5>Analyzing DEG on",control_experiment,"versus",treatment_experiment,"</H5>")
    # count_widget = count_widget + 1
    # cat("<H5>Analyzing DEG on",control_experiment,"versus",treatment_experiment,"</H5>")

    # Compute statistics on DEG between the two selected conditions
    results <- results(dsHTSeq, contrast = c( "population", treatment_experiment, control_experiment))
    rlog_results <- lfcShrink(dds=dsHTSeq, contrast=c("population",treatment_experiment,control_experiment), res=results)
    ordered_rlog_results = rlog_results[ order( rlog_results$padj, decreasing = FALSE),]

    # Determine the "Best" DE genes by using DE p-value and difference of expression log2 fold change
    genes_vector = row.names( ordered_rlog_results)
    mean_diff_set = ordered_rlog_results$log2FoldChange
    de_pval_set_adjustMultiComparison = -log10( ordered_rlog_results$padj)

    mean_thresholds = quantile(mean_diff_set, c(.05, .95))
    mean_thresholds_2 = c( mean( mean_diff_set) - sd( mean_diff_set), mean( mean_diff_set) + sd( mean_diff_set))
    mean_index_list = which( mean_diff_set <= mean_thresholds[1] | mean_diff_set >= mean_thresholds[2])
    mean_index_list_2 = which( mean_diff_set <= mean_thresholds_2[1] | mean_diff_set >= mean_thresholds_2[2])
    pval_index_list = which( de_pval_set_adjustMultiComparison >= -log10(0.05))
    final_index_list = intersect( mean_index_list, pval_index_list)
    final_index_list_2 = intersect( mean_index_list_2, pval_index_list)

    # Display a volcano plot of the analysis (if there is not too many genes)
    volcano_plot_data_frame=data.frame(row.names=genes_vector, MeanDiffSet=mean_diff_set, DEpVal=de_pval_set_adjustMultiComparison)
    volcano_plot_data_frame$Significant.DE.IQR <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list], TRUE, FALSE)
    volcano_plot_data_frame$Significant.DE.Zscore <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_2], TRUE, FALSE)
    best_genes_indexes = which( volcano_plot_data_frame$Significant.DE.Zscore == TRUE)[1:BEST_NUMBER]
    volcano_plot_data_frame$Best=rep( FALSE, nrow( volcano_plot_data_frame))
    volcano_plot_data_frame$Best[ best_genes_indexes] = TRUE

    BEST_DE_GENES_SET = append( BEST_DE_GENES_SET, row.names( volcano_plot_data_frame)[ best_genes_indexes])
    
    # volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = DEpVal)) +
    #           xlab( paste( "Diff of exprs means between", POPULATION_NAME_MAPPING[ control_experiment], "and", POPULATION_NAME_MAPPING[ treatment_experiment])) +
    #           ylab("-log10 p.value")+
    #           ggtitle(paste(control_experiment,"(overexprs left) vs",treatment_experiment,"(overexprs right)"))+
    #           geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
    #           geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
    #           geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
    #           scale_size_manual(values = c(3,1.5), guide="none") +
    #           geom_point(aes(color = !Significant.DE.Zscore),size=1) +
    #           scale_color_manual(values = c("firebrick2", "gray70"), name="Significant genes") +
    #           theme( plot.title = element_text( face="bold", size=12),
    #                  aspect.ratio=1,
    #                  legend.position="none") +
    #           geom_text_repel(
    #             data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Best==TRUE),
    #             aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Best==TRUE))),
    #             size = 3, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
    #             box.padding = unit(0.25, "lines"),
    #             point.padding = unit(0.25, "lines")) +
    #           geom_text_repel(
    #             data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Best==TRUE),
    #             aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Best==TRUE))),
    #             size = 3, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
    #             box.padding = unit(0.25, "lines"),
    #             point.padding = unit(0.25, "lines")) +
    #         scale_x_continuous(limits=c(-max(abs(volcano_plot_data_frame$MeanDiffSet)),max(abs(volcano_plot_data_frame$MeanDiffSet)) )) +
    #           scale_y_continuous(expand=c(0,0),limits=c(0,max(abs(volcano_plot_data_frame$DEpVal))+0.5))

    
    # widget_list[[ count_widget]] = volcano_plot_ggplot
    # count_widget = count_widget + 1
    #print(volcano_plot_ggplot)

    # Display the table of the first best DE genes
    display_columns = which( CONDITION_NAME_COLUMN_SET == control_experiment | CONDITION_NAME_COLUMN_SET == treatment_experiment)
    display_df = comparison_df[ row.names( volcano_plot_data_frame)[which( volcano_plot_data_frame$Best==TRUE)], display_columns]
    display_df$Log2FoldChange = volcano_plot_data_frame[ row.names( display_df), 1]
    # print( xtable( display_df, caption = paste( "\nBest", BEST_NUMBER,"differentially expressed genes")), type="latex", comment=FALSE, scalebox='0.75', caption.placement = 'top')
    widget_list[[ count_widget]] = datatable( display_df, caption = paste( "\nBest", BEST_NUMBER,"DEG between", control_experiment, "and", treatment_experiment))
    count_widget = count_widget + 1
  }

}

print( tagList( widget_list))

BEST_DE_GENES_SET = unique( BEST_DE_GENES_SET)

cat("<H5>List of all best differentially expressed genes<H5>")
datatable( NORMALIZED_COUNTS[ BEST_DE_GENES_SET, ], caption= "List of all best differentially expressed genes")
