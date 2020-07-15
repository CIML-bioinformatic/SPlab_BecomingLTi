# ###########################################################
# This script aims to produce DEG analysis between pair
# of subpopulations using the global normalization
# previously executed
# ###########################################################

## @knitr phase2_3_deseq

deg23_gene_inf = list()
deg23_gene_sup = list()
widget_list = list()
count_widget = 1
control_experiment_set = c( "EP2", "FLP2")
treatment_experiment_set = c( "EP3", "FLP3")
dsHTSeq_result_list = list()
for( i in 1:length( control_experiment_set)){

    # Get the name of the control and treatment experiments
    control_experiment = control_experiment_set[ i]
    treatment_experiment = treatment_experiment_set[ i]

    # Compute statistics on DEG between the two selected conditions
    results <- results(dsHTSeq, contrast = c( "population", treatment_experiment, control_experiment))
    rlog_results <- lfcShrink(dds=dsHTSeq, contrast=c("population",treatment_experiment,control_experiment), res=results)
    ordered_rlog_results = rlog_results[ order( rlog_results$padj, decreasing = FALSE),]
    
    dsHTSeq_result_list[[ paste0( control_experiment,"-",treatment_experiment)]] = rlog_results

    # Determine the "Best" DE genes by using DE p-value and difference of expression log2 fold change
    genes_vector = row.names( ordered_rlog_results)
    mean_diff_set = ordered_rlog_results$log2FoldChange
    de_pval_set_adjustMultiComparison = -log10( ordered_rlog_results$padj)

    mean_thresholds_2 = c( mean( mean_diff_set) - sd( mean_diff_set), mean( mean_diff_set) + sd( mean_diff_set))
    mean_index_list_2_inf = which( mean_diff_set <= mean_thresholds_2[1])
    mean_index_list_2_sup = which( mean_diff_set >= mean_thresholds_2[2])
    mean_index_list_2 = union( mean_index_list_2_inf, mean_index_list_2_sup)
    pval_index_list = which( de_pval_set_adjustMultiComparison >= -log10(0.05))
    final_index_list_2 = intersect( mean_index_list_2, pval_index_list)

    # Display a volcano plot of the analysis (if there is not too many genes)
    volcano_plot_data_frame=data.frame(row.names=genes_vector, MeanDiffSet=mean_diff_set, DEpVal=de_pval_set_adjustMultiComparison)
    volcano_plot_data_frame$Significant.DE.Zscore <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_2], TRUE, FALSE)
    best_genes_indexes = which( volcano_plot_data_frame$Significant.DE.Zscore == TRUE)[1:(BEST_NUMBER)]
    volcano_plot_data_frame$Best=rep( FALSE, nrow( volcano_plot_data_frame))
    volcano_plot_data_frame$Best[ best_genes_indexes] = TRUE

    deg23_gene_inf[[ i]] = intersect( rownames(volcano_plot_data_frame)[ mean_index_list_2_inf], row.names( volcano_plot_data_frame)[ best_genes_indexes])
    deg23_gene_sup[[ i]] = intersect( rownames(volcano_plot_data_frame)[ mean_index_list_2_sup], row.names( volcano_plot_data_frame)[ best_genes_indexes])

    volcano_plot_ggplot_base <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = DEpVal)) +
              xlab( paste( "Diff of exprs means between", POPULATION_NAME_MAPPING[ control_experiment], "and", POPULATION_NAME_MAPPING[ treatment_experiment])) +
              ylab("-log10 p.value")+
              ggtitle(paste( POPULATION_NAME_MAPPING[control_experiment],"(overexprs left) vs", POPULATION_NAME_MAPPING[ treatment_experiment],"(overexprs right)"))+
              geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
              geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
              geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
              scale_size_manual(values = c(3,1.5), guide="none") +
              geom_point(aes(color = !Significant.DE.Zscore),size=1) +
              scale_color_manual(values = c("firebrick2", "gray70"), name="Significant genes") +
              theme( plot.title = element_text( face="bold", size=12),
                     aspect.ratio=1,
                     legend.position="none") +
              scale_x_continuous(limits=c(-max(abs(volcano_plot_data_frame$MeanDiffSet)),max(abs(volcano_plot_data_frame$MeanDiffSet)) )) +
              scale_y_continuous(expand=c(0,0),limits=c(0,max(abs(volcano_plot_data_frame$DEpVal))+0.5)) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.line = element_line(colour = "black")) +
              theme( axis.text.x = element_text( size = AXIS_TICK_SIZE),
                     axis.text.y = element_text( size = AXIS_TICK_SIZE)) +
              theme( axis.title.x = element_text( size = AXIS_TITLE_SIZE),
                     axis.title.y = element_text( size = AXIS_TITLE_SIZE))+ 
              theme( legend.text = element_text( size=LEGEND_TEXT_SIZE))
    
    volcano_plot_ggplot = volcano_plot_ggplot_base + 
              geom_text_repel(
                data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Best==TRUE),
                aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Best==TRUE))),
                size = VOLCANO_PLOT_GENE_TEXT_SIZE, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
                box.padding = unit(0.25, "lines"),
                point.padding = unit(0.25, "lines")) +
              geom_text_repel(
                data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Best==TRUE),
                aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Best==TRUE))),
                size = VOLCANO_PLOT_GENE_TEXT_SIZE, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
                box.padding = unit(0.25, "lines"),
                point.padding = unit(0.25, "lines"))

    volcano_plot_ggplot_observed = volcano_plot_ggplot_base + 
              geom_text_repel(
                data = volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0 & row.names( volcano_plot_data_frame) %in% OBSERVED_GENES),],
                aes(label = rownames( volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0 & row.names( volcano_plot_data_frame) %in% OBSERVED_GENES),])),
                size = VOLCANO_PLOT_GENE_TEXT_SIZE, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
                box.padding = unit(0.25, "lines"),
                point.padding = unit(0.25, "lines")) +
              geom_text_repel(
                data = volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0 & row.names( volcano_plot_data_frame) %in% OBSERVED_GENES),],
                aes(label = rownames( volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0 & row.names( volcano_plot_data_frame) %in% OBSERVED_GENES),])),
                size = VOLCANO_PLOT_GENE_TEXT_SIZE, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
                box.padding = unit(0.25, "lines"),
                point.padding = unit(0.25, "lines"))
    
    print( volcano_plot_ggplot)
    print( volcano_plot_ggplot_observed)
    
    # Display the table of the first best DE genes
    display_columns = which( CONDITION_NAME_COLUMN_SET == control_experiment | CONDITION_NAME_COLUMN_SET == treatment_experiment)
    display_df = comparison_df[ row.names( volcano_plot_data_frame)[which( volcano_plot_data_frame$Best==TRUE)], display_columns]
    display_df$Log2FoldChange = volcano_plot_data_frame[ row.names( display_df), 1]
    
    widget_list[[ count_widget]] = htmltools::tags$h5( paste( "Comparing treatment=", treatment_experiment, " to control=", control_experiment))
    count_widget = count_widget + 1
    widget_list[[ count_widget]] = datatable( display_df, caption = paste( "\nBest", BEST_NUMBER,"DEG between", control_experiment, "and", treatment_experiment))
    count_widget = count_widget + 1

}

widget_list[[ count_widget]] = htmltools::tags$h5( "Genes commonly disregulated between phase 2 and phase 3 (in top 100 best disregulated")
count_widget = count_widget + 1

inf_genes_in_both = intersect( deg23_gene_inf[[ 1]], deg23_gene_inf[[ 2]])
sup_genes_in_both = intersect( deg23_gene_sup[[ 1]], deg23_gene_sup[[ 2]])

cat("<BR>Number of genes underexpressed in both phase 3 respect to phase 2:", length( inf_genes_in_both))
cat("<BR>Number of genes overexpressed in both phase 3 respect to phase 2:", length( sup_genes_in_both))

widget_list[[ count_widget]] = datatable( data.frame( gene_name = row.names( NORMALIZED_COUNTS[ inf_genes_in_both, ])), caption= "Genes underexpressed in both phase 3 respect to phase 2")
count_widget = count_widget + 1
widget_list[[ count_widget]] = datatable( data.frame( gene_name = row.names( NORMALIZED_COUNTS[ sup_genes_in_both, ])), caption= "Genes overexpressed in both phase 3 respect to phase 2")
count_widget = count_widget + 1

print( htmltools::tagList( widget_list))

