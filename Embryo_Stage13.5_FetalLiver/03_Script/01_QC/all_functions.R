# ############################################################
# All the functions which will be used in a scRNA-Seq Analysis
# ############################################################

## @knitr load_all_functions


library(gplots)
library(ggplot2)
library(gridExtra)
library(cowplot)

#pca
library( ade4)

# t-SNE (bh-SNE)
library(Rtsne)

# biclustering
library(pheatmap)

#pseudotime DDRTree
#library(monocle)

#volcano plot
#library(MAST) # MAST not present in the docker scrnaseq_seurat3_scran ... # via bioconductor
library(data.table)
library(plyr)
library(ggrepel)

# for Kernel Smooth
library(KernSmooth)


TWO_COLORS_VECTOR=c("#FBE9E7", "#FF3D00")
TWO_COLORS_GRADIANT=c("grey80",colorRampPalette(TWO_COLORS_VECTOR)(n = 2000))
color_NA_geom_point="grey50"

col_heatmap = c("grey",colorRampPalette(c("white",  "orange"))(n = 299))

THREE_COLORS_HEATMAP = colorRampPalette(c("blue", "white", "red"))(n = 299)

# #####################################################
# #####################################################
# #####################################################
# beautiful colors (same as in ggplot2)
# #####################################################
# #####################################################
# #####################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 200)[1:n]
}

# #####################################################
# #####################################################
# #####################################################
# pca function which generate only the histogram
# with the percentage of information per PC
# #####################################################
# #####################################################
# #####################################################
pca_hist_function <- function( raw_data_for_pca, PlotHist = FALSE, VectorToColorPCA = NULL, NameOfVectorToColor = "", 
                               dim_all=FALSE, dim_x="PC1", dim_y="PC2", bool_facet_wrap=FALSE, VectorFacetWrap=NA, title="",bool_scalePCA=FALSE, bool_centerPCA=TRUE,
                               bool_plot=TRUE){
  limits_allFigs=c(0,27)
  
  if(length(VectorFacetWrap)<=1){
    VectorFacetWrap = VectorToColorPCA
  }
  
  # special case if we want to color the pca by gene expression (continuous variable)
  if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
    VectorToColorPCA[which(VectorToColorPCA==0)]=NA
  }
  
  # raw_data_for_pca must have only the Gene Expression and no additional columns
  
  # Execute a Centered pca with ade4pca_function
  # ############################################
  result__pca <- dudi.pca( raw_data_for_pca, scannf = FALSE, nf = ncol(raw_data_for_pca), scale=bool_scalePCA, center=bool_centerPCA)
  
  # Generate the histogram ggplot with the percentage of information for each principal composants
  # ##############################################################################################
  percentageInfo_byPC_df = data.frame(PC=seq(1,length(result__pca$eig),1), PercentInfo=100*result__pca$eig/sum( result__pca$eig))
  
  if(PlotHist==TRUE){
    percentage_information_histogram <- ggplot(percentageInfo_byPC_df, aes(PC, PercentInfo)) +
      geom_bar(stat="identity", color="black", fill="gray90") +
      xlab("Principal Components") +
      ylab("Percentage of variants") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0,max(percentageInfo_byPC_df$PercentInfo)+2))
    print(percentage_information_histogram)
  }
  
  # #####################################
  # Information for the Sweave file
  # #####################################
  if(PlotHist==TRUE){
    cat("\n****************************************************")
    cat("\n Centered PCA:")
    cat("\n****************************************************")  
  }
  
  # Look at the most informative genes along the axes
  threshold_norm = 0.75
  threshold_percent = 0.05 # 0.05
  names_best_on_axis_all_axes = vector()
  for( cs_index in c(1,2,3)){ # c(1,2,3)
    # Compute on the axis, the best variables, i.e. the variable which have the max coordinates on that axis
    current_vector = result__pca$c1[[cs_index]]
    number_best = floor( length( current_vector) * threshold_percent)
    decreasing_order = order( abs(current_vector), decreasing = TRUE)
    names_best_on_axis = row.names( result__pca$c1)[decreasing_order][1: number_best]
    names_best_on_axis_all_axes = unique( append( names_best_on_axis_all_axes, names_best_on_axis))
    
    if(PlotHist==TRUE){
      cat("\n\nOn CS", cs_index,",", threshold_percent*100, "% of best variables are\n")
      print( names_best_on_axis)
      
      # Look at the variables with the maximal cumulative information on each axes
      cat("\nOn CS", cs_index, ", direct and cummulative information carried by best variables:\n")
      perc_info_vector=c() 
      genes_names=c() 
      perc_info_sum = 0
      for( gene_index in seq(1, length( current_vector))){
        perc_info = result__pca$c1[[ cs_index]][decreasing_order][gene_index]^2
        perc_info_sum = perc_info_sum + perc_info
        cat("\n", row.names( result__pca$c1)[decreasing_order][gene_index], "\t:\t", signif((perc_info)*100, 3), "%\t:\t", signif( perc_info_sum*100 , 3), " %")
        if( perc_info_sum >= threshold_norm){
          break 
        }
      }
    }
  }
  
  
  # ##########
  # plot PCA
  # ##########
  if(!is.null(VectorToColorPCA)){
    raw_data_for_pca$DataToColorPCA = VectorToColorPCA
    # take the Principal Composants Coordinates
    # ########################################
    raw_data_for_pca$PC1=result__pca$li$Axis1
    raw_data_for_pca$PC2=result__pca$li$Axis2
    raw_data_for_pca$PC3=result__pca$li$Axis3
    raw_data_for_pca$PC4=result__pca$li$Axis4
    raw_data_for_pca$PC5=result__pca$li$Axis5
    raw_data_for_pca$PC6=result__pca$li$Axis6
    pca_percentageInfoByAxis = as.integer(100*result__pca$eig/sum( result__pca$eig)+0.5)
    names(pca_percentageInfoByAxis)=paste0("PC",1:length(pca_percentageInfoByAxis))
    
    if(dim_all==FALSE){
      if(title!=""){
        title=paste(title,"\n")
      }
      PCX=dim_x
      PCY=dim_y
      # we change the scale of the plot's axis to have a perfect plot 
      distPCX= max(raw_data_for_pca[,PCX])-min(raw_data_for_pca[,PCX])
      distPCY= max(raw_data_for_pca[,PCY])-min(raw_data_for_pca[,PCY])
      addIntervalPCX=0
      addIntervalPCY=0
      if(distPCX>distPCY){
        addIntervalPCY=(distPCX-distPCY)/2
      } else {
        addIntervalPCX=(distPCY-distPCX)/2
      }
      if(bool_facet_wrap == FALSE){
        if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point() +
            coord_fixed() +
            scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA colored by",NameOfVectorToColor))
          print(plot_pca)
        }else{
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point() +
            coord_fixed() +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(title,paste("PCA colored by",NameOfVectorToColor))
          print(plot_pca)
        }
      }else{
        # we split the plot into multiple plots
        #raw_data_for_pca$DataToColorPCA2 = raw_data_for_pca$DataToColorPCA
        raw_data_for_pca$DataToFacetWrap = VectorFacetWrap
        raw_data_for_pca$DataToFacetWrap2 =raw_data_for_pca$DataToFacetWrap
        if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point(data=raw_data_for_pca[,-match("DataToFacetWrap" , colnames(raw_data_for_pca))],aes_string(x=PCX,y=PCY,group="DataToFacetWrap2"),colour="grey")+
            geom_point() +
            coord_fixed() +
            scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA splited by",NameOfVectorToColor)) + 
            #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none") +
            facet_wrap(~DataToFacetWrap)
          print(plot_pca)
        }else{
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point(data=raw_data_for_pca[,-match("DataToFacetWrap" , colnames(raw_data_for_pca))],aes_string(x=PCX,y=PCY,group="DataToFacetWrap2"),colour="grey")+
            geom_point() +
            coord_fixed() +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA splited by",NameOfVectorToColor)) + 
            #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none") +
            facet_wrap(~DataToFacetWrap)
          print(plot_pca)
        }
      }
      
    } else {
      ggplot_list = list()
      PC_v = c("PC1","PC2","PC3")
      index=0
      for(i in 1:(length(PC_v)-1)){
        for(j in (i+1):length(PC_v)){
          index=index+1
          gg_title=""
          if(index==1){
            gg_title=title
          }
          PCX=PC_v[i]
          PCY=PC_v[j]
          # we change the scale of the plot's axis to have a perfect plot 
          distPCX= max(raw_data_for_pca[,PCX])-min(raw_data_for_pca[,PCX])
          distPCY= max(raw_data_for_pca[,PCY])-min(raw_data_for_pca[,PCY])
          addIntervalPCX=0
          addIntervalPCY=0
          if(distPCX>distPCY){
            addIntervalPCY=(distPCX-distPCY)/2
          } else {
            addIntervalPCX=(distPCY-distPCX)/2
          }
          if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
            plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
              geom_point(size=0.6) +
              coord_fixed() +
              scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
              xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
              ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
              xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
              ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
              ggtitle(gg_title) +
              guides(colour=guide_legend(title=NameOfVectorToColor))
            #print(plot_pca)
            ggplot_list[[index]]=plot_pca
          }else{
            plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
              geom_point(size=0.6) +
              coord_fixed() +
              xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
              ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
              xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
              ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
              ggtitle(gg_title) +
              guides(colour=guide_legend(title=NameOfVectorToColor))
            #print(plot_pca)
            ggplot_list[[index]]=plot_pca
          }
        }
      }
      
      #gene_pca_list=pca_gene_information_function(result__pca)
      #plots_list=c(gene_pca_list[1],ggplot_list[1],gene_pca_list[2],ggplot_list[2],gene_pca_list[3],ggplot_list[3])
      #do.call("grid.arrange", c(plots_list, nrow=3,ncol=2))
      
      if(bool_plot){
        do.call("grid.arrange", c(ggplot_list, nrow=3))
      }
      #do.call("grid.arrange", c(ggplot_list, nrow=3))
    }
    
  }
  
  return(list(pca_info = result__pca, genes_information = names_best_on_axis_all_axes))
}


# #####################################################
# #####################################################
# #####################################################
# pca function which generate only the histogram
# with the percentage of information per PC
# MarkdownVersion: pas de grid.plot
# #####################################################
# #####################################################
# #####################################################
pca_hist_function_RmdVersion <- function( raw_data_for_pca, PlotHist = FALSE, VectorToColorPCA = NULL, NameOfVectorToColor = "", 
                               dim_all=FALSE, dim_x="PC1", dim_y="PC2", bool_facet_wrap=FALSE, VectorFacetWrap=NA, title="",bool_scalePCA=FALSE, bool_centerPCA=TRUE){
  limits_allFigs=c(0,27)
  
  if(length(VectorFacetWrap)<=1){
    VectorFacetWrap = VectorToColorPCA
  }
  
  # special case if we want to color the pca by gene expression (continuous variable)
  if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
    VectorToColorPCA[which(VectorToColorPCA==0)]=NA
  }
  
  # raw_data_for_pca must have only the Gene Expression and no additional columns
  
  # Execute a Centered pca with ade4pca_function
  # ############################################
  result__pca <- dudi.pca( raw_data_for_pca, scannf = FALSE, nf = ncol(raw_data_for_pca), scale=bool_scalePCA, center=bool_centerPCA)
  
  # Generate the histogram ggplot with the percentage of information for each principal composants
  # ##############################################################################################
  percentageInfo_byPC_df = data.frame(PC=seq(1,length(result__pca$eig),1), PercentInfo=100*result__pca$eig/sum( result__pca$eig))
  
  if(PlotHist==TRUE){
    percentage_information_histogram <- ggplot(percentageInfo_byPC_df, aes(PC, PercentInfo)) +
      geom_bar(stat="identity", color="black", fill="gray90") +
      xlab("Principal Components") +
      ylab("Percentage of variants") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0,max(percentageInfo_byPC_df$PercentInfo)+2))
    print(percentage_information_histogram)
  }

  # #####################################
  # Information for the Sweave file
  # #####################################
  if(PlotHist==TRUE){
    cat("\n****************************************************")
    cat("\n Centered PCA:")
    cat("\n****************************************************")
  }

  # Look at the most informative genes along the axes
  threshold_norm = 0.75
  threshold_percent = 0.05 # 0.05
  names_best_on_axis_all_axes = vector()
  for( cs_index in c(1,2,3)){ # c(1,2,3)
    # Compute on the axis, the best variables, i.e. the variable which have the max coordinates on that axis
    current_vector = result__pca$c1[[cs_index]]
    number_best = floor( length( current_vector) * threshold_percent)
    decreasing_order = order( abs(current_vector), decreasing = TRUE)
    names_best_on_axis = row.names( result__pca$c1)[decreasing_order][1: number_best]
    names_best_on_axis_all_axes = unique( append( names_best_on_axis_all_axes, names_best_on_axis))
    
    if(PlotHist==TRUE){
      cat("\n\nOn CS", cs_index,",", threshold_percent*100, "% of best variables are\n")
      print( names_best_on_axis)
      
      # Look at the variables with the maximal cumulative information on each axes
      cat("\nOn CS", cs_index, ", direct and cummulative information carried by best variables:\n")
      perc_info_vector=c() 
      genes_names=c() 
      perc_info_sum = 0
      for( gene_index in seq(1, length( current_vector))){
        perc_info = result__pca$c1[[ cs_index]][decreasing_order][gene_index]^2
        perc_info_sum = perc_info_sum + perc_info
        cat("\n", row.names( result__pca$c1)[decreasing_order][gene_index], "\t:\t", signif((perc_info)*100, 3), "%\t:\t", signif( perc_info_sum*100 , 3), " %")
        if( perc_info_sum >= threshold_norm){
          break 
        }
      }
    }
  }
  
  
  # ##########
  # plot PCA
  # ##########
  if(!is.null(VectorToColorPCA)){
    raw_data_for_pca$DataToColorPCA = VectorToColorPCA
    # take the Principal Composants Coordinates
    # ########################################
    raw_data_for_pca$PC1=result__pca$li$Axis1
    raw_data_for_pca$PC2=result__pca$li$Axis2
    raw_data_for_pca$PC3=result__pca$li$Axis3
    raw_data_for_pca$PC4=result__pca$li$Axis4
    raw_data_for_pca$PC5=result__pca$li$Axis5
    raw_data_for_pca$PC6=result__pca$li$Axis6
    pca_percentageInfoByAxis = as.integer(100*result__pca$eig/sum( result__pca$eig)+0.5)
    names(pca_percentageInfoByAxis)=paste0("PC",1:length(pca_percentageInfoByAxis))
    
    if(dim_all==FALSE){
      if(title!=""){
        title=paste(title,"\n")
      }
      PCX=dim_x
      PCY=dim_y
      # we change the scale of the plot's axis to have a perfect plot 
      distPCX= max(raw_data_for_pca[,PCX])-min(raw_data_for_pca[,PCX])
      distPCY= max(raw_data_for_pca[,PCY])-min(raw_data_for_pca[,PCY])
      addIntervalPCX=0
      addIntervalPCY=0
      if(distPCX>distPCY){
        addIntervalPCY=(distPCX-distPCY)/2
      } else {
        addIntervalPCX=(distPCY-distPCX)/2
      }
      if(bool_facet_wrap == FALSE){
        if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point() +
            coord_fixed() +
            scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA colored by",NameOfVectorToColor))
          print(plot_pca)
        }else{
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point() +
            coord_fixed() +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(title,paste("PCA colored by",NameOfVectorToColor))
          print(plot_pca)
        }
      }else{
        # we split the plot into multiple plots
        #raw_data_for_pca$DataToColorPCA2 = raw_data_for_pca$DataToColorPCA
        raw_data_for_pca$DataToFacetWrap = VectorFacetWrap
        raw_data_for_pca$DataToFacetWrap2 =raw_data_for_pca$DataToFacetWrap
        if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point(data=raw_data_for_pca[,-match("DataToFacetWrap" , colnames(raw_data_for_pca))],aes_string(x=PCX,y=PCY,group="DataToFacetWrap2"),colour="grey")+
            geom_point() +
            coord_fixed() +
            scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA splited by",NameOfVectorToColor)) + 
            #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none") +
            facet_wrap(~DataToFacetWrap)
          print(plot_pca)
        }else{
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point(data=raw_data_for_pca[,-match("DataToFacetWrap" , colnames(raw_data_for_pca))],aes_string(x=PCX,y=PCY,group="DataToFacetWrap2"),colour="grey")+
            geom_point() +
            coord_fixed() +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA splited by",NameOfVectorToColor)) + 
            #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none") +
            facet_wrap(~DataToFacetWrap)
          print(plot_pca)
        }
      }
      
    } else {
      ggplot_list = list()
      PC_v = c("PC1","PC2","PC3")
      index=0
      for(i in 1:(length(PC_v)-1)){
        for(j in (i+1):length(PC_v)){
          index=index+1
          gg_title=""
          if(index==1){
            gg_title=title
          }
          PCX=PC_v[i]
          PCY=PC_v[j]
          # we change the scale of the plot's axis to have a perfect plot 
          distPCX= max(raw_data_for_pca[,PCX])-min(raw_data_for_pca[,PCX])
          distPCY= max(raw_data_for_pca[,PCY])-min(raw_data_for_pca[,PCY])
          addIntervalPCX=0
          addIntervalPCY=0
          if(distPCX>distPCY){
            addIntervalPCY=(distPCX-distPCY)/2
          } else {
            addIntervalPCX=(distPCY-distPCX)/2
          }
          if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
            plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
              geom_point(size=0.6) +
              coord_fixed() +
              scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
              xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
              ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
              xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
              ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
              ggtitle(gg_title) +
              guides(colour=guide_legend(title=NameOfVectorToColor))
            print(plot_pca)
            ggplot_list[[index]]=plot_pca
          }else{
            plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
              geom_point(size=0.6) +
              coord_fixed() +
              xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
              ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
              xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
              ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
              ggtitle(gg_title) +
              guides(colour=guide_legend(title=NameOfVectorToColor))
            print(plot_pca)
            ggplot_list[[index]]=plot_pca
          }
        }
      }
      
      #gene_pca_list=pca_gene_information_function(result__pca)
      #plots_list=c(gene_pca_list[1],ggplot_list[1],gene_pca_list[2],ggplot_list[2],gene_pca_list[3],ggplot_list[3])
      #do.call("grid.arrange", c(plots_list, nrow=3,ncol=2))
      
      #do.call("grid.arrange", c(ggplot_list, nrow=3))
    }
    
  }
  
  return(list(pca_info = result__pca, genes_information = names_best_on_axis_all_axes))
}

# #####################################################
# #####################################################
# #####################################################
# pca function which generate plot for 2 specific PC
# #####################################################
# #####################################################
# #####################################################
PCA_specific_PC_function <- function(pca_info, PCx, PCy, VectorToColorPCA, threshold_norm=0.75, gg_title=""){
  pca=pca_info
  ggplot_list = list()
  threshold_norm_top_genes = 0.6
  
  variable_genes_withThreshold_v = c()
  #threshold_norm = 0.75
  # take the top gene which drive the PC1 and the PC2 from the pca report
  # #########################################################
  PC_i = PCx
  PC_j = PCy
  
  genes_vector = c()
  for(i in c(PC_i,PC_j)){
    current_vector = pca$c1[[i]]
    decreasing_order = order( abs(current_vector), decreasing = TRUE)
    for( gene_index in seq(1, length( current_vector))){
      current_norm =(sum( (pca$c1[[ i]][decreasing_order][1: gene_index])^2))
      genes_vector=c(genes_vector,row.names( pca$c1)[decreasing_order][gene_index])
      if( current_norm >= threshold_norm){
        break
      }
    }
  }
  genes_to_study = unique(genes_vector)
  
  pca_percentageInfoByAxis = as.integer(100*pca$eig/sum( pca$eig)+0.5)
  if(is.numeric(VectorToColorPCA)){
    VectorToColorPCA[VectorToColorPCA==0]=NA
    coord_pca_df = pca$li
    coord_pca_df$VectorToColorPCA = VectorToColorPCA
    coord_pca_plot = ggplot(coord_pca_df, aes_string(paste0("Axis",PC_i),paste0("Axis",PC_j), color="VectorToColorPCA"))+
      geom_point()+
      xlab(paste0("PC",PC_i," (",pca_percentageInfoByAxis[PC_i],"%)")) +
      ylab(paste0("PC",PC_j," (",pca_percentageInfoByAxis[PC_j],"%)")) +
      scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Gene") +
      ggtitle(gg_title)+
      coord_fixed()
    ggplot_list[[1]]=coord_pca_plot
  }else{
    coord_pca_df = pca$li
    coord_pca_df$VectorToColorPCA = VectorToColorPCA
    coord_pca_plot = ggplot(coord_pca_df, aes_string(paste0("Axis",PC_i),paste0("Axis",PC_j), color="VectorToColorPCA"))+
      geom_point()+
      xlab(paste0("PC",PC_i," (",pca_percentageInfoByAxis[PC_i],"%)")) +
      ylab(paste0("PC",PC_j," (",pca_percentageInfoByAxis[PC_j],"%)")) +
      ggtitle(gg_title)+
      coord_fixed()
    ggplot_list[[1]]=coord_pca_plot
  }
  
  
  # Create the article figure AR
  # #############################
  pca_to_use_df = pca$co[genes_to_study,]
  pca_to_use_df$Labels = rownames(pca_to_use_df)
  pca_genes_plot <- ggplot(pca_to_use_df, aes_string(paste0("Comp",PC_i),paste0("Comp",PC_j), label="Labels")) +
    geom_segment(aes_string(x=0, y=0, xend=paste0("Comp",PC_i), yend=paste0("Comp",PC_j)), arrow = arrow(length = unit(0.01, "npc")), color="grey" ) +
    #geom_text(size=3.5)
    #theme_bw(base_size = 1) + #theme(legend.position = "bottom") +
    #theme( plot.title = element_text( face="bold", size=1)) +
    geom_text_repel(segment.color = NA, size=2, fontface="italic")+#,
    #box.padding = unit(0.35, "lines"),
    #point.padding = unit(0.25, "lines")) +
    xlab(paste0("PC",PC_i," (",pca_percentageInfoByAxis[PC_i],"%)")) +
    ylab(paste0("PC",PC_j," (",pca_percentageInfoByAxis[PC_j],"%)")) +
    coord_fixed()
  #print(pca_genes_plot)
  ggplot_list[[2]]=pca_genes_plot
  
  variable_genes_withThreshold_v = c(variable_genes_withThreshold_v, rownames(pca_to_use_df))

  
  do.call("grid.arrange", c(ggplot_list, nrow=2))
  
  variable_genes_withThreshold_v = unique(variable_genes_withThreshold_v)
  #return(ggplot_list)
  return(variable_genes_withThreshold_v)
}

# #####################################################
# #####################################################
# #####################################################
# PCA genes information distribution
# #####################################################
# #####################################################
# #####################################################
top_genes_fromPCA_function <-function(pca_info, dim_min=1, dim_max=3, threshold_norm = 0.75){
  pca=pca_info
  variable_genes_withThreshold_v = c()
  #threshold_norm = 0.75
  # take the top gene which drive the PC1 and the PC2 from the pca report
  # #########################################################
  for(i in dim_min:dim_max){
    genes_vector = c()
    
    current_vector = pca$c1[[i]]
    decreasing_order = order( abs(current_vector), decreasing = TRUE)
    for( gene_index in seq(1, length( current_vector))){
      current_norm = (sum( (pca$c1[[ i]][decreasing_order][1: gene_index])^2))
      genes_vector=c(genes_vector,row.names( pca$c1)[decreasing_order][gene_index])
      if( current_norm >= threshold_norm){
        break
      }
    }
    variable_genes_withThreshold_v = c(variable_genes_withThreshold_v, unique(genes_vector) )
    
  }
  variable_genes_withThreshold_v = unique(variable_genes_withThreshold_v)
  return(variable_genes_withThreshold_v)
}

# #####################################################
# #####################################################
# #####################################################
# PCA genes information distribution
# #####################################################
# #####################################################
# #####################################################
pca_gene_information_function <-function(pca_info, dim_all=FALSE, dim_x=1, dim_y=2, threshold_norm = 0.75, title="",size_text=2){
  pca=pca_info
  ggplot_list = list()
  index=0
  threshold_norm_top_genes = 0.6
  
  variable_genes_withThreshold_v = c()
  #threshold_norm = 0.75
  # take the top gene which drive the PC1 and the PC2 from the pca report
  # #########################################################
  for(PC_i in 1:2){
    for(PC_j in (PC_i+1):3){
      index=index+1
      gg_title=""
      if(index==1){
        gg_title=title
      }
      genes_vector = c()
      for(i in c(PC_i,PC_j)){
        current_vector = pca$c1[[i]]
        decreasing_order = order( abs(current_vector), decreasing = TRUE)
        for( gene_index in seq(1, length( current_vector))){
          current_norm = (sum( (pca$c1[[ i]][decreasing_order][1: gene_index])^2))
          genes_vector=c(genes_vector,row.names( pca$c1)[decreasing_order][gene_index])
          if( current_norm >= threshold_norm){
            break
          }
        }
      }
      genes_to_study = unique(genes_vector)
      
      pca_percentageInfoByAxis = as.integer(100*pca$eig/sum( pca$eig)+0.5)
      
      # Create the article figure AR
      # #############################
      pca_to_use_df = pca$co[genes_to_study,]
      pca_to_use_df$Labels = rownames(pca_to_use_df)
      pca_genes_plot <- ggplot(pca_to_use_df, aes_string(paste0("Comp",PC_i),paste0("Comp",PC_j), label="Labels")) +
        geom_segment(aes_string(x=0, y=0, xend=paste0("Comp",PC_i), yend=paste0("Comp",PC_j)), arrow = arrow(length = unit(0.01, "npc")), color="grey" ) +
        #geom_text(size=3.5)
        #theme_bw(base_size = 1) + #theme(legend.position = "bottom") +
        #theme( plot.title = element_text( face="bold", size=1)) +
        geom_text_repel(segment.color = NA, size=size_text, fontface="italic")+#,
        #box.padding = unit(0.35, "lines"),
        #point.padding = unit(0.25, "lines")) +
        xlab(paste0("PC",PC_i," (",pca_percentageInfoByAxis[PC_i],"%)")) +
        ylab(paste0("PC",PC_j," (",pca_percentageInfoByAxis[PC_j],"%)")) +
        #coord_fixed()+
        theme(aspect.ratio=1)+
        #coord_equal()+
        ggtitle(gg_title)
      #print(pca_genes_plot)
      ggplot_list[[index]]=pca_genes_plot
      
      variable_genes_withThreshold_v = c(variable_genes_withThreshold_v, rownames(pca_to_use_df))
    }
  }
  
  if(FALSE){
    # Look at the genes with the maximal cumulative information on each axes
    threshold_norm = threshold_norm_top_genes
    genes_maxCumulativeInfo_vector = c()
    for(cs_index in c(1,2)){
      cat("\n\n Genes information on Principal Composant ",cs_index," :")
      cat("\n--------------------------------------------------")
      current_vector=pca$c1[[cs_index]]
      decreasing_order = order( abs(current_vector), decreasing = TRUE)
      perc_info_vector=c() 
      genes_names=c() 
      perc_info_sum = 0
      for( gene_index in seq(1, length( current_vector))){
        perc_info = pca$c1[[ cs_index]][decreasing_order][gene_index]^2
        perc_info_sum = perc_info_sum + perc_info
        cat("\n", row.names( pca$c1)[decreasing_order][gene_index], "\t:\t", signif((perc_info)*100, 3), "%\t:\t", signif( perc_info_sum*100 , 3), " %")
        if(! row.names( pca$c1)[decreasing_order][gene_index] %in% genes_maxCumulativeInfo_vector){
          genes_maxCumulativeInfo_vector = c(genes_maxCumulativeInfo_vector, row.names( pca$c1)[decreasing_order][gene_index])
        }
        if( perc_info_sum >= threshold_norm){
          break
        }
      }
    }
  }
  if(dim_all){
    do.call("grid.arrange", c(ggplot_list, nrow=3))
  }else{
    if(dim_x==1 & dim_y==2){
      print(ggplot_list[[1]])
    } else if(dim_x==1 & dim_y==3){
      print(ggplot_list[[2]])
    } else if(dim_x==2 & dim_y==3){
      print(ggplot_list[[3]])
    }
  }
  
  
  variable_genes_withThreshold_v = unique(variable_genes_withThreshold_v)
  #return(ggplot_list)
  return(variable_genes_withThreshold_v)
}


# #####################################################
# #####################################################
# #####################################################
# PCA genes information distribution
# MarkdownVersion: no grid.arrange
# #####################################################
# #####################################################
# #####################################################
pca_gene_information_function_RmdVersion <-function(pca_info, dim_all=FALSE, dim_x=1, dim_y=2, threshold_norm = 0.75, title="", size_text=2){
  pca=pca_info
  ggplot_list = list()
  index=0
  threshold_norm_top_genes = 0.6
  
  variable_genes_withThreshold_v = c()
  #threshold_norm = 0.75
  # take the top gene which drive the PC1 and the PC2 from the pca report
  # #########################################################
  for(PC_i in 1:2){
    for(PC_j in (PC_i+1):3){
      index=index+1
      gg_title=""
      if(index==1){
        gg_title=title
      }
      genes_vector = c()
      for(i in c(PC_i,PC_j)){
        current_vector = pca$c1[[i]]
        decreasing_order = order( abs(current_vector), decreasing = TRUE)
        for( gene_index in seq(1, length( current_vector))){
          current_norm = (sum( (pca$c1[[ i]][decreasing_order][1: gene_index])^2))
          genes_vector=c(genes_vector,row.names( pca$c1)[decreasing_order][gene_index])
          if( current_norm >= threshold_norm){
            break
          }
        }
      }
      genes_to_study = unique(genes_vector)
      
      pca_percentageInfoByAxis = as.integer(100*pca$eig/sum( pca$eig)+0.5)
      
      # Create the article figure AR
      # #############################
      pca_to_use_df = pca$co[genes_to_study,]
      pca_to_use_df$Labels = rownames(pca_to_use_df)
      pca_genes_plot <- ggplot(pca_to_use_df, aes_string(paste0("Comp",PC_i),paste0("Comp",PC_j), label="Labels")) +
        geom_segment(aes_string(x=0, y=0, xend=paste0("Comp",PC_i), yend=paste0("Comp",PC_j)), arrow = arrow(length = unit(0.01, "npc")), color="grey" ) +
        #geom_text(size=3.5)
        #theme_bw(base_size = 1) + #theme(legend.position = "bottom") +
        #theme( plot.title = element_text( face="bold", size=1)) +
        geom_text_repel(segment.color = NA, size=size_text, fontface="italic")+#,
        #box.padding = unit(0.35, "lines"),
        #point.padding = unit(0.25, "lines")) +
        xlab(paste0("PC",PC_i," (",pca_percentageInfoByAxis[PC_i],"%)")) +
        ylab(paste0("PC",PC_j," (",pca_percentageInfoByAxis[PC_j],"%)")) +
        coord_fixed()+
        ggtitle(gg_title)
      print(pca_genes_plot)
      ggplot_list[[index]]=pca_genes_plot
      
      variable_genes_withThreshold_v = c(variable_genes_withThreshold_v, rownames(pca_to_use_df))
    }
  }
  
  if(FALSE){
    # Look at the genes with the maximal cumulative information on each axes
    threshold_norm = threshold_norm_top_genes
    genes_maxCumulativeInfo_vector = c()
    for(cs_index in c(1,2)){
      cat("\n\n Genes information on Principal Composant ",cs_index," :")
      cat("\n--------------------------------------------------")
      current_vector=pca$c1[[cs_index]]
      decreasing_order = order( abs(current_vector), decreasing = TRUE)
      perc_info_vector=c() 
      genes_names=c() 
      perc_info_sum = 0
      for( gene_index in seq(1, length( current_vector))){
        perc_info = pca$c1[[ cs_index]][decreasing_order][gene_index]^2
        perc_info_sum = perc_info_sum + perc_info
        cat("\n", row.names( pca$c1)[decreasing_order][gene_index], "\t:\t", signif((perc_info)*100, 3), "%\t:\t", signif( perc_info_sum*100 , 3), " %")
        if(! row.names( pca$c1)[decreasing_order][gene_index] %in% genes_maxCumulativeInfo_vector){
          genes_maxCumulativeInfo_vector = c(genes_maxCumulativeInfo_vector, row.names( pca$c1)[decreasing_order][gene_index])
        }
        if( perc_info_sum >= threshold_norm){
          break
        }
      }
    }
  }
  #do.call("grid.arrange", c(ggplot_list, nrow=3))
  
  variable_genes_withThreshold_v = unique(variable_genes_withThreshold_v)
  #return(ggplot_list)
  return(variable_genes_withThreshold_v)
}



# #####################################################
# #####################################################
# #####################################################
# finding outlier genes in a gene loadings landscape
# #####################################################
# #####################################################
# #####################################################
genes_outliers_fromLoadingsPCA_function <- function(pca_info, Dim_v = 1:3, bool_plot = TRUE, title=""){
  
  # Calcul the center of gravity of the gene loadings on the PCs
  # Modify the genes loading to center them at the orgin
  # Then calcul a dist matrix between the genes
  # Then make a dist threshold (median) and for each genes count the number of neighbour which 
  # are at a distance smaller than this threshold
  # Watch the distribution of the neighbours for all the genes 
  # and take the genes which have less neighbour than "median-1.96*sd"
  centered_genes_loadings_df = pca_info$co[,Dim_v]
  colMeans_genes_loadings = colMeans(centered_genes_loadings_df)
  centered_genes_loadings_df = t(t(centered_genes_loadings_df) - colMeans_genes_loadings)
  
  gene_loading_dist = dist(centered_genes_loadings_df)
  median_dist = median(gene_loading_dist)
  mad_dist = mad(gene_loading_dist)
  cutoff_dist = median_dist #+ 2*mad_dist
  nbr_of_neighbours_v = apply(gene_loading_dist,2,function(x) sum(x<(cutoff_dist)))
  #hist(nbr_of_neighbours_v,breaks=40)
  #abline(v=median(nbr_of_neighbours_v)-2*mad(nbr_of_neighbours_v), col="red")
  
  CI_low = median( nbr_of_neighbours_v) - 1.96 * sd( nbr_of_neighbours_v)
  nbr_of_neighbours_v_bool <- nbr_of_neighbours_v < CI_low
  #abline(v=CI_low, col="green")
  # plutot que IC faire un density plot et faire le cutoff au min entre les 2 max
  
  centered_genes_loadings_df$Outliers = nbr_of_neighbours_v_bool
  Outlier_genes_v = names(nbr_of_neighbours_v_bool)[nbr_of_neighbours_v_bool]
  
  if(bool_plot){
    hist(nbr_of_neighbours_v,breaks=40)
    abline(v=CI_low, col="green")
    for (i in 1:(length(Dim_v)-1)){
      for(j in (i+1):length(Dim_v)){
        Comp_x = paste0("Comp", Dim_v[i])
        Comp_y = paste0("Comp", Dim_v[j])
        gg_plot = ggplot(centered_genes_loadings_df, aes_string(Comp_x,Comp_y, color="Outliers"))+
          geom_point()+ ggtitle(title)
        print(gg_plot)
      }
    }
  }
  
  return (Outlier_genes_v)
}

# #####################################################
# #####################################################
# #####################################################
# finding outlier genes in a gene loadings landscape
# based on genes_outliers_fromLoadingsPCA_function
# #####################################################
# #####################################################
# #####################################################
getOutliersFromMultiDim_function <- function(df, bool_hist = FALSE, bool_plot = TRUE, title="Defining the outliers via NiakiClust"){
  
  # Calcul the center of gravity of the gene loadings on the PCs
  # Modify the genes loading to center them at the orgin
  # Then calcul a dist matrix between the genes
  # Then make a dist threshold (median) and for each genes count the number of neighbour which 
  # are at a distance smaller than this threshold
  # Watch the distribution of the neighbours for all the genes 
  # and take the genes which have less neighbour than "median-1.96*sd"
  centered_genes_loadings_df = df
  #centered_genes_loadings_df=sweep(df, ncol(df), colMeans(df), FUN="-")
  Dim_v = colnames(df)
  colMeans_genes_loadings = colMeans(centered_genes_loadings_df)
  centered_genes_loadings_df =as.data.frame( t(t(centered_genes_loadings_df) - colMeans_genes_loadings) )
  
  gene_loading_dist = as.matrix(dist(centered_genes_loadings_df))
  median_dist = median(gene_loading_dist)
  #mad_dist = mad(gene_loading_dist)
  cutoff_dist = median_dist #+ 2*mad_dist
  nbr_of_neighbours_v = apply(gene_loading_dist,2,function(x) sum(x<(cutoff_dist)))
  #hist(nbr_of_neighbours_v,breaks=40)
  #abline(v=median(nbr_of_neighbours_v)-2*mad(nbr_of_neighbours_v), col="red")
  
  CI_low = median( nbr_of_neighbours_v) - 1.96 * sd( nbr_of_neighbours_v)
  nbr_of_neighbours_v_bool <- nbr_of_neighbours_v < CI_low
  #abline(v=CI_low, col="green")
  # plutot que IC faire un density plot et faire le cutoff au min entre les 2 max
  
  centered_genes_loadings_df$Outliers = sapply(nbr_of_neighbours_v_bool,function(x) 
    if(x){"outlier"}else{"grouped"})
  Outlier_genes_v = names(nbr_of_neighbours_v_bool)[nbr_of_neighbours_v_bool]
  
  result_v = sapply(rownames(centered_genes_loadings_df), function(x) 
    if(x %in% Outlier_genes_v){"outlier"}else{"grouped"})
  
  if(bool_hist){
    hist(nbr_of_neighbours_v,breaks=40)
    abline(v=CI_low, col="green")
  }
  if(bool_plot){
    for (i in 1:(length(Dim_v)-1)){
      for(j in (i+1):length(Dim_v)){
        Comp_x = Dim_v[i]
        Comp_y = Dim_v[j]
        gg_plot = ggplot(centered_genes_loadings_df, aes_string(Comp_x,Comp_y, color="Outliers"))+
          geom_point()+ ggtitle(title)
        print(gg_plot)
      }
    }
  }
  
  return (result_v)
}

# #####################################################
# #####################################################
# #####################################################
# radial coordinate of genes or cells loadings from PCA
# #####################################################
# #####################################################
# #####################################################
radialCoord_loadings_function <- function(x_v, y_v, Id_v, x_mean =mean(x_v), y_mean=mean(y_v) ){
  # x_v y_v and Id_v need to have the same order !
  # x_mean=mean(x_v)
  # y_mean=mean(y_v)
  x_v = x_v - x_mean
  y_v = y_v - y_mean
  radian_v=c()
  radius_v=c()
  for(i in 1:length(x_v)){
    x = x_v[i]
    y = y_v[i]
    radius=sqrt(x^2+y^2)
    if(x<0 && y<0){
      radian=atan(y/x)
    } else if(x>0 && y<0){
      radian=pi+atan(y/x)
    } else if(x>0 && y>0){
      radian=pi+atan(y/x)
    } else if(x<0 && y>0){
      radian=2*pi+atan(y/x)
    }
    radian_v=c(radian_v,radian)
    radius_v=c(radius_v,radius)
  }
  
  df= data.frame(ID= Id_v,
                 Radian = radian_v,
                 Radius = radius_v, stringsAsFactors = FALSE)
  return(df)
}

# #####################################################
# #####################################################
# #####################################################
# Function that execute a t-SNE (bh-SNE)
# #####################################################
# #####################################################
# #####################################################
tsne_function <- function( raw_data_for_tsne, VectorToColor = NULL, NameOfVectorToColor = "" ,perplexity_value=30, bool_facet_wrap=FALSE, VectorFacetWrap=NA,title="", bool_check_duplicates=TRUE){
  limits_allFigs=c(0,27)
  
  if(length(VectorFacetWrap)<=1){
    VectorFacetWrap = VectorToColor
  }
  
  # special case if we want to color the pca by gene expression (continuous variable)
  if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
    VectorToColor[which(VectorToColor==0)]=NA
  }
  
  #raw_data_for_tsne=raw_data_all_df[,GENES_TO_ANALYSE]
  if(min(rowSums(raw_data_for_tsne))==0){
    VectorToColor = VectorToColor[-which(rowSums(raw_data_for_tsne)==0)]
    raw_data_for_tsne = raw_data_for_tsne[-which(rowSums(raw_data_for_tsne)==0),]
  }
  
  
  # bh-SNE: t-SNE using Barnes-Hut implementation
  # #############################################
  set.seed(42) # Sets seed for reproducibility
  tsne_out <- Rtsne(as.matrix(raw_data_for_tsne),initial_dims = 20, perplexity=perplexity_value,check_duplicates=bool_check_duplicates) #initial_dims = 50
  
  # take the tSNE Coordinates
  # #########################
  raw_data_for_tsne$tsne_x = tsne_out$Y[,1]
  raw_data_for_tsne$tsne_y = tsne_out$Y[,2]
  raw_data_for_tsne$DataToColor = VectorToColor
  
  if(bool_facet_wrap==FALSE){
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
        
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    }
  }
  else{
    raw_data_for_tsne$DataToFacetWrap = VectorFacetWrap
    raw_data_for_tsne$DataToFacetWrap2 =raw_data_for_tsne$DataToFacetWrap
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    }
    
  }
}

# #####################################################
# #####################################################
# #####################################################
# Function that run t-SNE (bh-SNE)
# #####################################################
# #####################################################
# #####################################################
run_tsne_function <- function( raw_data_for_tsne ,perplexity_value=30, bool_check_duplicates=TRUE){
  
  #raw_data_for_tsne=raw_data_all_df[,GENES_TO_ANALYSE]
  if(min(rowSums(raw_data_for_tsne))==0){
    raw_data_for_tsne = raw_data_for_tsne[-which(rowSums(raw_data_for_tsne)==0),]
  }
  
  
  # bh-SNE: t-SNE using Barnes-Hut implementation
  # #############################################
  set.seed(42) # Sets seed for reproducibility
  tsne_out <- Rtsne(as.matrix(raw_data_for_tsne),initial_dims = 20, perplexity=perplexity_value,check_duplicates=bool_check_duplicates) #initial_dims = 50
  
  # take the tSNE Coordinates
  # #########################
  raw_data_for_tsne$tsne_x = tsne_out$Y[,1]
  raw_data_for_tsne$tsne_y = tsne_out$Y[,2]
  return(raw_data_for_tsne[,c("tsne_x","tsne_y")])
}
# #####################################################
# #####################################################
# #####################################################
# Function that execute a t-SNE (bh-SNE)
# #####################################################
# #####################################################
# #####################################################
plot_tsne_function <- function( raw_data_for_tsne, tsne_coord_df,VectorToColor = NULL, NameOfVectorToColor = "" , VectorFacetWrap=NA,title="", bool_plot=TRUE){
  limits_allFigs=c(0,27)
  
  #if(length(VectorFacetWrap)<=1){
  #  VectorFacetWrap = VectorToColor
  #}
  
  # special case if we want to color the pca by gene expression (continuous variable)
  if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
    VectorToColor[which(VectorToColor==0)]=NA
  }
  
  #raw_data_for_tsne=raw_data_all_df[,GENES_TO_ANALYSE]
  if(min(rowSums(raw_data_for_tsne))==0){
    VectorToColor = VectorToColor[-which(rowSums(raw_data_for_tsne)==0)]
    raw_data_for_tsne = raw_data_for_tsne[-which(rowSums(raw_data_for_tsne)==0),]
  }
  
  # take the tSNE Coordinates
  # #########################
  raw_data_for_tsne$tsne_x = tsne_coord_df$tsne_x
  raw_data_for_tsne$tsne_y = tsne_coord_df$tsne_y
  raw_data_for_tsne$DataToColor = VectorToColor
  
  if(length(VectorFacetWrap)<=1){
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        #theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
      
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        #theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    }
  }
  else{
    raw_data_for_tsne$DataToFacetWrap = VectorFacetWrap
    raw_data_for_tsne$DataToFacetWrap2 =raw_data_for_tsne$DataToFacetWrap
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        #theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        #theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    }
  }
  
  return(tSNE_plot)
}


# #####################################################
# #####################################################
# #####################################################
# volcano plot - combine LRT expression genes (modified)
# #####################################################
# #####################################################
# #####################################################
volcano_plot_LRT_function <- function(raw_data_for_volcano_plot, VectorToSplit, FigureLetterName=NULL){
  
  variable_to_return=list()
  
  cell_attribute_toSplit = "ParameterToSplit"
  genes_vector = colnames(raw_data_for_volcano_plot)
  
  # Build the dataframe with 1 gene per line to build the FluidigmAssay object from data
  # -- Build a basic data frame with the correct number of lines
  one_gene_per_line_df = data.frame( count=seq( 1,ncol(raw_data_for_volcano_plot)*nrow( raw_data_for_volcano_plot)), stringsAsFactors=FALSE)
  
  # -- Add the vector_to_split and UniqueCellID
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( VectorToSplit, each= ncol(raw_data_for_volcano_plot)))
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( rownames(raw_data_for_volcano_plot), each= ncol(raw_data_for_volcano_plot)))
  
  # -- Add the duplication of the genes names
  one_gene_per_line_df$Gene = rep( genes_vector, nrow( raw_data_for_volcano_plot))
  
  # -- Add the value of expressions
  gene_et = vector()
  for( row_index in seq( 1, nrow( raw_data_for_volcano_plot))){
    gene_et = append( gene_et, unlist( raw_data_for_volcano_plot[ row_index, ]))
  }
  one_gene_per_line_df$Et = (gene_et) #exp(gene_et) for normalized scRNAseq data
  
  # -- Set the name of the dataframe
  names( one_gene_per_line_df) = c("count",cell_attribute_toSplit,"UniqueCellID","Gene","Et")
  
  # Build the FluidigmAssay object from MAST library
  # using the header of dataframe
  vbeta.fa = FluidigmAssay(one_gene_per_line_df, idvars="UniqueCellID", 
                           primerid="Gene",
                           measurement="Et",
                           #ncells= "CellNumber",
                           geneid="Gene", 
                           cellvars=c("UniqueCellID"),
                           phenovars=c(cell_attribute_toSplit),
                           id= raw_data_filename)
  
  
  cell_attribute_toSplit_vector = sort(unique(VectorToSplit))
  
  for( index1 in 1:(length( cell_attribute_toSplit_vector)-1)){
    current_ref_state =  cell_attribute_toSplit_vector[ index1]
    
    lrt_test <- LRT(vbeta.fa, cell_attribute_toSplit, referent = current_ref_state)
    
    for( index2 in (index1+1):length( cell_attribute_toSplit_vector)){
      current_alt_state =  cell_attribute_toSplit_vector[ index2]
      # Get the separated two groups of cells
      vbeta.1cell.nozero.nobad.ref = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_ref_state)
      vbeta.1cell.nozero.nobad.alt = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_alt_state)
      
      
      # Compute the difference of mean expression on the two groups and build the mean and pvalue set for the volcano plot
      mean_diff_set = vector()
      #FAIRE DIFFERENCE ENTRE PERCENTAGE DEXPRESSION
      lrt_pval_set = vector()
      for( gene_name in genes_vector){
        mean_ref = mean( exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name])
        mean_alt = mean( exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name])
        pval = lrt_test[ which( lrt_test[[cell_attribute_toSplit]] == current_alt_state & lrt_test$primerid == gene_name),"p.value"]
        mean_diff_set = append( mean_diff_set, mean_ref - mean_alt)
        lrt_pval_set = append( lrt_pval_set, pval)
      }
      lrt_pval_set_adjustMultiComparison = p.adjust(lrt_pval_set, method="BH")
      lrt_pval_set_adjustMultiComparison = -log10(lrt_pval_set_adjustMultiComparison)
      lrt_pval_set = -log10(lrt_pval_set)
      mean_thresholds = quantile(mean_diff_set, c(.05, .95))
      mean_thresholds_2 = c( mean( mean_diff_set) - sd( mean_diff_set), mean( mean_diff_set) + sd( mean_diff_set))
      
      
      # Compute the list of genes that are the most extreme in volcano plot (upper left and upper right)    
      mean_index_list = which( mean_diff_set <= mean_thresholds[1] | mean_diff_set >= mean_thresholds[2])
      mean_index_list_2 = which( mean_diff_set <= mean_thresholds_2[1] | mean_diff_set >= mean_thresholds_2[2])
      pval_index_list = which( lrt_pval_set_adjustMultiComparison >= -log10(0.05))
      final_index_list = intersect( mean_index_list, pval_index_list)
      final_index_list_2 = intersect( mean_index_list_2, pval_index_list)
      
      
      # -- THE BRAND NEW & TOTALY AWESOME VOLCANO PLOT --
      # #################################################
      volcano_plot_data_frame=data.frame(row.names=genes_vector, MeanDiffSet=mean_diff_set, LRTpVal=lrt_pval_set_adjustMultiComparison) #lrt_pval_set
      
      volcano_plot_data_frame$Significant <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list], "< 0.05 Perc or > 0.95 Perc", "Not Sig")
      volcano_plot_data_frame$Significant2 <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_2], "abs( Z-score) >= 1", "Not Sig")
      
      if(TRUE){ #!is.null(FigureLetterName)
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = LRTpVal)) +
          #ggtitle(paste( "Volcano plot of expression mean difference vs LRT pvalue\n",current_ref_state, "vs", current_alt_state)) + 
          xlab(paste( "Diff of exprs means between", current_ref_state, "and", current_alt_state)) +
          ylab("-log10 p.value")+ 
          ggtitle(paste(current_alt_state,"(overexprs left) vs",current_ref_state,"(overexprs right)"))+
          geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds[1], colour="firebrick2", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds[2], colour="firebrick2", alpha=0.4) +
          geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
          geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
          #geom_point(aes(color = Significant, size =Significant)) +
          #scale_color_manual(values = c("firebrick2", "darkgrey")) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(aes(color = Significant2),size=3) +
          scale_color_manual(values = c("firebrick2", "gray70"), name="Significant genes") +
          #theme_bw(base_size = 12) + #theme(legend.position = "bottom") +
          theme( plot.title = element_text( face="bold", size=12),
                 aspect.ratio=1,
                 legend.position="none") +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
        scale_x_continuous(limits=c(-max(abs(volcano_plot_data_frame$MeanDiffSet)),max(abs(volcano_plot_data_frame$MeanDiffSet)) )) +
          scale_y_continuous(expand=c(0,0),limits=c(0,max(abs(volcano_plot_data_frame$LRTpVal))+0.5))
        print(volcano_plot_ggplot)
        #assign(paste0("Physio_",FigureLetterName,LETTERS[index1]), volcano_plot_ggplot)
        #save(list=paste0("Physio_",FigureLetterName,LETTERS[index1]), file=paste0(RDATA_DIR,"Physio-",FigureLetterName,LETTERS[index1],".RData"))  
        #ggsave(plot=volcano_plot_ggplot,height=6,width=6,dpi=200, filename=paste0( FIGURE_DIR, "Physio-",FigureLetterName,LETTERS[index1],".pdf"), useDingbats=FALSE)
      }
      
      if(FALSE){
        cat("\nMost significant genes in both mean difference and LRT p-value (", current_ref_state, "vs", current_alt_state, "):")
        cat( "\n\n lesser than 5 percentile or greater than 95 percentile:\n")
        print( genes_vector[final_index_list], quote=FALSE, row.names=FALSE)
        cat( "\n\n abs( Z-score) >= 1 :\n")
        print( genes_vector[final_index_list_2], quote=FALSE, row.names=FALSE)
      }
      
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]]=list()
      # put the genes which are more expressed in current_ref_state in one vector
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_ref_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]>0)]]
      # put the genes which are more expressed in current_alt_state in one vector  
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_alt_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]<0)]]
      
      
    }
  }
  
  for(i in 1:length(variable_to_return)){
    cat("\n\nMost significant genes in both mean difference and LRT p-value (", names(variable_to_return)[i], "):\n")
    for(j in 1:length(variable_to_return[[i]])){
      cat("Genes overexpressed", names(variable_to_return[[i]])[j], ":\n")
      print(variable_to_return[[i]][[j]])
    }
  }
  
  all_categories_v = sort(unique(VectorToSplit))
  if(length(all_categories_v)>2){
    cat("\n\n*****************************************************************************")
    cat("\n*****************************************************************************")
    cat("\n*****************************************************************************")
    cat("\n\nList of the genes which are overexpressed or underexpressed by CATEGORY vs ALL OTHER CATEGORIES:\n\n")
    variable_to_return
    diff_expr_list=list()
    
    for(i in 1:length(all_categories_v)){
      up_exprs_genes_v = c()
      down_exprs_genes_v = c()
      for(j in 1:length(variable_to_return)){
        if(all_categories_v[i]==names(variable_to_return[[j]])[1] || all_categories_v[i]==names(variable_to_return[[j]])[2]){
          if(length(up_exprs_genes_v)==0){
            up_exprs_genes_v = variable_to_return[[j]][[all_categories_v[i] ]]
          } else{
            up_exprs_genes_v = intersect(up_exprs_genes_v,variable_to_return[[j]][[all_categories_v[i] ]])
          }
          if(length(down_exprs_genes_v)==0){
            down_exprs_genes_v = variable_to_return[[j]][[which(names(variable_to_return[[j]])!=all_categories_v[i]) ]]
          } else{
            down_exprs_genes_v = intersect(down_exprs_genes_v,variable_to_return[[j]][[which(names(variable_to_return[[j]])!=all_categories_v[i]) ]])
          }
        }
      }
      
      diff_expr_list[[all_categories_v[i] ]] = list(UP=up_exprs_genes_v,DOWN=down_exprs_genes_v)
      
    }
    #print(diff_expr_list)
    for(i in 1:length(diff_expr_list)){
      cat("\n\n",names(diff_expr_list)[i],":")
      cat("\n**************************")
      for(j in 1:length(diff_expr_list[[i]])){
        cat("\n",names(diff_expr_list[[i]])[j],":\n")
        print(diff_expr_list[[i]][[j]])
      }
    }
  }
  
  return(variable_to_return)
}

# ######################################################
# ######################################################
# ######################################################
# volcano plot - combine LRT expression genes (modified)
# ######################################################
# ######################################################
# ######################################################
volcano_plot_LRT_function_v2 <- function(raw_data_for_volcano_plot, VectorToSplit, FigureLetterName=NULL){
  
  variable_to_return=list()
  
  cell_attribute_toSplit = "ParameterToSplit"
  genes_vector = colnames(raw_data_for_volcano_plot)
  
  # Build the dataframe with 1 gene per line to build the FluidigmAssay object from data
  # -- Build a basic data frame with the correct number of lines
  one_gene_per_line_df = data.frame( count=seq( 1,ncol(raw_data_for_volcano_plot)*nrow( raw_data_for_volcano_plot)), stringsAsFactors=FALSE)
  
  # -- Add the vector_to_split and UniqueCellID
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( VectorToSplit, each= ncol(raw_data_for_volcano_plot)))
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( rownames(raw_data_for_volcano_plot), each= ncol(raw_data_for_volcano_plot)))
  
  # -- Add the duplication of the genes names
  one_gene_per_line_df$Gene = rep( genes_vector, nrow( raw_data_for_volcano_plot))
  
  # -- Add the value of expressions
  gene_et = vector()
  for( row_index in seq( 1, nrow( raw_data_for_volcano_plot))){
    gene_et = append( gene_et, unlist( raw_data_for_volcano_plot[ row_index, ]))
  }
  one_gene_per_line_df$Et = gene_et
  
  # -- Set the name of the dataframe
  names( one_gene_per_line_df) = c("count",cell_attribute_toSplit,"UniqueCellID","Gene","Et")
  
  # Build the FluidigmAssay object from MAST library
  # using the header of dataframe
  vbeta.fa = FluidigmAssay(one_gene_per_line_df, idvars="UniqueCellID", 
                           primerid="Gene",
                           measurement="Et",
                           #ncells= "CellNumber",
                           geneid="Gene", 
                           cellvars=c("UniqueCellID"),
                           phenovars=c(cell_attribute_toSplit),
                           id= raw_data_filename)
  
  
  cell_attribute_toSplit_vector = sort(unique(VectorToSplit))
  
  for( index1 in 1:(length( cell_attribute_toSplit_vector)-1)){
    current_ref_state =  cell_attribute_toSplit_vector[ index1]
    
    lrt_test <- LRT(vbeta.fa, cell_attribute_toSplit, referent = current_ref_state)
    
    for( index2 in (index1+1):length( cell_attribute_toSplit_vector)){
      current_alt_state =  cell_attribute_toSplit_vector[ index2]
      # Get the separated two groups of cells
      vbeta.1cell.nozero.nobad.ref = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_ref_state)
      vbeta.1cell.nozero.nobad.alt = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_alt_state)
      
      
      # Compute the difference of mean expression on the two groups and build the mean and pvalue set for the volcano plot
      mean_diff_set_noZero = vector()
      mean_diff_set = vector()
      perc_diff_set = vector()
      #FAIRE DIFFERENCE ENTRE PERCENTAGE DEXPRESSION
      lrt_pval_set = vector()
      for( gene_name in genes_vector){
        mean_ref_noZero = mean( exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name][exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name]!=0])
        mean_alt_noZero = mean( exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name][exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name]!=0])
        if(is.nan(mean_ref)){mean_ref=0}
        if(is.nan(mean_alt)){mean_alt=0}
        perc_ref = sum(exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name]==0)/length(exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name])
        perc_alt = sum(exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name]==0)/length(exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name])
        mean_ref = mean( exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name])
        mean_alt = mean( exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name])
        pval = lrt_test[ which( lrt_test[[cell_attribute_toSplit]] == current_alt_state & lrt_test$primerid == gene_name),"p.value"]
        mean_diff_set_noZero = append( mean_diff_set_noZero, mean_ref_noZero - mean_alt_noZero)
        mean_diff_set = append( mean_diff_set, mean_ref - mean_alt)
        perc_diff_set = append( perc_diff_set, perc_ref - perc_alt)
        lrt_pval_set = append( lrt_pval_set, -log10(pval))
      }
      mean_thresholds_noZero = quantile(mean_diff_set_noZero, c(.05, .95))
      mean_thresholds = quantile(mean_diff_set, c(.05, .95))
      mean_thresholds_2 = c( mean( mean_diff_set) - sd( mean_diff_set), mean( mean_diff_set) + sd( mean_diff_set))
      
      
      # Compute the list of genes that are the most extreme in volcano plot (upper left and upper right)    
      mean_index_list = which( mean_diff_set <= mean_thresholds[1] | mean_diff_set >= mean_thresholds[2])
      mean_index_list_2 = which( mean_diff_set <= mean_thresholds_2[1] | mean_diff_set >= mean_thresholds_2[2])
      pval_index_list = which( lrt_pval_set >= -log10(0.05))
      final_index_list = intersect( mean_index_list, pval_index_list)
      final_index_list_2 = intersect( mean_index_list_2, pval_index_list)
      
      final_index_list_basedOnlyOnPval = pval_index_list
      
      # -- THE BRAND NEW & TOTALY AWESOME VOLCANO PLOT --
      # #################################################
      volcano_plot_data_frame=data.frame(row.names=genes_vector, MeanDiffSet=mean_diff_set_noZero, PercDiffSet=perc_diff_set, LRTpVal=lrt_pval_set)
      
      volcano_plot_data_frame$Significant <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list], "< 0.05 Perc or > 0.95 Perc", "Not Sig")
      volcano_plot_data_frame$Significant2 <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_2], "abs( Z-score) >= 1", "Not Sig")
      
      volcano_plot_data_frame$SignificantOnlyPval <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_basedOnlyOnPval], "pval>0.05","Not Sig")
      
      if(TRUE){ #!is.null(FigureLetterName)
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = PercDiffSet)) +
          #ggtitle(paste( "Volcano plot of expression mean difference vs LRT pvalue\n",current_ref_state, "vs", current_alt_state)) + 
          xlab(paste( "Diff of expr means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          #geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
          #scale_color_manual(values = c("firebrick2", "darkgrey")) +
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(aes(color = Significant2),size=3) +
          scale_color_manual(values = c("firebrick2", "gray70"), name="Significant genes") +
          #theme_bw(base_size = 12) + #theme(legend.position = "bottom") +
          theme( plot.title = element_text( face="bold", size=12),
                 aspect.ratio=1,
                 legend.position="none") +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) #+
          #scale_x_continuous(limits=c(-max(abs(volcano_plot_data_frame$MeanDiffSet)),max(abs(volcano_plot_data_frame$MeanDiffSet)) )) +
          #scale_y_continuous(expand=c(0,0),limits=c(0,max(abs(volcano_plot_data_frame$LRTpVal))+0.5))
        print(volcano_plot_ggplot)
        #assign(paste0("Physio_",FigureLetterName,LETTERS[index1]), volcano_plot_ggplot)
        #save(list=paste0("Physio_",FigureLetterName,LETTERS[index1]), file=paste0(RDATA_DIR,"Physio-",FigureLetterName,LETTERS[index1],".RData"))  
        #ggsave(plot=volcano_plot_ggplot,height=6,width=6,dpi=200, filename=paste0( FIGURE_DIR, "Physio-",FigureLetterName,LETTERS[index1],".pdf"), useDingbats=FALSE)
        
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = PercDiffSet)) +
          #ggtitle(paste( "Volcano plot of expression mean difference vs LRT pvalue\n",current_ref_state, "vs", current_alt_state)) + 
          xlab(paste( "Difference of expression means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          #geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
          #scale_color_manual(values = c("firebrick2", "darkgrey")) +
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(aes(color = SignificantOnlyPval),size=3) +
          scale_color_manual(values = c("gray70", "firebrick2"), name="Significant genes") +
          #theme_bw(base_size = 12) + #theme(legend.position = "bottom") +
          theme( plot.title = element_text( face="bold", size=12),
                 aspect.ratio=1,
                 legend.position="none") +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], SignificantOnlyPval=="pval>0.05"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], SignificantOnlyPval=="pval>0.05"))),
            size = 3, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], SignificantOnlyPval=="pval>0.05"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], SignificantOnlyPval=="pval>0.05"))),
            size = 3, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) #+
        #scale_x_continuous(limits=c(-max(abs(volcano_plot_data_frame$MeanDiffSet)),max(abs(volcano_plot_data_frame$MeanDiffSet)) )) +
        #scale_y_continuous(expand=c(0,0),limits=c(0,max(abs(volcano_plot_data_frame$LRTpVal))+0.5))
        print(volcano_plot_ggplot)
        
        
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = PercDiffSet, color=LRTpVal)) +
          #ggtitle(paste( "Volcano plot of expression mean difference vs LRT pvalue\n",current_ref_state, "vs", current_alt_state)) + 
          xlab(paste( "Difference of expression means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          #geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
          #scale_color_manual(values = c("firebrick2", "darkgrey")) +
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(size=3) +
          scale_colour_gradientn(colours = rainbow(7))
        print(volcano_plot_ggplot)
        
        
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = PercDiffSet, y = LRTpVal, color=LRTpVal)) +
          #ggtitle(paste( "Volcano plot of expression mean difference vs LRT pvalue\n",current_ref_state, "vs", current_alt_state)) + 
          xlab(paste( "Difference of expression means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          #geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
          #geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
          #scale_color_manual(values = c("firebrick2", "darkgrey")) +
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(size=3) +
          scale_colour_gradientn(colours = rainbow(7))
        print(volcano_plot_ggplot)
      }
      
      if(FALSE){
        cat("\nMost significant genes in both mean difference and LRT p-value (", current_ref_state, "vs", current_alt_state, "):")
        cat( "\n\n lesser than 5 percentile or greater than 95 percentile:\n")
        print( genes_vector[final_index_list], quote=FALSE, row.names=FALSE)
        cat( "\n\n abs( Z-score) >= 1 :\n")
        print( genes_vector[final_index_list_2], quote=FALSE, row.names=FALSE)
      }
      
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]]=list()
      # put the genes which are more expressed in current_ref_state in one vector
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_ref_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]>0)]]
      # put the genes which are more expressed in current_alt_state in one vector  
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_alt_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]<0)]]
      
      
    }
  }
  
  for(i in 1:length(variable_to_return)){
    cat("\n\nMost significant genes in both mean difference and LRT p-value (", names(variable_to_return)[i], "):\n")
    for(j in 1:length(variable_to_return[[i]])){
      cat("Genes overexpressed", names(variable_to_return[[i]])[j], ":\n")
      print(variable_to_return[[i]][[j]])
    }
  }
  
  return(variable_to_return)
}

# #####################################################
# #####################################################
# #####################################################
# Hierarchical clustering 
# #####################################################
# #####################################################
# #####################################################
hierarchical_clustering_function <- function (raw_data_for_biclustering, genes_for_clustering=ALL_FILTERED_GENES, side_bar_colored_variable="Plate", clust_method = "spearman"){
  
  colors_to_use=gg_color_hue(length(unique(raw_data_for_biclustering[,side_bar_colored_variable])))
  
  title_pheatmap=paste("Hierarchical clustering")
  
  # -- data for the study
  raw_data_transposed=t(raw_data_for_biclustering[, genes_for_clustering])
  # - colnames based on Pheno
  colnames(raw_data_transposed)=raw_data_for_biclustering$UniqueCellID[1:ncol(raw_data_transposed)]
  
  
  # we remove the genes with no variation, otherwise we can not do the hierarchical clustering on genes
  genes_no_variation=c()
  for (i in 1:nrow(raw_data_transposed)){
    sum_row=sum(raw_data_transposed[i,])
    if(sum_row==0){
      genes_no_variation=c(genes_no_variation,rownames(raw_data_transposed)[i])
    }
  }
  raw_data_transposed=raw_data_transposed[which(!rownames(raw_data_transposed) %in% genes_no_variation),]
  
  print(paste("We remove the genes with no variation: ",paste(genes_no_variation,collapse=" "),sep=""))
  
  cells_no_variation=c()
  index_cells_no_variation=c()
  for (i in 1:ncol(raw_data_transposed)){
    sum_col=sum(raw_data_transposed[,i])
    if(sum_col==0){
      cells_no_variation=c(cells_no_variation,colnames(raw_data_transposed)[i])
      index_cells_no_variation=c(index_cells_no_variation,i)
    }
  }
  raw_data_transposed=raw_data_transposed[,which(!colnames(raw_data_transposed) %in% cells_no_variation)]
  
  print(paste("We remove the cells with no variation: ",paste(cells_no_variation,collapse=" "),sep=""))
  
  
  
  # #######################
  # hierarchical clustering
  # #######################
  
  # the new heatmap.2 HIERARCHICAL CLUSTERING
  # #############################################
  sidebar_color_vector=raw_data_for_biclustering[,side_bar_colored_variable]
  if(length(index_cells_no_variation)>0){
    sidebar_color_vector=raw_data_for_biclustering[-index_cells_no_variation,side_bar_colored_variable]
  } 
  unique_elements=unique(raw_data_for_biclustering[,side_bar_colored_variable])
  #vector_colors=colors_to_use #c("black",gg_color_hue(length(unique_elements)-1))
  for(i in 1:length(sidebar_color_vector)){
    sidebar_color_vector[i]= colors_to_use[which(unique_elements %in% c(sidebar_color_vector[i]))]
  }
  
  
  if(clust_method=="euclidean"){
    heatmap.2(raw_data_transposed, 
              trace="none", 
              col=c("grey80",colorRampPalette(c("white", "blue"))(n = 400)),
              #Rowv=as.dendrogram(hr), 
              #Colv=as.dendrogram(hc),
              labCol = FALSE,
              margins = c(2, 20),
              ColSideColors=sidebar_color_vector) 
    par(lend = 1)           # square line ends for the color legend
    legend("topright",      # location of the legend on the heatmap plot
           legend = unique_elements, # category labels
           col =colors_to_use,  # color key
           lty= 1,             # line style
           lwd = 10            # line width
    )
  } else {
    hr <- hclust(as.dist(1-cor(t(raw_data_transposed), method=clust_method)), method="complete")
    hc <- hclust(as.dist(1-cor(raw_data_transposed, method=clust_method)), method="complete")
    heatmap.2(raw_data_transposed, 
              trace="none", 
              col=c("grey80",colorRampPalette(c("white", "blue"))(n = 400)),
              Rowv=as.dendrogram(hr), 
              Colv=as.dendrogram(hc),
              labCol = FALSE,
              margins = c(2, 20),
              ColSideColors=sidebar_color_vector) 
    par(lend = 1)           # square line ends for the color legend
    legend("topright",      # location of the legend on the heatmap plot
           legend = unique_elements, # category labels
           col =colors_to_use,  # color key
           lty= 1,             # line style
           lwd = 10            # line width
    )
  }
  
}

# #####################################################
# #####################################################
# #####################################################
# Biclustering backSPIN
# #####################################################
# #####################################################
# #####################################################
biclustering_function <- function (raw_data_for_biclustering, side_bar_colored_variable="Plate", genes_for_clustering=ALL_FILTERED_GENES, fig_title=""){
  
  colors_to_use=gg_color_hue(length(unique(raw_data_for_biclustering[,side_bar_colored_variable])))
  
  
  title_ggplot=paste("Biclustering:",fig_title)
  
  # -- data for the study
  raw_data_transposed=t(raw_data_for_biclustering[, genes_for_clustering])
  # - colnames based on Pheno
  colnames(raw_data_transposed)=raw_data_for_biclustering$UniqueCellID[1:ncol(raw_data_transposed)]
  
  
  # we remove the genes with no variation, otherwise we can not do the hierarchical clustering on genes
  genes_no_variation=c()
  for (i in 1:nrow(raw_data_transposed)){
    sum_row=sum(raw_data_transposed[i,])
    if(sum_row==0){
      genes_no_variation=c(genes_no_variation,rownames(raw_data_transposed)[i])
    }
  }
  raw_data_transposed=raw_data_transposed[which(!rownames(raw_data_transposed) %in% genes_no_variation),]
  
  print(paste("We remove the genes with no variation: ",paste(genes_no_variation,collapse=" "),sep=""))
  
  cells_no_variation=c()
  index_cells_no_variation=c()
  for (i in 1:ncol(raw_data_transposed)){
    sum_col=sum(raw_data_transposed[,i])
    if(sum_col==0){
      cells_no_variation=c(cells_no_variation,colnames(raw_data_transposed)[i])
      index_cells_no_variation=c(index_cells_no_variation,i)
    }
  }
  raw_data_transposed=raw_data_transposed[,which(!colnames(raw_data_transposed) %in% cells_no_variation)]
  
  print(paste("We remove the cells with no variation: ",paste(cells_no_variation,collapse=" "),sep=""))
  
  # #####################
  # backSPIN biclustering
  # https://github.com/linnarsson-lab/BackSPIN
  # http://hms-dbmi.github.io/scw/analysis-of-heterogeneity-and-subpopulations.html
  # #####################
  
  # set the backSpin algorithm dir
  backSPIN_algorithm_python="~/workspace/Logiciels/BackSPIN-1.0/backSPIN.py"
  
  data_dir="~/workspace/"
  # Delete CEF files if they already exist
  file_backSPIN=paste(data_dir,"datafile_backSPIN.cef",sep="")
  if (file.exists(file_backSPIN)){
    file.remove(file_backSPIN)
  }
  file_backSPIN_clustered=paste(data_dir,"datafile_backSPIN.clustered.cef",sep="")
  if (file.exists(file_backSPIN_clustered)){
    file.remove(file_backSPIN_clustered)
  }
  
  # creation of the CEF file (the file used by backSPIN)
  # ####################################################
  write(file=file_backSPIN, paste("CEF\t0\t1\t1", nrow(raw_data_transposed), ncol(raw_data_transposed), "0", sep="\t"))
  write(file=file_backSPIN, paste("\tcell", paste(colnames(raw_data_transposed), collapse="\t"), sep='\t'), append=T)
  x <- cbind(rep("",nrow(raw_data_transposed)), raw_data_transposed)
  write(file=file_backSPIN, paste(c("gene"), sep="\t"), append=T)
  write.table(file=file_backSPIN, x, append=T, col.names=F, row.names=T, quote=F, sep="\t")
  
  # call the backSPIN algorithm in python which create a new CEF file where cells and genes are biclustered
  # #######################################################################################################
  system(paste(backSPIN_algorithm_python," -i ", file_backSPIN, " -o ", file_backSPIN_clustered," -v -d 3",sep=""))
  
  # retrieve the data (to create the heatmap) generated by the backSPIN algorithm 
  # #############################################################################
  backSPIN_clustered_data_frame <- read.delim(file_backSPIN_clustered,stringsAsFactors=F,sep="\t",skip=1,header=F)
  index_first_group_gene=match(0,backSPIN_clustered_data_frame[-c(1,2,3,4,5,6),3])
  index_second_group_gene=match(1,backSPIN_clustered_data_frame[-c(1,2,3,4,5,6),3])
  index_first_group_cell=match(0,backSPIN_clustered_data_frame[3,-c(1,2,3,4,5,6)])
  index_second_group_cell=match(1,backSPIN_clustered_data_frame[3,-c(1,2,3,4,5,6)])
  
  names_row=backSPIN_clustered_data_frame[-c(1,2,3,4,5,6),1] # lines specific to our type of data (check the raw CEF file in the Terminal to understand: "< datafile_backSPIN.clustered.cef cef view")
  names_col=backSPIN_clustered_data_frame[1,-c(1,2,3,4,5,6)]
  backSPIN_clustered_data_frame <- backSPIN_clustered_data_frame[-c(1:6),-c(1:6)]
  backSPIN_clustered_data_frame=matrix(as.numeric(unlist(backSPIN_clustered_data_frame)),nrow=nrow(backSPIN_clustered_data_frame))
  rownames(backSPIN_clustered_data_frame)=names_row
  colnames(backSPIN_clustered_data_frame)=names_col
  
  
  # #######################################
  # --- PLOT USING ggplot2: + geom_raster()
  # #######################################
  
  # reorganization of the data for ggplot
  
  # Build the dataframe with 1 gene per line to build the FluidigmAssay object from data
  # -- Build a basic data frame with the correct number of lines
  one_gene_per_line_df = data.frame( count=seq( 1,nrow(backSPIN_clustered_data_frame)*ncol(backSPIN_clustered_data_frame) ), stringsAsFactors=FALSE)
  
  # -- Add the duplication of the genes names
  #one_gene_per_line_df$Genes = rep( rownames(backSPIN_clustered_data_frame), ncol(backSPIN_clustered_data_frame))
  #one_gene_per_line_df$Genes = factor(rep( rownames(backSPIN_clustered_data_frame), ncol(backSPIN_clustered_data_frame)), levels=rev(unique(rep( rownames(backSPIN_clustered_data_frame), ncol(backSPIN_clustered_data_frame)))))
  one_gene_per_line_df$Genes = structure( rep(rev(1:nrow(backSPIN_clustered_data_frame)),ncol(backSPIN_clustered_data_frame)), .Label=rev(rownames(backSPIN_clustered_data_frame)), class="factor" )
  
  # -- Add the duplication of the genes names
  cell_names_repetition_vector=c()
  for(i in 1:ncol(backSPIN_clustered_data_frame)){
    cell_names_repetition_vector=c(cell_names_repetition_vector,rep(colnames(backSPIN_clustered_data_frame)[i],times=nrow(backSPIN_clustered_data_frame)))
  }
  #one_gene_per_line_df$Cells = factor(cell_names_repetition_vector, levels=unique(cell_names_repetition_vector))
  one_gene_per_line_df$UniqueCellID = structure( rep(1:ncol(backSPIN_clustered_data_frame),each=nrow(backSPIN_clustered_data_frame) ), .Label=colnames(backSPIN_clustered_data_frame), class="factor")
  
  
  # -- Add the value of expressions
  gene_et = vector()
  for( i in 1:ncol(backSPIN_clustered_data_frame)){
    gene_et = append( gene_et, unlist( backSPIN_clustered_data_frame[ , i ]))
  }
  maxEt=max(gene_et)
  for(i in 1:length(gene_et)){
    if(gene_et[i]==0){
      gene_et[i]=NA
    }
  }
  
  # split the continuous data to make discrete data
  gene_et_cut = cut(gene_et, breaks= seq(0,maxEt+1, by=(maxEt+1)/6),labels=1:6)
  
  one_gene_per_line_df$Et = gene_et_cut
  
  #data frame for the side bar colored by the phenotype
  side_bar_df = data.frame(count=1:ncol(backSPIN_clustered_data_frame))
  side_bar_df$UniqueCellID =  structure(1:ncol(backSPIN_clustered_data_frame),  .Label=colnames(backSPIN_clustered_data_frame), class="factor")
  side_bar_df$Y = rep(-1,times=ncol(backSPIN_clustered_data_frame))
  side_bar_df = merge(side_bar_df, raw_data_for_biclustering[,c("UniqueCellID",side_bar_colored_variable)],by="UniqueCellID", all=FALSE)
  
  biclust_ggplot<-ggplot() +
    geom_raster(data=one_gene_per_line_df,aes(x=UniqueCellID,y=Genes,fill=Et)) +
    geom_raster(data=side_bar_df,aes_string(x="UniqueCellID",y="Y",fill=side_bar_colored_variable)) + #facet_wrap(c(side_bar_colored_variable)) +
    scale_fill_manual(values=c(scales::seq_gradient_pal("#132B43","#56B1F7","Lab")(seq(0,1,length.out=6)), "gray18",gg_color_hue(length(unique(as.vector(t(side_bar_df[side_bar_colored_variable]))))-1) ),na.value="grey50",breaks=c("1","2","3","4","5","6",sort(unique(as.vector(t(side_bar_df[side_bar_colored_variable]))))), labels=c("0 to 4.3","4.4 to 8.7","8.8 to 13.0","13.1 to 17.3","17.4 to :21.6","21.7 to 26", sort(unique(as.vector(t(side_bar_df[side_bar_colored_variable]))))) )+
    geom_rect(size=0.5,fill=NA,colour="black",aes(xmin=index_first_group_cell-0.5,xmax=(index_second_group_cell-1)+0.5, ymin=(nrow(backSPIN_clustered_data_frame)-(index_second_group_gene-1)+1)-0.5,ymax=(nrow(backSPIN_clustered_data_frame)-index_first_group_gene+1)+0.5)) +
    geom_rect(size=0.5,fill=NA,colour="black",aes(xmin=index_second_group_cell-0.5,xmax=ncol(backSPIN_clustered_data_frame)+0.5, ymin=1-0.5,ymax=(nrow(backSPIN_clustered_data_frame)-index_second_group_gene+1)+0.5)) +
    #geom_rect(aes(xmin=500,xmax=1000, ymin=-1,ymax=-3)) +
    theme(axis.text.x  = element_blank(), axis.text.y  = element_text(size=7)) + #element_text(angle=90, size=5)
    ggtitle(title_ggplot)
  
  print(biclust_ggplot)
  
  
}

# #####################################################
# #####################################################
# #####################################################
# Pseudotime
# #####################################################
# #####################################################
# #####################################################
pseudotimeDDRTree_function <- function(raw_data_for_pseudotimeDDRTree, genes_to_study_v, cell_attribute_colorBy=c("GroupPheno"), cell_attribute_facetwrapBy=NA ,bool_facet_wrap=FALSE, DDRTree_max_components=2, plot_title=""){
  #DDRTree_max_components=2
  
  
  bool_change_cellAttributeFaceWrap=FALSE
  if(!is.na(cell_attribute_facetwrapBy)){
    bool_change_cellAttributeFaceWrap=TRUE
  }
  
  # creation of a cellData for the CellDateSet
  raw_data_transposed=t(raw_data_for_pseudotimeDDRTree[, genes_to_study_v])
  colnames(raw_data_transposed)=raw_data_for_pseudotimeDDRTree$UniqueCellID
  
  # creation of a phenoData for the CellDateSet
  phenoData_data_frame=data.frame(UniqueCellID=raw_data_for_pseudotimeDDRTree[,"UniqueCellID"])
  for(gene in genes_to_study_v){
    phenoData_data_frame[[gene]]=as.vector(raw_data_transposed[gene,])
  }
  for(i in 1:length(cell_attribute_colorBy)){
    if(!cell_attribute_colorBy[i] %in% c("State","Pseudotime")){
      phenoData_data_frame[[cell_attribute_colorBy[i] ]] = raw_data_for_pseudotimeDDRTree[,cell_attribute_colorBy[i]]
    }
  }
  
  row.names(phenoData_data_frame)=colnames(raw_data_transposed)
  phenoData_annotated <- new("AnnotatedDataFrame", data = phenoData_data_frame)
  
  
  #cds=new ("CellDataSet",exprs=raw_data_transposed)
  cds=newCellDataSet(as.matrix(raw_data_transposed),phenoData=phenoData_annotated, expressionFamily=gaussianff())
  
  ggplot_DDRTree_title=paste("DDRTree Pseudotime:",plot_title, "(max CP=",DDRTree_max_components,")")
  
  set.seed(45)
  # -- DDRTree reduce dimension
  cds_reduce_dim_DDRTree=reduceDimension(cds, max_components=DDRTree_max_components, reduction_method="DDRTree",norm_method="none",pseudo_expr=0) #DDRTree
  tree_MST_niaki_DDRTree=orderCells(cds_reduce_dim_DDRTree)
  
  # -- ggplot of the DDRTree
  
  for(cell_attribute_colorBy_value in cell_attribute_colorBy){
    if(bool_change_cellAttributeFaceWrap){
      cell_attribute_facetwrapBy=cell_attribute_colorBy_value
    }
    ggplot_cell_trajectory<-plot_cell_trajectory(tree_MST_niaki_DDRTree,cell_size=1, color_by=cell_attribute_colorBy_value)#+ facet_wrap(~State)#, color_by="SortPheno") #, markers="MKI67")
    ggplot_cell_trajectory<-ggplot_cell_trajectory + ggtitle(ggplot_DDRTree_title)
    ggplot_cell_trajectory<-ggplot_cell_trajectory+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),line = element_blank())
    #ggplot_cell_trajectory
    #print(ggplot_cell_trajectory)
    
    
    
    #grid.arrange(ggplot_cell_trajectory, ggplot_cell_trajectory_multi, nrow = 2, heights = c(0.5, 0.5))
    if(bool_facet_wrap==TRUE && !(identical(cell_attribute_colorBy_value,"Pseudotime")) ){
      ggplot_cell_trajectory_multi<-plot_cell_trajectory(tree_MST_niaki_DDRTree,cell_size=1, show_branch_points = FALSE, color_by=cell_attribute_colorBy_value)+ facet_wrap(cell_attribute_facetwrapBy)#, color_by="SortPheno") #, markers="MKI67")
      ggplot_cell_trajectory_multi<-ggplot_cell_trajectory_multi + ggtitle(ggplot_DDRTree_title)
      ggplot_cell_trajectory_multi<-ggplot_cell_trajectory_multi+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),line = element_blank())
      #ggplot_cell_trajectory
      #print(ggplot_cell_trajectory_multi)
      #grid.arrange(ggplot_cell_trajectory, ggplot_cell_trajectory_multi, nrow = 2, heights = c(0.5, 0.5))
      print(ggplot_cell_trajectory)
      print(ggplot_cell_trajectory_multi)
    }else{
      print(ggplot_cell_trajectory)
    }
  }
  
  # - recovery of the pseudotime values for each single-cell
  pseudotime_data_Monocle_DDRTree=t(pData(tree_MST_niaki_DDRTree)["Pseudotime"])
  
  # - recovery of the DDRTree state for each single-cell:
  DDRTree_state=t(pData(tree_MST_niaki_DDRTree)["State"])
  
  variables_to_return= list()
  variables_to_return[["DDRTree_pseudotime"]]=pseudotime_data_Monocle_DDRTree
  variables_to_return[["DDRTree_state"]]=DDRTree_state
  variables_to_return[["CellDataSet_DDRTree"]]=tree_MST_niaki_DDRTree
  
  return(variables_to_return)
}

# #####################################################
# #####################################################
# #####################################################
# Pseudotime CALCUL
# #####################################################
# #####################################################
# #####################################################
pseudotimeDDRTree_function_calcul <- function(raw_data_for_pseudotimeDDRTree, genes_to_study_v, DDRTree_max_components=3){
  #DDRTree_max_components=2
  
  
  
  # creation of a cellData for the CellDateSet
  raw_data_transposed=t(raw_data_for_pseudotimeDDRTree[, genes_to_study_v])
  colnames(raw_data_transposed)=raw_data_for_pseudotimeDDRTree$UniqueCellID
  
  # creation of a phenoData for the CellDateSet
  phenoData_data_frame=data.frame(UniqueCellID=rownames(raw_data_for_pseudotimeDDRTree))
  for(gene in genes_to_study_v){
    phenoData_data_frame[[gene]]=as.vector(raw_data_transposed[gene,])
  }
  cell_attribute_colorBy=colnames(raw_data_for_pseudotimeDDRTree)
  for(i in 1:length(cell_attribute_colorBy)){
    if(!cell_attribute_colorBy[i] %in% c("State","Pseudotime")){
      phenoData_data_frame[[cell_attribute_colorBy[i] ]] = raw_data_for_pseudotimeDDRTree[,cell_attribute_colorBy[i]]
    }
  }
  
  row.names(phenoData_data_frame)=colnames(raw_data_transposed)
  phenoData_annotated <- new("AnnotatedDataFrame", data = phenoData_data_frame)
  
  
  #cds=new ("CellDataSet",exprs=raw_data_transposed)
  cds=newCellDataSet(as.matrix(raw_data_transposed),phenoData=phenoData_annotated, expressionFamily=gaussianff()) #negbinomial.size())
  
  
  set.seed(45)
  # -- DDRTree reduce dimension
  cds_reduce_dim_DDRTree=reduceDimension(cds, max_components=DDRTree_max_components, reduction_method="DDRTree",norm_method="none",pseudo_expr=0) #DDRTree
  tree_MST_niaki_DDRTree=orderCells(cds_reduce_dim_DDRTree)
  
  # - recovery of the pseudotime values for each single-cell
  pseudotime_data_Monocle_DDRTree=t(pData(tree_MST_niaki_DDRTree)["Pseudotime"])
  
  # - recovery of the DDRTree state for each single-cell:
  DDRTree_state=t(pData(tree_MST_niaki_DDRTree)["State"])
  
  variables_to_return= list()
  variables_to_return[["DDRTree_pseudotime"]]=pseudotime_data_Monocle_DDRTree
  variables_to_return[["DDRTree_state"]]=DDRTree_state
  variables_to_return[["CellDataSet_DDRTree"]]=tree_MST_niaki_DDRTree
  
  return(variables_to_return)
}

# #####################################################
# #####################################################
# #####################################################
# Pseudotime PLOTING
# #####################################################
# #####################################################
# #####################################################
pseudotimeDDRTree_function_plot <- function(tree_MST_niaki_DDRTree, cell_attribute_colorBy=c("GroupPheno"), cell_attribute_facetwrapBy=NA ,bool_facet_wrap=FALSE,  plot_title="",DDRTree_max_components=3){
  
  
  bool_change_cellAttributeFaceWrap=FALSE
  if(!is.na(cell_attribute_facetwrapBy)){
    bool_change_cellAttributeFaceWrap=TRUE
  }
  
  ggplot_DDRTree_title=paste("DDRTree Pseudotime:",plot_title, "(max CP=",DDRTree_max_components,")")
  
  # -- ggplot of the DDRTree
  
  for(cell_attribute_colorBy_value in cell_attribute_colorBy){
    if(bool_change_cellAttributeFaceWrap){
      cell_attribute_facetwrapBy=cell_attribute_colorBy_value
    }
    ggplot_cell_trajectory<-plot_cell_trajectory(tree_MST_niaki_DDRTree,cell_size=1, color_by=cell_attribute_colorBy_value)#+ facet_wrap(~State)#, color_by="SortPheno") #, markers="MKI67")
    ggplot_cell_trajectory<-ggplot_cell_trajectory + ggtitle(ggplot_DDRTree_title)
    ggplot_cell_trajectory<-ggplot_cell_trajectory+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),line = element_blank())
    #ggplot_cell_trajectory
    #print(ggplot_cell_trajectory)
    
    
    
    #grid.arrange(ggplot_cell_trajectory, ggplot_cell_trajectory_multi, nrow = 2, heights = c(0.5, 0.5))
    if(bool_facet_wrap==TRUE && !(identical(cell_attribute_colorBy_value,"Pseudotime")) ){
      ggplot_cell_trajectory_multi<-plot_cell_trajectory(tree_MST_niaki_DDRTree,cell_size=1, show_branch_points = FALSE, color_by=cell_attribute_colorBy_value)+ facet_wrap(cell_attribute_facetwrapBy)#, color_by="SortPheno") #, markers="MKI67")
      ggplot_cell_trajectory_multi<-ggplot_cell_trajectory_multi + ggtitle(ggplot_DDRTree_title)
      ggplot_cell_trajectory_multi<-ggplot_cell_trajectory_multi+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),line = element_blank())
      #ggplot_cell_trajectory
      #print(ggplot_cell_trajectory_multi)
      #grid.arrange(ggplot_cell_trajectory, ggplot_cell_trajectory_multi, nrow = 2, heights = c(0.5, 0.5))
      print(ggplot_cell_trajectory)
      print(ggplot_cell_trajectory_multi)
    }else{
      print(ggplot_cell_trajectory)
    }
  }
  
  # - recovery of the pseudotime values for each single-cell
  pseudotime_data_Monocle_DDRTree=t(pData(tree_MST_niaki_DDRTree)["Pseudotime"])
  
  # - recovery of the DDRTree state for each single-cell:
  DDRTree_state=t(pData(tree_MST_niaki_DDRTree)["State"])
  
  variables_to_return= list()
  variables_to_return[["DDRTree_pseudotime"]]=pseudotime_data_Monocle_DDRTree
  variables_to_return[["DDRTree_state"]]=DDRTree_state
  variables_to_return[["CellDataSet_DDRTree"]]=tree_MST_niaki_DDRTree
  
  return(variables_to_return)
}

# #####################################################
# #####################################################
# #####################################################
# function LEGEND for heatmap.3 
# #####################################################
# #####################################################
# #####################################################
legend_for_heatmap.3 <- function(raw_data_df,
                                 ColSide_name_v,
                                 position_legend ="topright"){
  
  # deals of the column side bars
  # #############################
  column_annotation_matrix = matrix(ncol=length(ColSide_name_v),nrow=nrow(raw_data_df))
  legend_side_bars_names_v = c() #for the legend
  legend_side_bars_colors_v =c() #for the legend
  for(index_column_name in 1:length(ColSide_name_v)){
    colors_to_use = gg_color_hue(length(unique(raw_data_df[,ColSide_name_v[index_column_name] ])))
    names(colors_to_use) = sort(unique(raw_data_df[,ColSide_name_v[index_column_name] ]))
    sidebar_color_vector = raw_data_df[,ColSide_name_v[index_column_name] ]
    sidebar_color_vector = sapply(sidebar_color_vector, function(xx) colors_to_use[xx] )
    column_annotation_matrix[,index_column_name]=as.vector(sidebar_color_vector)
    legend_side_bars_names_v = c( legend_side_bars_names_v, names(colors_to_use), "")
    legend_side_bars_colors_v = c( legend_side_bars_colors_v, colors_to_use, "white")
  }
  colnames(column_annotation_matrix) <- ColSide_name_v
  ColSideColors = column_annotation_matrix
  
  legend(position_legend,legend=legend_side_bars_names_v, fill=legend_side_bars_colors_v, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
}


# #####################################################
# #####################################################
# #####################################################
# function heatmap.3 used to generate a heatmap.2 with 
# multiple color side bars
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
# #####################################################
# #####################################################
# #####################################################
heatmap.3 <- function(raw_data_df,
                      genes_v,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      #ColSideColors,
                      ColSide_name_v,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      titleColorGradient = "Value",
                      KeyValueName="",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x = t(raw_data_df[,genes_v])
  
  # deals of the column side bars
  # #############################
  column_annotation_matrix = matrix(ncol=length(ColSide_name_v),nrow=nrow(raw_data_df))
  legend_side_bars_names_v = c() #for the legend
  legend_side_bars_colors_v =c() #for the legend
  for(index_column_name in 1:length(ColSide_name_v)){
    colors_to_use = gg_color_hue(length(unique(raw_data_df[,ColSide_name_v[index_column_name] ])))
    names(colors_to_use) = sort(unique(raw_data_df[,ColSide_name_v[index_column_name] ]))
    sidebar_color_vector = raw_data_df[,ColSide_name_v[index_column_name] ]
    sidebar_color_vector = sapply(sidebar_color_vector, function(xx) colors_to_use[xx] )
    column_annotation_matrix[,index_column_name]=as.vector(sidebar_color_vector)
    legend_side_bars_names_v = c( legend_side_bars_names_v, names(colors_to_use), "")
    legend_side_bars_colors_v = c( legend_side_bars_colors_v, colors_to_use, "white")
  }
  colnames(column_annotation_matrix) <- ColSide_name_v
  ColSideColors = column_annotation_matrix
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(xpd=TRUE,mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(xpd=TRUE,mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(xpd=TRUE,mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol, font=3) #inaki modification
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow, font=3) #inaki modification
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex, font=3) #inaki
  par(xpd=TRUE,mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(xpd=TRUE,mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    par(xpd=TRUE)
    #legend(nrow(raw_data_df)+10,200,legend=legend_side_bars_names_v, fill=legend_side_bars_colors_v, border=FALSE, bty="n", cex=1)
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(xpd=TRUE,mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else {
      title(titleColorGradient) #inaki
      
    }
    
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
  #par(xpd=TRUE)
  #legend(1,0,legend=legend_side_bars_names_v, fill=legend_side_bars_colors_v, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
}

# #####################################################
# #####################################################
# #####################################################
# isOutlier_MAD_function same function as "isOutlier" in scater package 
# #####################################################
# #####################################################
# #####################################################
isOutlier_MAD_function <- function (x_v,  nmads=2, type="both", log=TRUE){
  vector_v=x_v
  if(log==TRUE){
    vector_v=log10(x_v)
  }
  mad_value = mad(vector_v)
  bool_v=c()
  if(type=="both"){
    bool_v = c( (vector_v < (median(vector_v) - nmads*mad_value)) | (vector_v > (median(vector_v) + nmads*mad_value)) )
  } else if(type=="lower"){
    bool_v = c( vector_v < (median(vector_v) - nmads*mad_value) )
  } else if(type=="higher"){
    bool_v = c( vector_v > (median(vector_v) + nmads*mad_value) )
  }
  return(bool_v)
}