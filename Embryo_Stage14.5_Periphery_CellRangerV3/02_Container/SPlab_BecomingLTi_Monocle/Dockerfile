# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rocker/tidyverse:3.6.1

MAINTAINER RFenouil (fenouil@ciml.univ-mrs.fr)




# ##################
# OTHER DEPENDENCIES
# ##################

#### Python
# pip
RUN apt-get update && apt-get install -y \
    python3-pip \
    && pip3 install --upgrade --user pip

# Umap processing library (python)
RUN python3 -m pip install 'umap-learn==0.3.10'


# ##########
# R PACKAGES 
# ##########

#### Figures & layout management
# ggplot2
RUN Rscript -e 'install.packages( "ggplot2")'
RUN Rscript -e 'install.packages( "cowplot")'        # plot_grid, themes, ...
RUN Rscript -e 'install.packages( "ggpubr")'         # add_summary, geom_signif, ...
RUN Rscript -e 'install.packages( "ggrepel")'        # geom_text_repel, geom_label_repel
RUN Rscript -e 'install.packages( "gridExtra")'      # grid.arrange, ...
RUN Rscript -e 'BiocManager::install( "patchwork")'  # +/ operators for ggplots

# plotly
RUN Rscript -e 'install.packages( "plotly")'

# general
RUN Rscript -e 'install.packages( "gplots")'         # heatmap.2
RUN Rscript -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN Rscript -e 'BiocManager::install( "iheatmapr")'  # iheatmap (interactive, uses plotly), dependencies OK with BiocManager
RUN Rscript -e 'install.packages( "pheatmap")'       # pheatmap


#### Reporting
RUN Rscript -e 'install.packages( "DT")'             # datatable
RUN Rscript -e 'install.packages( "pander")'         # pander


#### General
RUN Rscript -e 'install.packages( "funr")'           # get_script_path
RUN Rscript -e 'install.packages( "reshape")'        # melt

RUN Rscript -e 'install.packages( "umap")'

RUN  apt-get update --fix-missing \ 
  && export DEBIAN_FRONTEND=noninteractive \
  && apt-get -y install libx11-dev freeglut3 freeglut3-dev libpng-dev libxml2-dev libssl-dev libcurl4-openssl-dev

RUN sudo apt-get update \
  && apt-get -y install build-essential cmake gsl-bin libgsl0-dev libeigen3-dev libboost-all-dev libssl-dev libcurl4-openssl-dev libssl-dev libcairo2-dev libxt-dev libgtk2.0-dev libcairo2-dev xvfb xauth xfonts-base

RUN Rscript -e 'install.packages( "ade4")'

RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install("pheatmap", ask=FALSE)'
RUN Rscript -e 'BiocManager::install("monocle", ask=FALSE)'
RUN Rscript -e 'BiocManager::install("BiocGenerics", ask=FALSE)'

# Package Pagoda2
RUN Rscript -e 'library(devtools);install_github("igraph/rigraph")'
RUN Rscript -e 'library(devtools);install_github("jkrijthe/Rtsne",ref="openmp")'
RUN Rscript -e 'install.packages(c("Cairo","urltools"))'
#RUN Rscript -e 'library(devtools);install_github("hms-dbmi/pagoda2")'
RUN Rscript -e 'install.packages( "d3heatmap")'




# #################################################
# DECLARE EXPECTED PORTS TO BE EXPOSED
# #################################################

# Rstudio listens on port 8787 by default
# Container typically started with `-p 8787:8787`
EXPOSE 8787





