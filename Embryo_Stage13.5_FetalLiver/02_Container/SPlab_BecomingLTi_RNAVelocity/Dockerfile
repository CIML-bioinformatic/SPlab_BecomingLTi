# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rocker/tidyverse:3.4

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
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite( "patchwork")'  # +/ operators for ggplots

# plotly
RUN Rscript -e 'install.packages( "plotly")'

# general
RUN Rscript -e 'install.packages( "gplots")'         # heatmap.2
RUN Rscript -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite( "iheatmapr")'  # iheatmap (interactive, uses plotly), dependencies OK with BiocManager
RUN Rscript -e 'install.packages( "pheatmap")'       # pheatmap

#### Reporting
RUN Rscript -e 'install.packages( "DT")'             # datatable
RUN Rscript -e 'install.packages( "pander")'         # pander

#### General
RUN Rscript -e 'install.packages( "funr")'           # get_script_path
RUN Rscript -e 'install.packages( "reshape")'        # melt

#### package Seurat
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite( "Seurat")'
RUN Rscript -e 'install.packages( "umap")'

#### Package velocyto
#RUN apt-get update && apt-get install -y libhdf5-dev libboost-all-dev
#RUN Rscript -e 'install.packages("BiocManager")'
#RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocGenerics")'
#RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("pcaMethods", ask=FALSE)'
#RUN Rscript -e 'install.packages( "h5")'
#RUN Rscript -e 'library(devtools);install_github("velocyto-team/velocyto.R")'

#### Package Pagoda2
#RUN Rscript -e 'library(devtools);install_github("igraph/rigraph")'
#RUN Rscript -e 'library(devtools);install_github("jkrijthe/Rtsne",ref="openmp")'
#RUN Rscript -e 'install.packages(c("Cairo","urltools"))'
#RUN Rscript -e 'library(devtools);install_github("hms-dbmi/pagoda2")'

RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  libbamtools-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  libxml2-dev \
  libz-dev \
  curl \
  libhdf5-cpp-100 \ 
  libarmadillo7 \
  libarmadillo-dev

RUN Rscript -e 'chooseCRANmirror(ind=52); install.packages(c("devtools", "Rcpp","RcppArmadillo", "Matrix", "mgcv", "abind","igraph","h5","Rtsne","cluster","data.table"))'
RUN Rscript -e 'install.packages( "dplyr")'
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("pcaMethods","edgeR","Rsamtools","GenomicAlignments","GenomeInfoDb","Biostrings"),suppressAutoUpdate=TRUE,ask=FALSE,suppressUpdates=TRUE)'

RUN Rscript -e 'library( devtools); install_github("velocyto-team/velocyto.R")'

#RUN mkdir /velocyto
#WORKDIR /velocyto

#RUN \
#  git clone https://github.com/velocyto-team/velocyto.R && \
#  mkdir -p /velocyto/R/x86_64-redhat-linux-gnu-library/3.4

#RUN \
#  echo '.libPaths(c("/velocyto/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths()))' > .Rprofile && \
#  R -e 'devtools::install_local("/velocyto/velocyto.R/",dep=T,upgrade_dependencies=F)'
#
#ENV  LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH \
#  R_PROFILE=/velocyto/.Rprofile



# #################################################
# DECLARE EXPECTED PORTS TO BE EXPOSED
# #################################################

# Rstudio listens on port 8787 by default
# Container typically started with `-p 8787:8787`
EXPOSE 8787





