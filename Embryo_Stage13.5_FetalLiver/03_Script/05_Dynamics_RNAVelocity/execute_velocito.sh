#!/bin/bash

# This script aims to execute the velocyto program on LysoDC data.
# It has to be used in the Docker scrnaseq_rnavelocity docker container
# 

# Provide the folder where the 10x data are stored
export CELL_RANGER_PATH=/mnt/NAS7/SPlab/BIOINFO_PROJECT/BecomingLTi/Embryo_Stage13.5_FetalLiver/00_RawData/CellRanger_Counts/outs
# Provide the path to the genome annotation file
export GENOME_ANNOTATION_PATH=/mnt/NAS7/SPlab/BIOINFO_PROJECT/BecomingLTi/Embryo_Stage13.5_FetalLiver/01_Reference/Genome/GRCm38/annotation/Mus_musculus.GRCm38.90.gtf
# Provide the folder where to put the Velocyto results (loom file)
export OUTPUT_PATH=/mnt/NAS7/SPlab/BIOINFO_PROJECT/BecomingLTi/Embryo_Stage13.5_FetalLiver/05_Output/05_Dynamics_RNAVelocity
# The sample ID
export SAMPLE_ID=11230811


# Ensure Samtools is in the PATH
export PATH=/samtools/bin:$PATH

# Create output folder
mkdir -p $OUTPUT_PATH

# Execute Velocyto on 10x data
cd $OUTPUT_PATH
velocyto run -o $OUTPUT_PATH -e $SAMPLE_ID -c False -U False -u "no" $CELL_RANGER_PATH/possorted_genome_bam.bam $GENOME_ANNOTATION_PATH

