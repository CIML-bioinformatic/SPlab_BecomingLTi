#!/bin/bash

# This script aims to execute the velocyto program on LysoDC data.
# It has to be used in the Docker scrnaseq_rnavelocity docker container
# 

# Provide the path to the bam file produced by CellRanger count
export BAM_FILE_PATH=$1
# Provide the path to the genome annotation file
export GENOME_ANNOTATION_PATH=$2
# Provide the folder where to put the Velocyto results (loom file)
export OUTPUT_PATH=$(dirname "$3")
# The sample ID
export SAMPLE_ID=11230812

echo BAM_FILE_PATH=$BAM_FILE_PATH
echo GENOME_ANNOTATION_PATH=$GENOME_ANNOTATION_PATH
echo OUTPUT_PATH=$OUTPUT_PATH

# Ensure Samtools is in the PATH
export PATH=/samtools/bin:$PATH

# Create output folder
mkdir -p $OUTPUT_PATH

# Execute Velocyto on 10x data
velocyto run --outputfolder $OUTPUT_PATH --sampleid $SAMPLE_ID $BAM_FILE_PATH $GENOME_ANNOTATION_PATH

