#!/bin/bash

# Set the variables
INPUT_DIR="../pe_mini"
CONFIG="../pe_mini_config.json"
STAR_INDEX="../pe_mini_genome_index"

PREPROCESS_OUTPUT_DIR="./test_run/pe_mini/preprocess_output"
STAR_OUTPUT_DIR="./test_run/pe_mini/STAR_output" 
COUNT_OUTPUT_DIR="./test_run/pe_mini/count_output"
# PREPROCESS_OUTPUT_DIR and COUNT_OUTPUT_DIR are allowed to be the same directory, and STAR_OUTPUT_DIR can be a subdirectory of it. 

N_THREADS=2

mkdir -p $PREPROCESS_OUTPUT_DIR $STAR_OUTPUT_DIR $COUNT_OUTPUT_DIR

# Preprocess the SE mini data
bcwithqc preprocess \
$INPUT_DIR \
--config=$CONFIG \
--output-dir=$PREPROCESS_OUTPUT_DIR \
--threads=$N_THREADS \
-vvv

# Align first pair with STAR
STAR \
--runThreadN $N_THREADS \
--genomeDir $STAR_INDEX \
--readFilesIn $PREPROCESS_OUTPUT_DIR/sans_bc_pe_mini_1_r1.fq $PREPROCESS_OUTPUT_DIR/sans_bc_pe_mini_1_r2.fq \
--outFileNamePrefix $STAR_OUTPUT_DIR/pe_mini_pair1_ \
--outFilterMultimapNmax 1 \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI AS nM GX GN

# Align second pair with STAR
STAR \
--runThreadN $N_THREADS \
--genomeDir $STAR_INDEX \
--readFilesIn $PREPROCESS_OUTPUT_DIR/sans_bc_pe_mini_2_r1.fq $PREPROCESS_OUTPUT_DIR/sans_bc_pe_mini_2_r2.fq \
--outFileNamePrefix $STAR_OUTPUT_DIR/pe_mini_pair2_ \
--outFilterMultimapNmax 1 \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI AS nM GX GN

# Count the results
bcwithqc count \
$PREPROCESS_OUTPUT_DIR \
--STAR-output-dir=$STAR_OUTPUT_DIR \
--config=$CONFIG \
--output-dir=$COUNT_OUTPUT_DIR \
--threads=$N_THREADS \
--keep-intermediary \
-vvv

# remove keep-intermediary for automatic cleanup. 
# Will remove STAR alignment results if they are in a subdirectory of COUNT_OUTPUT_DIR