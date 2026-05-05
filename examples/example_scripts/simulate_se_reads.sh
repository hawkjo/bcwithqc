#!/bin/bash

# Set the variables
CONFIG="../simulate_se_mini_config.json"
STAR_INDEX="../se_mini_genome_index"

SIMULATE_OUTPUT_DIR="./test_run/simulate_se/simulated_fastqs"
PREPROCESS_OUTPUT_DIR="./test_run/simulate_se/preprocess_output"
STAR_OUTPUT_DIR="./test_run/simulate_se/STAR_output/" #Trailing slash is required
COUNT_OUTPUT_DIR="./test_run/simulate_se/count_output"
# PREPROCESS_OUTPUT_DIR and COUNT_OUTPUT_DIR are allowed to be the same directory, and STAR_OUTPUT_DIR can be a subdirectory of it. 

N_THREADS=1

mkdir -p $SIMULATE_OUTPUT_DIR $PREPROCESS_OUTPUT_DIR $STAR_OUTPUT_DIR $COUNT_OUTPUT_DIR

# Simulate SE reads
bcwithqc simulate_reads \
--config=$CONFIG \
--output-dir=$SIMULATE_OUTPUT_DIR \
--nreads=10000 \
--error-probability=0.1 \
--substitution-probability=0.7 \
-vvv

# Preprocess the simulated reads
bcwithqc preprocess \
$SIMULATE_OUTPUT_DIR \
--config=$CONFIG \
--output-dir=$PREPROCESS_OUTPUT_DIR \
--threads=$N_THREADS \
-vvv

# Align with STAR
STAR \
--runThreadN $N_THREADS \
--genomeDir $STAR_INDEX \
--readFilesIn $PREPROCESS_OUTPUT_DIR/sans_bc_se_mini.fq \
--outFileNamePrefix $STAR_OUTPUT_DIR \
--outFilterMultimapNmax 1 \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI AS nM GX GN

# Call bcwithqc count for the final results
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
