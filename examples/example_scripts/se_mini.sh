#!/bin/bash

# Set the variables
INPUT_DIR="../se_mini"
CONFIG="../se_mini_config.json"
STAR_INDEX="../se_mini_genome_index"

PREPROCESS_OUTPUT_DIR="./test_run/se_mini/preprocess_output"
STAR_OUTPUT_DIR="./test_run/se_mini/STAR_output/" #Trailing slash is required
COUNT_OUTPUT_DIR="./test_run/se_mini/count_output"
# PREPROCESS_OUTPUT_DIR and COUNT_OUTPUT_DIR are allowed to be the same directory, and STAR_OUTPUT_DIR can be a subdirectory of it. 

mkdir -p "$PREPROCESS_OUTPUT_DIR" "$STAR_OUTPUT_DIR" "$COUNT_OUTPUT_DIR"

# Preprocess the SE mini data
bcwithqc preprocess \
"$INPUT_DIR" \
--config="$CONFIG" \
--output-dir="$PREPROCESS_OUTPUT_DIR" \
-v

# Align with STAR
STAR \
--runThreadN 1 \
--genomeDir "$STAR_INDEX" \
--readFilesIn "$PREPROCESS_OUTPUT_DIR/sans_bc_se_mini.fq" \
--outFileNamePrefix "$STAR_OUTPUT_DIR" \
--outFilterMultimapNmax 1 \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI AS nM GX GN

# Count the results
bcwithqc count \
"$PREPROCESS_OUTPUT_DIR" \
--STAR-output-dir="$STAR_OUTPUT_DIR" \
--config="$CONFIG" \
--output-dir="$COUNT_OUTPUT_DIR" \
--keep-intermediary \
-v

# remove keep-intermediary for automatic cleanup. 
# Will remove STAR alignment results if they are in a subdirectory of COUNT_OUTPUT_DIR