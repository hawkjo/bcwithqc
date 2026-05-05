#!/bin/bash

# Set the variables
CDNA_INPUT="../gDNA_fastqs_paired_end"
CDNA_CONFIG="../gDNA_paired_end_config.json"
STAR_INDEX="../SDR001_REF_index"

PREPROCESS_OUTPUT_DIR="./test_run/gDNA_paired_end/preprocess_output"
STAR_OUTPUT_DIR="./test_run/gDNA_paired_end/STAR_output/" #Trailing slash is required
COUNT_OUTPUT_DIR="./test_run/gDNA_paired_end/count_output"
# PREPROCESS_OUTPUT_DIR and COUNT_OUTPUT_DIR are allowed to be the same directory, and STAR_OUTPUT_DIR can be a subdirectory of it. 

N_THREADS=1

mkdir -p "$PREPROCESS_OUTPUT_DIR" "$STAR_OUTPUT_DIR" "$COUNT_OUTPUT_DIR"

# Preprocess gDNA paired-end
bcwithqc preprocess \
"$CDNA_INPUT" \
--config="$CDNA_CONFIG" \
--output-dir="$PREPROCESS_OUTPUT_DIR" \
--threads="$N_THREADS" \
-vvv

# Align gDNA with STAR
STAR \
--runThreadN "$N_THREADS" \
--genomeDir "$STAR_INDEX" \
--readFilesIn "$PREPROCESS_OUTPUT_DIR/sans_bc_gDNA_1_sequence.fq" "$PREPROCESS_OUTPUT_DIR/sans_bc_gDNA_2_sequence.fq" \
--outFileNamePrefix "$STAR_OUTPUT_DIR" \
--outFilterMultimapNmax 1 \
--outSAMtype BAM Unsorted \
--outSAMattributes NH HI AS nM GX GN

# Count the results for gDNA
bcwithqc count \
"$PREPROCESS_OUTPUT_DIR" \
--STAR-output-dir="$STAR_OUTPUT_DIR" \
--config="$CDNA_CONFIG" \
--output-dir="$COUNT_OUTPUT_DIR" \
--threads="$N_THREADS" \
--keep-intermediary \
-vvv

# remove keep-intermediary for automatic cleanup. 
# Will remove STAR alignment results if they are in a subdirectory of COUNT_OUTPUT_DIR