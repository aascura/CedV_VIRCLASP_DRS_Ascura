#!/bin/bash

# Title: NanoCount - Estimation of transcript abundance and calculation of TPM (transcripts per million) for SPEN and STM2457 samples
# Author: Andrea Ascura
# Created: 09Jun2026
# Last Update: 10Jun2026

# SPEN samples NanoCount processing
# Specify location of aligned reads output from Dorado_Minimap2_Alignment.sh (Input)
SORTED_INDEXED_OUTPUT_SPEN="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/SPEN/sorted_indexed_output"

# Specify location of sorted and indexed reads generated from this script (Output)
NANOCOUNT_OUTPUT_SPEN="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/SPEN/nanocount_output"

# Execute nanocount to generate table containing transcript counts and TPM
NanoCount -i "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_control_1.bam" -o "${NANOCOUNT_OUTPUT_SPEN}/transcript_counts_cedv_control_1.tsv"
NanoCount -i "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_control_2.bam" -o "${NANOCOUNT_OUTPUT_SPEN}/transcript_counts_cedv_control_2.tsv"
NanoCount -i "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_spen_1.bam" -o "${NANOCOUNT_OUTPUT_SPEN}/transcript_counts_cedv_spen_1.tsv"
NanoCount -i "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_spen_2.bam" -o "${NANOCOUNT_OUTPUT_SPEN}/transcript_counts_cedv_spen_2.tsv"

# STM samples NanoCount processing
# Specify location of aligned reads output from Dorado_Minimap2_Alignment.sh (Input)
SORTED_INDEXED_OUTPUT_STM="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/STM/sorted_indexed_output"

# Specify location of sorted and indexed reads generated from this script (Output)
NANOCOUNT_OUTPUT_STM="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/STM/nanocount_output"

# Execute nanocount to generate table containing transcript counts and TPM
NanoCount -i "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_control_1.bam" -o "${NANOCOUNT_OUTPUT_STM}/transcript_counts_cedv_control_3.tsv"
NanoCount -i "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_control_2.bam" -o "${NANOCOUNT_OUTPUT_STM}/transcript_counts_cedv_control_4.tsv"
NanoCount -i "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_stm_1.bam" -o "${NANOCOUNT_OUTPUT_STM}/transcript_counts_cedv_stm_1.tsv"
NanoCount -i "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_stm_2.bam" -o "${NANOCOUNT_OUTPUT_STM}/transcript_counts_cedv_stm_2.tsv"

