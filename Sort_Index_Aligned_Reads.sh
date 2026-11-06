#!/bin/bash

# Title: Oxford Dorado DRS Analysis - BAM Sorting and Indexing Aligned Reads using Samtools (Post-processing Alignment)
# Purpose: Sort and index aligned reads using samtools for both SPEN and STM2457 experiments
# Author: Andrea Ascura
# Created: 03Sep2025
# Last Update: 04Nov2026

# SPEN sorting and indexing
# Specify location of aligned reads output from Dorado_Minimap2_Alignment.sh (Input)
ALIGNED_OUTPUT_SPEN="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/SPEN/aligned_output"

# Specify location of sorted and indexed reads generated from this script (Output)
SORTED_INDEXED_OUTPUT_SPEN="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/SPEN/sorted_indexed_output"

# 1 - Sort alignments by reference position to create a sorted BAM
samtools sort -o "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_control_1.bam" "${ALIGNED_OUTPUT_SPEN}/called_control_1.bam"
samtools sort -o "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_control_2.bam" "${ALIGNED_OUTPUT_SPEN}/called_control_2.bam"
samtools sort -o "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_spen_1.bam" "${ALIGNED_OUTPUT_SPEN}/called_spen_1.bam"
samtools sort -o "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_spen_2.bam" "${ALIGNED_OUTPUT_SPEN}/called_spen_2.bam"

# 2 - Indexed sorted alignments
samtools index "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_control_1.bam"
samtools index "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_control_2.bam"
samtools index "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_spen_1.bam"
samtools index "${SORTED_INDEXED_OUTPUT_SPEN}/sorted_txn_aligned_cedv_spen_2.bam"

# STM2457 sorting and indexing
# Specify location of aligned reads output from Dorado_Minimap2_Alignment.sh (Input)
ALIGNED_OUTPUT_STM="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/STM/aligned_output"

# Specify location of sorted and indexed reads generated from this script (Output)
SORTED_INDEXED_OUTPUT_STM="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/STM/sorted_indexed_output"

# 1 - Sort alignments by reference position to create a sorted BAM
samtools sort -o "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_control_1.bam" "${ALIGNED_OUTPUT_STM}/calls_control_1.bam"
samtools sort -o "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_control_2.bam" "${ALIGNED_OUTPUT_STM}/calls_control_2.bam"
samtools sort -o "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_stm_1.bam" "${ALIGNED_OUTPUT_STM}/calls_stm_1.bam"
samtools sort -o "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_stm_2.bam" "${ALIGNED_OUTPUT_STM}/calls_stm_2.bam"

# 2 - Indexed sorted alignments
samtools index "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_control_1.bam"
samtools index "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_control_2.bam"
samtools index "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_stm_1.bam"
samtools index "${SORTED_INDEXED_OUTPUT_STM}/sorted_txn_aligned_cedv_stm_2.bam"
