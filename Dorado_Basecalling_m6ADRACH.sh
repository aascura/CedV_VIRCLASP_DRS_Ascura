#!/bin/bash

# Title: Oxford Dorado DRS Analysis - Simplex Basecalling with Modified Bases
# Purpose: Use Dorado to basecall Oxford Nanopore DRS raw data (POD5 files) for SPEN and STM2457 experiments. 
# Author: Andrea Ascura
# Updated: 04Nov2026

# Basecall Control_1, SPEN_1, Control_2, and SPEN_2 POD5 files.
# Input directory - Specify location of POD5 files
POD5_CONTROL_1="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/074_CedV_mRNA_Rep_1_merged_Control_SPEN/Control_1"
POD5_SPEN_1="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/074_CedV_mRNA_Rep_1_merged_Control_SPEN/SPEN_1"
POD5_CONTROL_2="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/074_CedV_mRNA_Rep_2_Control_SPEN/Control_2/20250813_1329_P2S-02001-A_PBC38893_14fae2ab/pod5"
POD5_SPEN_2="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/074_CedV_mRNA_Rep_2_Control_SPEN/SPEN_2/20250813_1329_P2S-02001-B_PBC48374_4dabf104/pod5"

# Output directory - Specify location of basecalled files (these files will be in the input for alignment)
BASECALLED_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/basecalled_output"

# Run Dorado basecaller - Simplex Basecalling with Modified Bases
dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_CONTROL_1" \
    --device cuda:0 \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/calls_control_1.bam"

dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_SPEN_1" \
    --device cuda:all \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/ calls_spen_1.bam"

dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_CONTROL_2" \
    --device cuda:all \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/calls_control_2.bam"

dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_SPEN_2" \
    --device cuda:all \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/calls_spen_2.bam"


# Basecall Control_1, STM_1, Control_2, and STM_2 POD5 files. STM = STM2457
# Input directory - Specify location of POD5 files
POD5_CONTROL_1="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Control/20260416_1722_P2S-02001-A_PAY38847_4aa25d20/pod5"
POD5_STM_1="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/STM/20260416_1722_P2S-02001-B_PAY27500_18e73064/pod5"
POD5_CONTROL_2="/home/ascanolab/Ascura/DirectRNAseq_Ascura/260421_Control_STM2457_Rep2/Control/20260421_1647_P2S-02001-A_PAY38847_0fd7357e/pod5"
POD5_STM_2="/home/ascanolab/Ascura/DirectRNAseq_Ascura/260421_Control_STM2457_Rep2/STM2457/20260421_1647_P2S-02001-B_PAY27500_98b03d1e/pod5"

# Output directory - Specify location of basecalled files (these files will be in the input for alignment)
BASECALLED_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/basecalled_output"
BASECALLED_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/basecalled_output"

# Run Dorado basecaller - Simplex Basecalling with Modified Bases
dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_CONTROL_1" \
    --device cuda:all \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/calls_control_1.bam"

dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_STM_1" \
    --device cuda:all \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/calls_stm_1.bam"

dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_CONTROL_2" \
    --device cuda:all \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/calls_control_2.bam"

dorado basecaller sup@v5.2.0,m6A_DRACH@v1 "$POD5_STM_2" \
    --device cuda:all \
    --min-qscore 7 \
    --output-dir "${BASECALLED_OUTPUT}/calls_stm_2.bam"




