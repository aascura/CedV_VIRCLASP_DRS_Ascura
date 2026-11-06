#!/bin/bash

# Title: Oxford Modkit DRS Analysis - Compress and Indexing Modkit Pileup
# Purpose: Compress and Tabix Index Bedmethyl Files to Prep for Merging in Modkit_Merge_BedMethyl.sh
# Reference: https://nanoporetech.github.io/modkit/intro_dmr.html
# Author: Andrea Ascura
# Created: 05Sep2025
# Last Update: 04Nov2026

# SPEN Compress and Tabix Index
# 1 - Specify location of bedmethyl outputs from Modkit_Modbam_Pileup.sh (Input)
PILEUP_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/pileup_output_V2"

# 2 - Prepare input data by compressing and indexing bedMethyl files produced from modkit pileup
bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_control_1.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_control_1.bed.gz"

bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_control_2.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_control_2.bed.gz"

bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_spen_1.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_spen_1.bed.gz"

bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_spen_2.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.9_cedv_spen_2.bed.gz"

# STM2457 Compress and Tabix Index
# 1 - Specify location of bedmethyl outputs from Modkit_Modbam_Pileup.sh (Input)
PILEUP_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/pileup_output"

# 2 - Prepare input data by compressing and indexing bedMethyl files produced from modkit pileup
bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed.gz"

bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed.gz"

bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_1.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_1.bed.gz"

bgzip -k "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_2.bed"
tabix -p bed "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_2.bed.gz"