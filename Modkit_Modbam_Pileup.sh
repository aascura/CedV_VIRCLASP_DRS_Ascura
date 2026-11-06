#!/bin/bash

# Title: Oxford Modkit DRS Analysis - Pileup of Modified bases
# Purpose: Tabulate counts of base modifications from every sequencing read over each aligned reference genomic position and apply modification threshold
# Guide Reference: https://nanoporetech.github.io/modkit/quick_start.html
# Author: Andrea Ascura
# Created: 04Sep2025
# Last Update: 04Nov2026

# SPEN Modkit Pileup
# 1 - Specify location of aligned, sorted and indexed modBAM files output from Sort_Index_Aligned_Reads.sh (Input)
SORTED_INDEXED_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/sorted_indexed_output"

# 2 - Specify location of desired output location (pileup bed files)
PILEUP_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/pileup_output"

# 3 - Create bedmethyl file to tabulate counts of base modifications (base modification threshold = 0.7)
modkit pileup --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_control_1.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed" --log-filepath pileup_bmt0.7_cedv_control_1.log
modkit pileup --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_control_2.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed" --log-filepath pileup_bmt0.7_cedv_control_2.log
modkit pileup --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_spen_1.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_spen_1.bed" --log-filepath pileup_bmt0.7_cedv_spen_1.log
modkit pileup --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_spen_2.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_spen_2.bed" --log-filepath pileup_bmt0.7_cedv_spen_2.log

# STM2457 Modkit Pileup
# 1 - Specify location of aligned, sorted and indexed modBAM files output from Sort_Index_Aligned_Reads.sh (Input)
SORTED_INDEXED_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/sorted_indexed_output"

# 2 - Specify location of desired output location (pileup bed files)
PILEUP_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/pileup_output"

# 3 - Create bedmethyl file to tabulate counts of base modifications (base modification threshold = 0.7)
modkit pileup --filter-threshold A:0.9 --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_control_1.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed" --log-filepath "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.log"
modkit pileup --filter-threshold A:0.9 --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_control_2.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed" --log-filepath "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.log"
modkit pileup --filter-threshold A:0.9 --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_stm_1.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_1.bed" --log-filepath "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_1.log"
modkit pileup --filter-threshold A:0.9 --mod-thresholds a:0.7 "${SORTED_INDEXED_OUTPUT}/sorted_aligned_cedv_stm_2.bam" "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_2.bed" --log-filepath "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_2.log"