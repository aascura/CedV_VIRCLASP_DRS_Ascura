#!/bin/bash

# Title: Oxford Modkit DRS Analysis - Differential Methylation Scoring of Individual Replicates
# Purpose: Perform differential methylation analysis (WT vs SPEN KD/METTL3i) between biological replicates and sample types at single base positions with basecalled A (aka without region segmentation) for downstream correlation of individual replicates
# Reference: https://nanoporetech.github.io/modkit/intro_dmr.html
# Author: Andrea Ascura
# Created: 05Sep2025
# Last Update: 10Nov2025

# SPEN Individual Replicates
# 1 - Specify location of aligned (indexed and sorted) modBAMs
PILEUP_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/m6A_DRACH/pileup_output_V2"

# 2 - Specify location of reference genomes
#FASTA_HG38="/home/ascanolab/seq_pipelines/parpipe2_v0.2/genome/GRCh38.primary_assembly.genome.fa"
FASTA_CEDV="/home/ascanolab/Ascura/DirectRNAseq_Ascura/reference_genomes/GCF_000924595.1_ViralProj264921_genomic.fna"

# 3 - Specify location of results - output
DMR_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/m6A_DRACH/dmr_output/individual_replicates"

# 4 - Perform differential methylation analysis between biological replicates and identical sample conditions
# Control_1 vs Control_2 OR SPEN_1 vs SPEN_2
modkit dmr pair \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed.gz" \
  -o "${DMR_OUTPUT}/cedv_bmt0.7_dmr_control_1_control_2.bed" \
  --ref "$FASTA_CEDV" \
  --base A \
  --threads 8 \
  --log-filepath "${DMR_OUTPUT}/cedv_dmr_control_replicates.log"

modkit dmr pair \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_spen_1.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_spen_2.bed.gz" \
  -o "${DMR_OUTPUT}/cedv_bmt0.7_dmr_spen_1_spen_2.bed" \
  --ref "$FASTA_CEDV" \
  --base A \
  --threads 8 \
  --log-filepath "${DMR_OUTPUT}/cedv_dmr_spen_replicates.log"  

# STM2457 Individual Replicates
# 1 - Specify location of aligned (indexed and sorted) modBAMs
PILEUP_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/pileup_output"

# 2 - Specify location of reference genomes
#FASTA_HG38="/home/ascanolab/seq_pipelines/parpipe2_v0.2/genome/GRCh38.primary_assembly.genome.fa"
FASTA_CEDV="/home/ascanolab/Ascura/DirectRNAseq_Ascura/reference_genomes/GCF_000924595.1_ViralProj264921_genomic.fna"

# 3 - Specify location of results - output
DMR_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/dmr_output/individual_replicates"

# 4 - Perform differential methylation analysis between biological replicates and identical sample conditions
# Control_1 vs Control_2 OR STM_1 vs STM_2
modkit dmr pair \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed.gz" \
  -o "${DMR_OUTPUT}/cedv_bmt0.7_dmr_control_1_control_2.bed" \
  --ref "$FASTA_CEDV" \
  --base A \
  --threads 8 \
  --log-filepath "${DMR_OUTPUT}/cedv_dmr_control_replicates.log"

modkit dmr pair \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_1.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_2.bed.gz" \
  -o "${DMR_OUTPUT}/cedv_bmt0.7_dmr_stm_1_stm_2.bed" \
  --ref "$FASTA_CEDV" \
  --base A \
  --threads 8 \
  --log-filepath "${DMR_OUTPUT}/cedv_dmr_stm_replicates.log"  
  