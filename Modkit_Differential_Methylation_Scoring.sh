#!/bin/bash

# Title: Oxford Modkit DRS Analysis - Differential Methylation Scoring
# Purpose: Perform differential methylation analysis (WT vs SPEN KD; WT vs METTL3i aka STM2457) at single base positions with basecalled A (aka without region segmentation)
# Reference: https://nanoporetech.github.io/modkit/intro_dmr.html
# Author: Andrea Ascura
# Created: 05Sep2025
# Last Update: 04Nov2026

# SPEN Differential Methylation Scoring
# 1 - Specify location of aligned (indexed and sorted) modBAMs
PILEUP_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/pileup_output_V2"

# 2 - Specify location of reference genomes
#FASTA_HG38="/home/ascanolab/seq_pipelines/parpipe2_v0.2/genome/GRCh38.primary_assembly.genome.fa"
FASTA_CEDV="/home/ascanolab/Ascura/DirectRNAseq_Ascura/reference_genomes/GCF_000924595.1_ViralProj264921_genomic.fna"

# 3 - Specify location of results - output
DMR_OUTPUT="/media/ascanolab/Sequence_nvme1n1/Andrea_Ascura/Direct_RNAseq/dmr_output"

# 4 - Perform differential methylation analysis
modkit dmr pair \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed.gz" \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_spen_1.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_spen_2.bed.gz" \
  -o "${DMR_OUTPUT}/cedv_bmt0.7_dmr_control_spen.bed" \
  --ref "$FASTA_CEDV" \
  --base A \
  --threads 8 \
  --log-filepath "${DMR_OUTPUT}/cedv_dmr_control_spen.log"

# STM2457 Differential Methylation Scoring
# 1 - Specify location of aligned (indexed and sorted) modBAMs
PILEUP_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/pileup_output"

# 2 - Specify location of reference genomes
#FASTA_HG38="/home/ascanolab/seq_pipelines/parpipe2_v0.2/genome/GRCh38.primary_assembly.genome.fa"
FASTA_CEDV="/home/ascanolab/Ascura/DirectRNAseq_Ascura/reference_genomes/GCF_000924595.1_ViralProj264921_genomic.fna"

# 3 - Specify location of results - output
DMR_OUTPUT="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/dmr_output"

# 4 - Perform differential methylation analysis
modkit dmr pair \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_1.bed.gz" \
  -a "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_control_2.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_1.bed.gz" \
  -b "${PILEUP_OUTPUT}/pileup_bmt0.7_cedv_stm_2.bed.gz" \
  -o "${DMR_OUTPUT}/cedv_bmt0.7_dmr_control_stm.bed" \
  --ref "$FASTA_CEDV" \
  --base A \
  --threads 8 \
  --log-filepath "${DMR_OUTPUT}/cedv_dmr_control_stm.log"


