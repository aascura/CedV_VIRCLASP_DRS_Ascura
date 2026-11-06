#!/bin/bash

# Title: Oxford Dorado DRS Analysis - Alignment
# Purpose:Use Dorado to align basecalled Oxford Nanopore DRS data to CedV transcriptome
# Author: Andrea Ascura
# Created: 02Sep2025
# Last Update: 04Nov2026

# Specify location of CedV transcriptome and Hg38 reference used for reads from STM2457 and SPEN experiments
FASTA_CEDV="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/cedv_transcriptome.fna"
GTF_CEDV="/home/ascanolab/Ascura/DirectRNAseq_Ascura/reference_genomes/genomic.gtf"
FASTA_HG38="/home/ascanolab/seq_pipelines/parpipe2_v0.2/genome/GRCh38.primary_assembly.genome.fa"
GTF_HG38="/home/ascanolab/seq_pipelines/parpipe2_v0.2/genome/gencode.v45.primary_assembly.basic.annotation.gtf"

# SPEN Minimap2 Alignment:
# Specify location of basecalled BAM files for SPEN experiments in bam_pass folder (input)
CALLED_OUTPUT_SPEN="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Oct2025_SPEN/m6A_DRACH/basecalled_output/bam_pass"

# Specify file names of input files
BAM_CONTROL_1="$CALLED_OUTPUT_SPEN/called_control_1.bam"
BAM_CONTROL_2="$CALLED_OUTPUT_SPEN/called_control_2.bam"
BAM_SPEN_1="$CALLED_OUTPUT_SPEN/called_spen_1.bam"
BAM_SPEN_2="$CALLED_OUTPUT_SPEN/called_spen_2.bam"

# Specify location of aligned output files (output)
ALIGNED_OUTPUT_SPEN="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/SPEN/aligned_output"

# Run minimap2 aligner for basecalled reads against CedV reference for basecalled Control_1, SPEN_1, Control_2 and SPEN_2
dorado aligner "$FASTA_CEDV" "$BAM_CONTROL_1" --output-dir "$ALIGNED_OUTPUT_SPEN" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""
dorado aligner "$FASTA_CEDV" "$BAM_CONTROL_2" --output-dir "$ALIGNED_OUTPUT_SPEN" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""
dorado aligner "$FASTA_CEDV" "$BAM_SPEN_1" --output-dir "$ALIGNED_OUTPUT_SPEN" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""
dorado aligner "$FASTA_CEDV" "$BAM_SPEN_2" --output-dir "$ALIGNED_OUTPUT_SPEN" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""

# Run minimap2 aligner for basecalled reads against Hg38 reference genome for basecalled Control_1, SPEN_1, Control_2 and SPEN_2
dorado aligner "$FASTA_HG38" "$BAM_CONTROL_1" --output-dir "$ALIGNED_OUTPUT_SPEN/hg38_aligned_control_1" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""
dorado aligner "$FASTA_HG38" "$BAM_CONTROL_2" --output-dir "$ALIGNED_OUTPUT_SPEN/hg38_aligned_control_2" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""
dorado aligner "$FASTA_HG38" "$BAM_SPEN_1" --output-dir "$ALIGNED_OUTPUT_SPEN/hg38_aligned_spen_1" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""
dorado aligner "$FASTA_HG38" "$BAM_SPEN_2" --output-dir "$ALIGNED_OUTPUT_SPEN/hg38_aligned_spen_2" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""

# STM2457 Minimap2 Alignment:
# Specify location of basecalled BAM files for STM2457 experiments in bam_pass folder (input)
CALLED_OUTPUT_STM="/home/ascanolab/Ascura/DirectRNAseq_Ascura/Apr2026_STM/basecalled_output/bam_pass"

# Specify file names of input files
BAM_CONTROL_3="$CALLED_OUTPUT_STM/calls_control_1.bam"
BAM_CONTROL_4="$CALLED_OUTPUT_STM/calls_control_2.bam"
BAM_STM_1="$CALLED_OUTPUT_STM/calls_stm_1.bam"
BAM_STM_2="$CALLED_OUTPUT_STM/calls_stm_2.bam"

# Specify location of aligned output files (output)
ALIGNED_OUTPUT_STM="/home/ascanolab/Ascura/DirectRNAseq_Ascura/June2026_Normalization/STM/aligned_output"

# Run minimap2 aligner for basecalled reads against CedV reference for basecalled Control_1, STM_1, Control_2 and STM_2
dorado aligner "$FASTA_CEDV" "$BAM_CONTROL_3" --output-dir "$ALIGNED_OUTPUT_STM" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""
dorado aligner "$FASTA_CEDV" "$BAM_CONTROL_4" --output-dir "$ALIGNED_OUTPUT_STM" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""
dorado aligner "$FASTA_CEDV" "$BAM_STM_1" --output-dir "$ALIGNED_OUTPUT_STM" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""
dorado aligner "$FASTA_CEDV" "$BAM_STM_2" --output-dir "$ALIGNED_OUTPUT_STM" --emit-summary --mm2-opts "--junc-bed "$GTF_CEDV""

# Run minimap2 aligner for basecalled reads against Hg38 reference genome for basecalled Control_1, STM_1, Control_2 and STM_2
dorado aligner "$FASTA_HG38" "$BAM_CONTROL_3" --output-dir "$ALIGNED_OUTPUT_STM/hg38_aligned_control_3" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""
dorado aligner "$FASTA_HG38" "$BAM_CONTROL_4" --output-dir "$ALIGNED_OUTPUT_STM/hg38_aligned_control_4" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""
dorado aligner "$FASTA_HG38" "$BAM_STM_1" --output-dir "$ALIGNED_OUTPUT_STM/hg38_aligned_stm_1" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""
dorado aligner "$FASTA_HG38" "$BAM_STM_2" --output-dir "$ALIGNED_OUTPUT_STM/hg38_aligned_stm_2" --emit-summary --mm2-opts "--junc-bed "$GTF_HG38""

