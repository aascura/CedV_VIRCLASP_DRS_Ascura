# Load in packages
library(dplyr)
library(ggplot2)
library(tidyr)

# Title: Total A vs m6A Count
# Purpose: Calculate m6A-to-A conversion rate per transcript (Horizontal m6A rate per transcript A's), related to Table S4
# Author: Andrea Ascura
# Date Updated: 03 Oct 2026


#### SPEN TPM RDMR Analysis ####
# 1 - Load in number of total A's per CedV transcript
# Reference Source: https://www.ncbi.nlm.nih.gov/nuccore/NC_025351.1
# Instructions: Number of A's per transcript were determining by count number of A occurrences in corresponding highlighted region of gene
N_A_count <- 635
P_A_count <- 927
M_A_count <- 576
F_A_count <- 743
G_A_count <- 740
L_A_count <- 2820

# 1 - Load in results for CedV Control and SPEN modA and total A counts
dmr_spen <- read.csv('./Nanopore_DRS/dmr_output/cedv_spen_control_rdmr_annotated_regions.csv')

tpm_dmr_spen <- read.csv('./Nanocount/SPEN/spen_control_dmr_tpm_fdr.csv')
tpm_dmr_spen_sig <- tpm_dmr_spen %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::rename(gene = gene_control)

# 2 - Determine # of A's and % of A's in given transcript that are m6A and those that are significant 
m6A_A <- function(dmr_df, transcript, gene_A_count, tpm_dmr_df) {
  
  condition_gene <- dmr_df %>%
    dplyr::filter(gene == transcript)
  m6A_transcript <- nrow(condition_gene) / gene_A_count * 100

  condition_gene_tpm <- tpm_dmr_df %>%
    dplyr::filter(gene == transcript)
  m6A_transcript_sig <- nrow(condition_gene_tpm) / nrow(condition_gene) * 100
  
  RDMR_min <- min(condition_gene_tpm$normalized_modA_diff, na.rm = FALSE) # Min value = more hypomethylated
  RDMR_max <- max(condition_gene_tpm$normalized_modA_diff, na.rm = FALSE) # Max value = less hypomethylated
  
  df <- data.frame(
    gene_name = transcript,
    num_m6A_A = nrow(condition_gene),
    num_A = gene_A_count,
    perc_m6A_A = m6A_transcript,
    num_sig_tpm_rdmr = nrow(condition_gene_tpm),
    perc_sig_tpm_rdmr = m6A_transcript_sig, # Percentage of total # of m6A
    max_rdmr = RDMR_min,
    min_rdmr = RDMR_max
  )
  return(df)
}

spen_N_m6A <- m6A_A(dmr_spen, 'N', N_A_count, tpm_dmr_spen_sig)
spen_P_m6A <- m6A_A(dmr_spen, 'P', P_A_count, tpm_dmr_spen_sig)
spen_M_m6A <- m6A_A(dmr_spen, 'M', M_A_count, tpm_dmr_spen_sig)
spen_F_m6A <- m6A_A(dmr_spen, 'F', F_A_count, tpm_dmr_spen_sig)
spen_G_m6A <- m6A_A(dmr_spen, 'G', G_A_count, tpm_dmr_spen_sig)
spen_L_m6A <- m6A_A(dmr_spen, 'L', L_A_count, tpm_dmr_spen_sig)

# 3 - Recombine individual transcript metrics
spen_m6A_all <- rbind(spen_N_m6A, spen_P_m6A, spen_M_m6A, spen_F_m6A, spen_G_m6A, spen_L_m6A)

# 4 - Export combined datafile 
write.csv(spen_m6A_all, './Nanopore_DRS/m6A_A_Transcript_Analysis/spen_m6A_A_transcripts.csv', row.names = FALSE)

#### STM2457 TPM RDMR Analysis ####

# 1 - Load in results for CedV Control and STM modA and total A counts
dmr_stm <- read.csv('./Nanopore_DRS/Apr2026_STM_outputs_only/dmr_output/cedv_stm_control_rdmr_annotated_regions.csv')

tpm_dmr_stm <- read.csv('./Nanocount/STM/stm_control_dmr_tpm_fdr.csv')
tpm_dmr_stm_sig <- tpm_dmr_stm %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::rename(gene = gene_control)

# 2 - Determine # of A's and % of A's in given transcript that are m6A and those that are significant 
# RUN SAME FUNCTION AS ABOVE IN STEP 2

stm_N_m6A <- m6A_A(dmr_stm, 'N', N_A_count, tpm_dmr_stm_sig)
stm_P_m6A <- m6A_A(dmr_stm, 'P', P_A_count, tpm_dmr_stm_sig)
stm_M_m6A <- m6A_A(dmr_stm, 'M', M_A_count, tpm_dmr_stm_sig)
stm_F_m6A <- m6A_A(dmr_stm, 'F', F_A_count, tpm_dmr_stm_sig)
stm_G_m6A <- m6A_A(dmr_stm, 'G', G_A_count, tpm_dmr_stm_sig)
stm_L_m6A <- m6A_A(dmr_stm, 'L', L_A_count, tpm_dmr_stm_sig)

# 3 - Recombine individual transcript metrics
stm_m6A_all <- rbind(stm_N_m6A, stm_P_m6A, stm_M_m6A, stm_F_m6A, stm_G_m6A, stm_L_m6A)

# 4 - Export combined datafile 
write.csv(stm_m6A_all, './Nanopore_DRS/m6A_A_Transcript_Analysis/stm_m6A_A_transcripts.csv', row.names = FALSE)

