# Load in packages
library(dplyr)
library(ggplot2)

# Title: TPM_Normalized_RDMR_m6A_Location.R
# Purpose: Evaluate % of TPM-normalized RDMR sites located in each CedV transcript, related to Table S4
# Author: Andrea Ascura
# Last Update: 03 Oct 2026

# 1 - Read in TPM-normalized RDMR values
stm_sig <- read.csv('./NanoCount/STM/stm_control_dmr_tpm_fdr_sig.csv')
spen_sig <- read.csv('./NanoCount/SPEN/spen_control_dmr_tpm_fdr_sig.csv')

# 2 - Define function to compute proportion of modification sites in given mRNA and proportion with 10% m6A reduction
calc_perc_rdmr <- function(gene_name, df_total) {
  gene_count <- df_total %>%
    dplyr::filter(gene_control == gene_name) # Number of sig RDMRs located in a gene
  
  gene_count_rt_1 <- gene_count %>%
    dplyr::filter(normalized_modA_diff <= -0.1) # Number of sig RDMRs located in a gene that meet reduction threshold
  
  gene_count_rt_05 <- gene_count %>%
    dplyr::filter(normalized_modA_diff <= -0.05) # Number of sig RDMRs located in a gene that meet reduction threshold
  
  mod_summary <- data.frame(
    gene = gene_name,
    num_sig_tpm_rdmr = nrow(gene_count),
    reduction_threshold = c(-0.1, -0.05),
    m6A_sites_reduction_threshold = c(nrow(gene_count_rt_1), nrow(gene_count_rt_05))
  )

  print(mod_summary)
  }

# 3 - Run function for STM data and combine 
stm_N <- calc_perc_rdmr('N', stm_sig)
stm_P <- calc_perc_rdmr('P', stm_sig)
stm_M <- calc_perc_rdmr('M', stm_sig)
stm_F <- calc_perc_rdmr('F', stm_sig)
stm_G <- calc_perc_rdmr('G', stm_sig)
stm_L <- calc_perc_rdmr('L', stm_sig)

stm_all <- rbind(stm_N, stm_P, stm_M, stm_F, stm_G, stm_L)

# 4 - Run function for SPEN data and combine
spen_N <- calc_perc_rdmr('N', spen_sig)
spen_P <- calc_perc_rdmr('P', spen_sig)
spen_M <- calc_perc_rdmr('M', spen_sig)
spen_F <- calc_perc_rdmr('F', spen_sig)
spen_G <- calc_perc_rdmr('G', spen_sig)
spen_L <- calc_perc_rdmr('L', spen_sig)

spen_all <- rbind(spen_N, spen_P, spen_M, spen_F, spen_G, spen_L)

# 5 - Write output file for future referencing - ENSURE % REDUCTION FILTER IS SPECIFIED IN FILENAME
write.csv(stm_all, './NanoCount/STM/stm_rdmr_locations_perc_rdmr_reduction_5_10.csv', row.names = FALSE)
write.csv(spen_all, './NanoCount/SPEN/spen_rdmr_locations_perc_rdmr_reduction_5_10.csv', row.names = FALSE)
