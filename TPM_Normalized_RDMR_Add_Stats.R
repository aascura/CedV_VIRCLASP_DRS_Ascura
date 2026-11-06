# Load in packages
library(dplyr)
library(ggplot2)
library(readxl)

# Title: TPM_Normalized_RDMR_Add_Stats.R
# Purpose: Calculate differential modification rates (modA_ratio_spen - modA_ratio_control) and compute FDR for each differentially modified region (single nucleotide position), related to Figure 6D-E
# Author: Andrea Ascura
# Last Update: 12 June 2026

# sample_a = Control
# sample_b = SPEN

# 1 - Load in differential analysis (dmr_output) from Modkit_Differential_Methylation_Scoring.sh and remove unnecessary strings from numeric columns
# Modkit descriptions to dmr single base analysis column names: https://nanoporetech.github.io/modkit/intro_dmr.html

#### SPEN FDR Calculations ####
spen_tpm_df_raw <- read.csv('./NanoCount/SPEN/spen_control_tpm_normalized.csv')

spen_tpm_df <- spen_tpm_df_raw %>%
  dplyr::mutate(mod_diff = (modA_ratio_spen - modA_ratio_control)) %>%
  dplyr::mutate(normalized_modA_diff = (normalized_modA_spen - normalized_modA_control))

# 2 - Calculate p-value and adjusted P-value (FDR)
dmr_fdr_spen <- spen_tpm_df %>%
  # Calculate standard error
  mutate(se_diff = sqrt((normalized_modA_spen * (1 - normalized_modA_spen)) / (sample_b_total)) + ((normalized_modA_control * (1 - normalized_modA_control)) / (sample_a_total))) %>%
  # Remove cases where standard error is zero (proportions are 0 or 1)
  dplyr::filter(se_diff != 0) %>%
  # Calculate z-score
  mutate(z_score = normalized_modA_diff / se_diff) %>%
  # Calculate the two-tailed p-value from the z-score
  # This is 2 * P(Z > |z_score|)
  mutate(p_value = 2 * pnorm(-abs(z_score))) %>%
  # Calculate FDR from p-value using Benjamini-Hochberg correction
  mutate(fdr = p.adjust(p_value, method="BH"))

# 3 - Filter for only significant sites (FDR < 0.05)
dmr_fdr_spen_sig <- dmr_fdr_spen %>%
  dplyr::filter(fdr < 0.05)

# 4 - Export Normalized SPEN RDMR with stats
write.csv(dmr_fdr_spen, './NanoCount/SPEN/spen_control_dmr_tpm_fdr.csv', row.names = FALSE)
write.csv(dmr_fdr_spen_sig, './NanoCount/SPEN/spen_control_dmr_tpm_fdr_sig.csv', row.names = FALSE)

#### STM FDR Calculations ####

# 1 - Load in differential analysis (dmr_output) from Modkit_Differential_Methylation_Scoring.sh
stm_tpm_df_raw <- read.csv('./NanoCount/STM/stm_control_tpm_normalized.csv')

stm_tpm_df <- stm_tpm_df_raw %>%
  dplyr::mutate(mod_diff = (modA_ratio_stm - modA_ratio_control)) %>%
  dplyr::mutate(normalized_modA_diff = (normalized_modA_stm - normalized_modA_control))

# 2 - Calculate p-value and adjusted P-value (FDR)
dmr_fdr_stm <- stm_tpm_df %>%
  # Calculate standard error
  mutate(se_diff = sqrt((normalized_modA_stm * (1 - normalized_modA_stm)) / (sample_b_total)) + ((normalized_modA_control * (1 - normalized_modA_control)) / (sample_a_total))) %>%
  # Remove cases where standard error is zero (proportions are 0 or 1)
  dplyr::filter(se_diff != 0) %>%
  # Calculate z-score
  mutate(z_score = normalized_modA_diff / se_diff) %>%
  # Calculate the two-tailed p-value from the z-score
  # This is 2 * P(Z > |z_score|)
  mutate(p_value = 2 * pnorm(-abs(z_score))) %>%
  # Calculate FDR from p-value using Benjamini-Hochberg correction
  mutate(fdr = p.adjust(p_value, method="BH"))

# 3 - Filter for only significant sites (FDR < 0.05)
dmr_fdr_stm_sig <- dmr_fdr_stm %>%
  dplyr::filter(fdr < 0.05)

# 4 - Export Normalized STM RDMR with stats
write.csv(dmr_fdr_stm, './NanoCount/STM/stm_control_dmr_tpm_fdr.csv', row.names = FALSE)
write.csv(dmr_fdr_stm_sig, './NanoCount/STM/stm_control_dmr_tpm_fdr_sig.csv', row.names = FALSE)
