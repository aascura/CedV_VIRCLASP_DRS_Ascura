# Load in packages
library(dplyr)
library(data.table)
library(stringr)
library(corrplot)
library(ggplot2)
library(ggpubr)

# Title: ModRate_Correlation_Plots.R
# Purpose: Assess and plot correlation for SPEN KD/METTL3 inhibition and Control Replicates for RDMR from Modkit
# Author: Andrea Ascura
# Last Update: 26 June 2026

#### SPEN Correlation Analysis ####
# 1 - Load in CedV DMRs for SPEN KD and Control samples (Control_1 vs SPEN_1)
dmr_spen_control_1_raw <- fread('./Nanopore_DRS/dmr_output/individual_replicates/cedv_bmt0.7_dmr_control_spen_rep1.csv')
dmr_spen_control_2_raw <- fread('./Nanopore_DRS/dmr_output/individual_replicates/cedv_bmt0.7_dmr_control_spen_rep2.csv')

# 2A - Format datasets and select necessary columns for correlation plots for Replicate 1
dmr_spen_control_1 <- dmr_spen_control_1_raw %>%
  dplyr::select('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14') %>%
  dplyr::rename(chromosome = V1,
                start = V2,
                end = V3,
                name = V4,
                score = V5,
                strand = V6,
                sample_a_counts = V7,
                sample_a_total = V8,
                sample_b_counts = V9,
                sample_b_total = V10,
                sample_a_percents = V11,
                sample_b_percents = V12,
                sample_a_fraction_modified = V13,
                sample_b_fraction_modified = V14)

dmr_spen_control_1 <- dmr_spen_control_1 %>%
  mutate(sample_a_percents = str_remove(sample_a_percents,'a:')) %>%
  mutate(sample_a_counts = str_remove(sample_a_counts,'a:')) %>%
  mutate(sample_b_percents = str_remove(sample_b_percents,'a:')) %>%
  mutate(sample_b_counts = str_remove(sample_b_counts,'a:'))

dmr_spen_control_1$sample_a_percents <- as.numeric(as.character(dmr_spen_control_1$sample_a_percents))
dmr_spen_control_1$sample_b_percents <- as.numeric(as.character(dmr_spen_control_1$sample_b_percents))
dmr_spen_control_1$sample_a_counts <- as.numeric(as.character(dmr_spen_control_1$sample_a_counts))
dmr_spen_control_1$sample_b_counts <- as.numeric(as.character(dmr_spen_control_1$sample_b_counts))

# 2B - Format datasets and select necessary columns for correlation plots for Replicate 2
dmr_spen_control_2 <- dmr_spen_control_2_raw %>%
  dplyr::select('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14') %>%
  dplyr::rename(chromosome = V1,
                start = V2,
                end = V3,
                name = V4,
                score = V5,
                strand = V6,
                sample_a_counts = V7,
                sample_a_total = V8,
                sample_b_counts = V9,
                sample_b_total = V10,
                sample_a_percents = V11,
                sample_b_percents = V12,
                sample_a_fraction_modified = V13,
                sample_b_fraction_modified = V14)

dmr_spen_control_2 <- dmr_spen_control_2 %>%
  mutate(sample_a_percents = str_remove(sample_a_percents,'a:')) %>%
  mutate(sample_a_counts = str_remove(sample_a_counts,'a:')) %>%
  mutate(sample_b_percents = str_remove(sample_b_percents,'a:')) %>%
  mutate(sample_b_counts = str_remove(sample_b_counts,'a:'))

dmr_spen_control_2$sample_a_percents <- as.numeric(as.character(dmr_spen_control_2$sample_a_percents))
dmr_spen_control_2$sample_b_percents <- as.numeric(as.character(dmr_spen_control_2$sample_b_percents))
dmr_spen_control_2$sample_a_counts <- as.numeric(as.character(dmr_spen_control_2$sample_a_counts))
dmr_spen_control_2$sample_b_counts <- as.numeric(as.character(dmr_spen_control_2$sample_b_counts))

# 3A - Calculate modification ratio difference (mod_diff): modA_ratio_spen - modA_ratio_control for Replicate 1
dmr_spen_control_1_mod_diff <- dmr_spen_control_1 %>%
  mutate(modA_ratio_control = sample_a_counts/sample_a_total) %>% # p1; sample_a = control
  mutate(modA_ratio_spen = sample_b_counts/sample_b_total) %>% # p2; sample_b = spen
  mutate(mod_diff = modA_ratio_spen - modA_ratio_control) %>%
  dplyr::select(start, end, sample_a_counts, sample_a_total, sample_b_counts, sample_b_total, modA_ratio_control, modA_ratio_spen, mod_diff)

# 3B - Calculate modification ratio difference (mod_diff): modA_ratio_spen - modA_ratio_control for Replicate 2
dmr_spen_control_2_mod_diff <- dmr_spen_control_2 %>%
  mutate(modA_ratio_control = sample_a_counts/sample_a_total) %>% # p1; sample_a = control
  mutate(modA_ratio_spen = sample_b_counts/sample_b_total) %>% # p2; sample_b = spen
  mutate(mod_diff = modA_ratio_spen - modA_ratio_control) %>%
  dplyr::select(start, end, sample_a_counts, sample_a_total, sample_b_counts, sample_b_total, modA_ratio_control, modA_ratio_spen, mod_diff)

# 4 - Export formatted DMR table for future reference and analysis
write.csv(dmr_spen_control_1_mod_diff,
          './Nanopore_DRS/dmr_output/individual_replicates/dmr_spen_control_1_mod_diff.csv',
          row.names = FALSE)
write.csv(dmr_spen_control_2_mod_diff,
          './Nanopore_DRS/dmr_output/individual_replicates/dmr_spen_control_2_mod_diff.csv',
          row.names = FALSE)

# 5 - Import formatted DMR df for replicate comparisons
# Correlation comparisons and resulting value ( = # represents correlation coefficient):
# Control_1 modA_ratio_control vs Control_2 modA_ratio_control = 0.96
# SPEN_1 modA_ratio_spen vs SPEN_2 modA_ratio_spen = 0.83

spen_control_1 <- fread('./Nanopore_DRS/dmr_output/individual_replicates/dmr_spen_control_1_mod_diff.csv') %>%
  dplyr::filter(sample_a_total >= 50 & sample_b_total >= 50) %>%
  dplyr::select(start,modA_ratio_spen) # Replace with _control or _spen
spen_control_2 <- fread('./Nanopore_DRS/dmr_output/individual_replicates/dmr_spen_control_2_mod_diff.csv') %>%
  dplyr::filter(sample_a_total >= 50 & sample_b_total >= 50) %>%
  dplyr::select(start,modA_ratio_spen) # Replace with _control or _spen

# 6 - Join replicate columns by common start positions
combined_reps <- dplyr::inner_join(spen_control_2,spen_control_1,by = join_by(start)) %>%
  dplyr::select(modA_ratio_spen.x, modA_ratio_spen.y) # Replace with control or spen

# 7 - Evaluate and plot correlation using corrplot
combined_reps_corrplt <- ggscatter(combined_reps, x = 'modA_ratio_spen.x', y = 'modA_ratio_spen.y',
                                   add = "reg.line", 
                                   add.params = list(color = '#1D4497FF'),
                                   palette = 'npg',
                                   size = 2,
                                   xlab = 'ModA Ratio (SPEN KD Rep. 1)', ylab = 'ModA Ratio (SPEN KD Rep. 2)',
                                   font.x = c(12,'bold'), font.y = c(12,'bold')
)+ stat_cor(method = "pearson")
ggpar(combined_reps_corrplt)

# 8 - Export correlation plot image
ggsave(combined_reps_corrplt, filename = "./Nanopore_DRS/plot_outputs/corrplt_modA_ratio_spen.pdf", width = 4, height = 4) # Replace with control or spen

#### STM Correlation Analysis ####
# 1 - Load in CedV DMRs for SPEN KD and Control samples (Control_1 vs STM_1)
dmr_control_raw <- fread('./Individual_Replicates_STM/cedv_bmt0.7_dmr_control_1_control_2.csv')
dmr_stm_raw <- fread('./Individual_Replicates_STM/cedv_bmt0.7_dmr_stm_1_stm_2.csv')

# 2 - Format datasets and select necessary columns for correlation plots
dmr_format <-function(dmr_df){
  dmr_dataset <- dmr_df %>%
  dplyr::select('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14') %>%
  dplyr::rename(chromosome = V1,
                start = V2,
                end = V3,
                name = V4,
                score = V5,
                strand = V6,
                sample_a_counts = V7,
                sample_a_total = V8,
                sample_b_counts = V9,
                sample_b_total = V10,
                sample_a_percents = V11,
                sample_b_percents = V12,
                sample_a_fraction_modified = V13,
                sample_b_fraction_modified = V14)
  dmr_dataset_2 <- dmr_dataset %>%
    mutate(sample_a_percents = str_remove(sample_a_percents,'a:')) %>%
    mutate(sample_a_counts = str_remove(sample_a_counts,'a:')) %>%
    mutate(sample_b_percents = str_remove(sample_b_percents,'a:')) %>%
    mutate(sample_b_counts = str_remove(sample_b_counts,'a:'))
  
  dmr_dataset_2$sample_a_percents <- as.numeric(as.character(dmr_dataset_2$sample_a_percents))
  dmr_dataset_2$sample_b_percents <- as.numeric(as.character(dmr_dataset_2$sample_b_percents))
  dmr_dataset_2$sample_a_counts <- as.numeric(as.character(dmr_dataset_2$sample_a_counts))
  dmr_dataset_2$sample_b_counts <- as.numeric(as.character(dmr_dataset_2$sample_b_counts))
  
  return(dmr_dataset_2)
}
dmr_stm <- dmr_format(dmr_stm_raw)

dmr_stm_filt <- dmr_stm %>%
  dplyr::filter(sample_a_total >= 50 & sample_b_total >= 50) %>%
  dplyr::filter(sample_a_fraction_modified != 1 & sample_b_fraction_modified != 1) %>%
  dplyr::select(sample_a_fraction_modified, sample_b_fraction_modified)

# 3 - Plot correlations
stm_corrplt <- ggscatter(dmr_stm_filt, x = 'sample_a_fraction_modified', y = 'sample_b_fraction_modified',
                                   add = "reg.line", 
                                   add.params = list(color = '#1D4497FF'),
                                   palette = 'npg',
                                   size = 2,
                                   xlab = 'ModA Ratio (STM2457 Rep. 1)', ylab = 'ModA Ratio (STM2457 Rep. 2)',
                                   font.x = c(12,'bold'), font.y = c(12,'bold')
)+ stat_cor(method = "pearson") + xlim(0.0,0.04) + ylim(0.0,0.04)
ggpar(stm_corrplt)

# 4 - Export correlation plot image
ggsave(stm_corrplt, filename = "./Individual_Replicates_STM/plot_outputs/corrplt_modA_ratio_stm.pdf", width = 4, height = 4) # Replace with control or spen
