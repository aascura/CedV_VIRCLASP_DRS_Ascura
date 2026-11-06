# Load in packages
library(dplyr)
library(ggplot2)

# Title: Absolute_ModRate_Metagene.R
# Purpose: Visualize absolute modification counts (no TPM_normalization) across CedV genome in Control and SPEN KD/METTL3i mRNA transcripts as a metagene-like plot, related to Figure S6H
# Author: Andrea Ascura
# Last Update: 22 June 2026


#### SPEN Metagene ####
# 1 - Load in results for CedV Control and SPEN modA and total A counts
dmr_spen <- read.csv('./NanoCount/SPEN/spen_control_dmr_tpm_fdr.csv')

# 2 - Reformat dataframe to enable plotting of SPEN and Control on same graph 
cedv_control_1_2 <- dmr_spen %>%
  dplyr::select(start, sample_a_counts) %>%
  dplyr::rename('sample_counts' = 'sample_a_counts')
cedv_control_1_2$sample <- 'control'

cedv_spen <- dmr_spen %>%
  dplyr::select(start, sample_b_counts) %>%
  dplyr::rename('sample_counts' = 'sample_b_counts') 
cedv_spen$sample <- 'spen'

# 3 - Bind categorized datasets together
cedv_spen_bind <- rbind(cedv_control_1_2,cedv_spen)

# 4 - Grouped Plots: Visualize modification ratios for both Control and SPEN sample on entire CedV genome 
# 5 - Define CedV genomic boundaries (gene start + end) (nucleotide positions) to subset data into regions and annotate full DMR
# https://www.ncbi.nlm.nih.gov/nuccore/JQ001776
# N: 56...2010
# N CDS: 144...1676
# P: 2014...4517
# P CDS: 2112...4325
# M: 4521...6125
# M CDS: 4635...5717
# F: 6129...8166
# F CDS: 6405...8078
# G: 8170...10275
# G CDS: 8268...10136
# L: 10279...18140
# L CDS: 10572...18077

# 6 - Plot modification ratios for Control and SPEN sample with annotated gene and CDS boundaries
p1 <- ggplot(cedv_spen_bind, aes(x = start, y = log10(sample_counts))) +
  #geom_bar(stat = 'identity', aes(color = sample)) + 
  geom_smooth(aes(color = sample), se = FALSE, span = 0.20) +
  #ylim(5,10) +
  #geom_vline(xintercept = merged_smooth_line_high$x, linetype = 'dashed', color = 'purple') +
  #geom_vline(xintercept = merged_smooth_line_low$x, linetype = 'dashed', color = 'green') +
  
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  geom_vline(xintercept = 56, linetype = 'dashed', color = 'grey25') + # N gene start
  geom_vline(xintercept = 2010, linetype = 'dashed', color = 'grey25') + # N gene end
  geom_vline(xintercept = 144, linetype = 'dashed', color = 'lightgrey') + # N CDS start
  geom_vline(xintercept = 1676, linetype = 'dashed', color = 'lightgrey') + # N CDS end
  
  geom_vline(xintercept = 2014, linetype = 'dashed', color = 'grey25') + # P gene start
  geom_vline(xintercept = 4517, linetype = 'dashed', color = 'grey25') + # P gene end
  geom_vline(xintercept = 2112, linetype = 'dashed', color = 'lightgrey') + # P CDS start
  geom_vline(xintercept = 4325, linetype = 'dashed', color = 'lightgrey') + # P CDS end
  
  geom_vline(xintercept = 4521, linetype = 'dashed', color = 'grey25') + # M gene start
  geom_vline(xintercept = 6125, linetype = 'dashed', color = 'grey25') + # M gene end
  geom_vline(xintercept = 4635, linetype = 'dashed', color = 'lightgrey') + # M CDS start
  geom_vline(xintercept = 5717, linetype = 'dashed', color = 'lightgrey') + # M CDS end
  
  geom_vline(xintercept = 6129, linetype = 'dashed', color = 'grey25') + # F gene start
  geom_vline(xintercept = 8166, linetype = 'dashed', color = 'grey25') + # F gene end
  geom_vline(xintercept = 6405, linetype = 'dashed', color = 'lightgrey') + # F CDS start
  geom_vline(xintercept = 8078, linetype = 'dashed', color = 'lightgrey') + # F CDS end
  
  geom_vline(xintercept = 8170, linetype = 'dashed', color = 'grey25') + # G gene start
  geom_vline(xintercept = 10275, linetype = 'dashed', color = 'grey25') + # G gene end
  geom_vline(xintercept = 8268, linetype = 'dashed', color = 'lightgrey') + # G CDS start
  geom_vline(xintercept = 10136, linetype = 'dashed', color = 'lightgrey') + # G CDS end
  
  geom_vline(xintercept = 10279, linetype = 'dashed', color = 'grey25') + # L gene start
  geom_vline(xintercept = 18140, linetype = 'dashed', color = 'grey25') + # L gene end
  geom_vline(xintercept = 10572, linetype = 'dashed', color = 'lightgrey') + # L CDS start
  geom_vline(xintercept = 18077, linetype = 'dashed', color = 'lightgrey') + # L CDS end
  xlab('Nucleotide Position') + 
  ylab('Absolute Modification (log10)') +
  theme_classic()
p1 + labs(title = 'siSPEN CedV mRNA')

#### STM Metagene ####
# 1 - Load in results for CedV Control and SPEN modA and total A counts
dmr_stm <- read.csv('./NanoCount/STM/stm_control_dmr_tpm_fdr.csv')

# 2 - Reformat dataframe to enable plotting of SPEN and Control on same graph 
cedv_control_3_4 <- dmr_stm %>%
  dplyr::select(start, sample_a_counts) %>%
  dplyr::rename('sample_counts' = 'sample_a_counts')
cedv_control_3_4$sample <- 'control'

cedv_stm <- dmr_stm %>%
  dplyr::select(start, sample_b_counts) %>%
  dplyr::rename('sample_counts' = 'sample_b_counts')
cedv_stm$sample <- 'stm'

# 3 - Bind categorized datasets together
cedv_stm_bind <- rbind(cedv_control_3_4,cedv_stm)

# 4 - Grouped Plots: Visualize modification ratios for both Control and SPEN sample on entire CedV genome 
# 5 - Define CedV genomic boundaries (gene start + end) (nucleotide positions) to subset data into regions and annotate full DMR
# https://www.ncbi.nlm.nih.gov/nuccore/JQ001776
# N: 56...2010
# N CDS: 144...1676
# P: 2014...4517
# P CDS: 2112...4325
# M: 4521...6125
# M CDS: 4635...5717
# F: 6129...8166
# F CDS: 6405...8078
# G: 8170...10275
# G CDS: 8268...10136
# L: 10279...18140
# L CDS: 10572...18077

# 6 - Plot modification ratios for Control and STM sample with annotated gene and CDS boundaries
p2 <- ggplot(cedv_stm_bind, aes(x = start, y = log10(sample_counts))) +
  #geom_bar(stat = 'identity', aes(color = sample)) + 
  geom_smooth(aes(color = sample), se = FALSE, span = 0.2, method = 'loess') +
  #ylim(5,10) +
  #geom_vline(xintercept = merged_smooth_line_high$x, linetype = 'dashed', color = 'purple') +
  #geom_vline(xintercept = merged_smooth_line_low$x, linetype = 'dashed', color = 'green') +
  
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  geom_vline(xintercept = 56, linetype = 'dashed', color = 'grey25') + # N gene start
  geom_vline(xintercept = 2010, linetype = 'dashed', color = 'grey25') + # N gene end
  geom_vline(xintercept = 144, linetype = 'dashed', color = 'lightgrey') + # N CDS start
  geom_vline(xintercept = 1676, linetype = 'dashed', color = 'lightgrey') + # N CDS end
  
  geom_vline(xintercept = 2014, linetype = 'dashed', color = 'grey25') + # P gene start
  geom_vline(xintercept = 4517, linetype = 'dashed', color = 'grey25') + # P gene end
  geom_vline(xintercept = 2112, linetype = 'dashed', color = 'lightgrey') + # P CDS start
  geom_vline(xintercept = 4325, linetype = 'dashed', color = 'lightgrey') + # P CDS end
  
  geom_vline(xintercept = 4521, linetype = 'dashed', color = 'grey25') + # M gene start
  geom_vline(xintercept = 6125, linetype = 'dashed', color = 'grey25') + # M gene end
  geom_vline(xintercept = 4635, linetype = 'dashed', color = 'lightgrey') + # M CDS start
  geom_vline(xintercept = 5717, linetype = 'dashed', color = 'lightgrey') + # M CDS end
  
  geom_vline(xintercept = 6129, linetype = 'dashed', color = 'grey25') + # F gene start
  geom_vline(xintercept = 8166, linetype = 'dashed', color = 'grey25') + # F gene end
  geom_vline(xintercept = 6405, linetype = 'dashed', color = 'lightgrey') + # F CDS start
  geom_vline(xintercept = 8078, linetype = 'dashed', color = 'lightgrey') + # F CDS end
  
  geom_vline(xintercept = 8170, linetype = 'dashed', color = 'grey25') + # G gene start
  geom_vline(xintercept = 10275, linetype = 'dashed', color = 'grey25') + # G gene end
  geom_vline(xintercept = 8268, linetype = 'dashed', color = 'lightgrey') + # G CDS start
  geom_vline(xintercept = 10136, linetype = 'dashed', color = 'lightgrey') + # G CDS end
  
  geom_vline(xintercept = 10279, linetype = 'dashed', color = 'grey25') + # L gene start
  geom_vline(xintercept = 18140, linetype = 'dashed', color = 'grey25') + # L gene end
  geom_vline(xintercept = 10572, linetype = 'dashed', color = 'lightgrey') + # L CDS start
  geom_vline(xintercept = 18077, linetype = 'dashed', color = 'lightgrey') + # L CDS end
  xlab('Nucleotide Position') + 
  ylab('Absolute Modification (log10)') +
  theme_classic()
p2 + labs(title = 'METTL3i CedV mRNA')

# Combine SPEN and STM absolute modification ratio plots
# 1 - Combine SPEN and STM datasets
cedv_all_bind <- rbind(cedv_spen_bind,cedv_stm_bind) %>%
  dplyr::filter(sample_counts > 1)

# 2 - Plot combined abolsute modification ratios for SPEN and STM2457 datasets
p3 <- ggplot(cedv_all_bind, aes(x = start, y = log10(sample_counts))) +
  #geom_bar(stat = 'identity', aes(color = sample)) + 
  geom_smooth(aes(color = sample), se = FALSE, span = 0.15, method = 'loess') +
  #geom_vline(xintercept = merged_smooth_line_high$x, linetype = 'dashed', color = 'purple') +
  #geom_vline(xintercept = merged_smooth_line_low$x, linetype = 'dashed', color = 'green') +
  
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  geom_vline(xintercept = 56, linetype = 'dashed', color = 'grey25') + # N gene start
  geom_vline(xintercept = 2010, linetype = 'dashed', color = 'grey25') + # N gene end
  geom_vline(xintercept = 144, linetype = 'dashed', color = 'lightgrey') + # N CDS start
  geom_vline(xintercept = 1676, linetype = 'dashed', color = 'lightgrey') + # N CDS end
  
  geom_vline(xintercept = 2014, linetype = 'dashed', color = 'grey25') + # P gene start
  geom_vline(xintercept = 4517, linetype = 'dashed', color = 'grey25') + # P gene end
  geom_vline(xintercept = 2112, linetype = 'dashed', color = 'lightgrey') + # P CDS start
  geom_vline(xintercept = 4325, linetype = 'dashed', color = 'lightgrey') + # P CDS end
  
  geom_vline(xintercept = 4521, linetype = 'dashed', color = 'grey25') + # M gene start
  geom_vline(xintercept = 6125, linetype = 'dashed', color = 'grey25') + # M gene end
  geom_vline(xintercept = 4635, linetype = 'dashed', color = 'lightgrey') + # M CDS start
  geom_vline(xintercept = 5717, linetype = 'dashed', color = 'lightgrey') + # M CDS end
  
  geom_vline(xintercept = 6129, linetype = 'dashed', color = 'grey25') + # F gene start
  geom_vline(xintercept = 8166, linetype = 'dashed', color = 'grey25') + # F gene end
  geom_vline(xintercept = 6405, linetype = 'dashed', color = 'lightgrey') + # F CDS start
  geom_vline(xintercept = 8078, linetype = 'dashed', color = 'lightgrey') + # F CDS end
  
  geom_vline(xintercept = 8170, linetype = 'dashed', color = 'grey25') + # G gene start
  geom_vline(xintercept = 10275, linetype = 'dashed', color = 'grey25') + # G gene end
  geom_vline(xintercept = 8268, linetype = 'dashed', color = 'lightgrey') + # G CDS start
  geom_vline(xintercept = 10136, linetype = 'dashed', color = 'lightgrey') + # G CDS end
  
  geom_vline(xintercept = 10279, linetype = 'dashed', color = 'grey25') + # L gene start
  geom_vline(xintercept = 18140, linetype = 'dashed', color = 'grey25') + # L gene end
  geom_vline(xintercept = 10572, linetype = 'dashed', color = 'lightgrey') + # L CDS start
  geom_vline(xintercept = 18077, linetype = 'dashed', color = 'lightgrey') + # L CDS end
  xlab('Nucleotide Position') + 
  ylab('Absolute Modification Count (log10)') +
  theme_classic()
  #theme(legend.position = 'none')
p3 + labs(title = 'METTL3i & siSPEN CedV mRNA')

# 3 - Export plot and specify dimensions
ggsave(p3, filename = "./NanoCount/plot_output/cedv_absolute_count_metagene_spen_stm.pdf", 
       width = 10, height = 3)

