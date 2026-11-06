# Load in packages
library(dplyr)
library(ggplot2)

# Title: TPM_Normalized_RDMR_ModRate_Differences_Individual_Transcripts.R
# Purpose: Annotate positions and visualize relative differentially methylated regions (RDMR) using modification rate differences on individual transcripts, related to Figure 6D-E
# Author: Andrea Ascura
# Last Update: 01 Oct 2026

#### SPEN ####
# 1 - Load in CedV single nt positions with MRDs and FDRs
dmr_spen_sig_raw <- read.csv('./NanoCount/SPEN/spen_control_dmr_tpm_fdr_sig.csv')

# 2 - Assign color value to significant MRD directions (Positive = Red; Negative = Blue)
dmr_spen_sig <- dmr_spen_sig_raw
dmr_spen_sig$mod_direction <- ''
dmr_spen_sig$mod_direction <- ifelse(dmr_spen_sig$normalized_modA_diff > 0, 'pos',
                                     ifelse(dmr_spen_sig$normalized_modA_diff < 0, 'neg', 'unchanged'))

#### N mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_spen_sig_transcript <- dmr_spen_sig %>%
  dplyr::filter(gene_control == 'N')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
spen_rdmr_transcript <- ggplot(dmr_spen_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  geom_vline(xintercept = 56, linetype = 'dashed', color = 'grey25') + # N gene start
  geom_vline(xintercept = 2010, linetype = 'dashed', color = 'grey25') + # N gene end
  geom_vline(xintercept = 144, linetype = 'dashed', color = 'lightgrey') + # N CDS start
  geom_vline(xintercept = 1676, linetype = 'dashed', color = 'lightgrey') + # N CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (SPEN KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
spen_rdmr_transcript + labs(title = 'SPEN KD CedV N mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(spen_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/SPEN/cedv_spen_rdmr_curve_N.pdf", #EDIT FILENAME!
       width = 1.45, height = 2.5)

#### P mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_spen_sig_transcript <- dmr_spen_sig %>%
  dplyr::filter(gene_control == 'P')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
spen_rdmr_transcript <- ggplot(dmr_spen_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 2014, linetype = 'dashed', color = 'grey25') + # P gene start
  geom_vline(xintercept = 4517, linetype = 'dashed', color = 'grey25') + # P gene end
  geom_vline(xintercept = 2112, linetype = 'dashed', color = 'lightgrey') + # P CDS start
  geom_vline(xintercept = 4325, linetype = 'dashed', color = 'lightgrey') + # P CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (SPEN KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
spen_rdmr_transcript + labs(title = 'SPEN KD CedV P mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(spen_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/SPEN/cedv_spen_rdmr_curve_P.pdf", #EDIT FILENAME!
       width = 2.18, height = 2.5)

#### M mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_spen_sig_transcript <- dmr_spen_sig %>%
  dplyr::filter(gene_control == 'M')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
spen_rdmr_transcript <- ggplot(dmr_spen_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 4521, linetype = 'dashed', color = 'grey25') + # M gene start
  geom_vline(xintercept = 6125, linetype = 'dashed', color = 'grey25') + # M gene end
  geom_vline(xintercept = 4635, linetype = 'dashed', color = 'lightgrey') + # M CDS start
  geom_vline(xintercept = 5717, linetype = 'dashed', color = 'lightgrey') + # M CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (SPEN KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
spen_rdmr_transcript + labs(title = 'SPEN KD CedV M mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(spen_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/SPEN/cedv_spen_rdmr_curve_M.pdf", #EDIT FILENAME!
       width = 1, height = 2.5)

#### F mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_spen_sig_transcript <- dmr_spen_sig %>%
  dplyr::filter(gene_control == 'F')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
spen_rdmr_transcript <- ggplot(dmr_spen_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  

  geom_vline(xintercept = 6129, linetype = 'dashed', color = 'grey25') + # F gene start
  geom_vline(xintercept = 8166, linetype = 'dashed', color = 'grey25') + # F gene end
  geom_vline(xintercept = 6405, linetype = 'dashed', color = 'lightgrey') + # F CDS start
  geom_vline(xintercept = 8078, linetype = 'dashed', color = 'lightgrey') + # F CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (SPEN KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
spen_rdmr_transcript + labs(title = 'SPEN KD CedV F mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(spen_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/SPEN/cedv_spen_rdmr_curve_F.pdf", #EDIT FILENAME!
       width = 1.5, height = 2.5)

#### G mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_spen_sig_transcript <- dmr_spen_sig %>%
  dplyr::filter(gene_control == 'G')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
spen_rdmr_transcript <- ggplot(dmr_spen_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 8170, linetype = 'dashed', color = 'grey25') + # G gene start
  geom_vline(xintercept = 10275, linetype = 'dashed', color = 'grey25') + # G gene end
  geom_vline(xintercept = 8268, linetype = 'dashed', color = 'lightgrey') + # G CDS start
  geom_vline(xintercept = 10136, linetype = 'dashed', color = 'lightgrey') + # G CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (SPEN KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
spen_rdmr_transcript + labs(title = 'SPEN KD CedV G mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(spen_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/SPEN/cedv_spen_rdmr_curve_G.pdf", #EDIT FILENAME!
       width = 1.68, height = 2.5)

#### L mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_spen_sig_transcript <- dmr_spen_sig %>%
  dplyr::filter(gene_control == 'L')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
spen_rdmr_transcript <- ggplot(dmr_spen_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 10279, linetype = 'dashed', color = 'grey25') + # L gene start
  geom_vline(xintercept = 18140, linetype = 'dashed', color = 'grey25') + # L gene end
  geom_vline(xintercept = 10572, linetype = 'dashed', color = 'lightgrey') + # L CDS start
  geom_vline(xintercept = 18077, linetype = 'dashed', color = 'lightgrey') + # L CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (SPEN KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
spen_rdmr_transcript + labs(title = 'SPEN KD CedV L mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(spen_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/SPEN/cedv_spen_rdmr_curve_L.pdf", #EDIT FILENAME!
       width = 6.18, height = 2.5)

#### STM ####
# 1 - Load in CedV single nt positions with MRDs and FDRs
dmr_stm_sig_raw <- read.csv('./NanoCount/STM/stm_control_dmr_tpm_fdr_sig.csv')

# 2 - Assign color value to significant MRD directions (Positive = Red; Negative = Blue)
dmr_stm_sig <- dmr_stm_sig_raw
dmr_stm_sig$mod_direction <- ''
dmr_stm_sig$mod_direction <- ifelse(dmr_stm_sig$normalized_modA_diff > 0, 'pos',
                                     ifelse(dmr_stm_sig$normalized_modA_diff < 0, 'neg', 'unchanged'))

#### N mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_stm_sig_transcript <- dmr_stm_sig %>%
  dplyr::filter(gene_control == 'N')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
stm_rdmr_transcript <- ggplot(dmr_stm_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 56, linetype = 'dashed', color = 'grey25') + # N gene start
  geom_vline(xintercept = 2010, linetype = 'dashed', color = 'grey25') + # N gene end
  geom_vline(xintercept = 144, linetype = 'dashed', color = 'lightgrey') + # N CDS start
  geom_vline(xintercept = 1676, linetype = 'dashed', color = 'lightgrey') + # N CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (stm KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
stm_rdmr_transcript + labs(title = 'stm KD CedV N mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(stm_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/stm/cedv_stm_rdmr_curve_N.pdf", #EDIT FILENAME!
       width = 1.45, height = 2.5)

#### P mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_stm_sig_transcript <- dmr_stm_sig %>%
  dplyr::filter(gene_control == 'P')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
stm_rdmr_transcript <- ggplot(dmr_stm_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 2014, linetype = 'dashed', color = 'grey25') + # P gene start
  geom_vline(xintercept = 4517, linetype = 'dashed', color = 'grey25') + # P gene end
  geom_vline(xintercept = 2112, linetype = 'dashed', color = 'lightgrey') + # P CDS start
  geom_vline(xintercept = 4325, linetype = 'dashed', color = 'lightgrey') + # P CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (stm KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
stm_rdmr_transcript + labs(title = 'stm KD CedV P mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(stm_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/stm/cedv_stm_rdmr_curve_P.pdf", #EDIT FILENAME!
       width = 2.18, height = 2.5)

#### M mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_stm_sig_transcript <- dmr_stm_sig %>%
  dplyr::filter(gene_control == 'M')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
stm_rdmr_transcript <- ggplot(dmr_stm_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 4521, linetype = 'dashed', color = 'grey25') + # M gene start
  geom_vline(xintercept = 6125, linetype = 'dashed', color = 'grey25') + # M gene end
  geom_vline(xintercept = 4635, linetype = 'dashed', color = 'lightgrey') + # M CDS start
  geom_vline(xintercept = 5717, linetype = 'dashed', color = 'lightgrey') + # M CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (stm KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
stm_rdmr_transcript + labs(title = 'stm KD CedV M mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(stm_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/stm/cedv_stm_rdmr_curve_M.pdf", #EDIT FILENAME!
       width = 1, height = 2.5)

#### F mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_stm_sig_transcript <- dmr_stm_sig %>%
  dplyr::filter(gene_control == 'F')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
stm_rdmr_transcript <- ggplot(dmr_stm_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  
  geom_vline(xintercept = 6129, linetype = 'dashed', color = 'grey25') + # F gene start
  geom_vline(xintercept = 8166, linetype = 'dashed', color = 'grey25') + # F gene end
  geom_vline(xintercept = 6405, linetype = 'dashed', color = 'lightgrey') + # F CDS start
  geom_vline(xintercept = 8078, linetype = 'dashed', color = 'lightgrey') + # F CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (stm KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
stm_rdmr_transcript + labs(title = 'stm KD CedV F mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(stm_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/stm/cedv_stm_rdmr_curve_F.pdf", #EDIT FILENAME!
       width = 1.5, height = 2.5)

#### G mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_stm_sig_transcript <- dmr_stm_sig %>%
  dplyr::filter(gene_control == 'G')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
stm_rdmr_transcript <- ggplot(dmr_stm_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 8170, linetype = 'dashed', color = 'grey25') + # G gene start
  geom_vline(xintercept = 10275, linetype = 'dashed', color = 'grey25') + # G gene end
  geom_vline(xintercept = 8268, linetype = 'dashed', color = 'lightgrey') + # G CDS start
  geom_vline(xintercept = 10136, linetype = 'dashed', color = 'lightgrey') + # G CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (stm KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
stm_rdmr_transcript + labs(title = 'stm KD CedV G mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(stm_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/stm/cedv_stm_rdmr_curve_G.pdf", #EDIT FILENAME!
       width = 1.68, height = 2.5)

#### L mRNA ####
# 3 - Filter dataset for individual transcript - EDIT GENE FILTER!
dmr_stm_sig_transcript <- dmr_stm_sig %>%
  dplyr::filter(gene_control == 'L')

# 4 - Plot significant DMRs (FDR < 0.05) on entire genome as visualization reference
stm_rdmr_transcript <- ggplot(dmr_stm_sig_transcript, aes(x = start, y = normalized_modA_diff, fill = mod_direction)) + 
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  # EDIT THIS PORTION OF THE PLOT FOR EACH TRANSCRIPT!
  # Refer to CedV_Genome_modRate_Metagene.R for xintercept positions
  geom_hline(yintercept = 0, linetype = 'solid', color = 'black') + # Y = 0 graph axis
  
  geom_vline(xintercept = 10279, linetype = 'dashed', color = 'grey25') + # L gene start
  geom_vline(xintercept = 18140, linetype = 'dashed', color = 'grey25') + # L gene end
  geom_vline(xintercept = 10572, linetype = 'dashed', color = 'lightgrey') + # L CDS start
  geom_vline(xintercept = 18077, linetype = 'dashed', color = 'lightgrey') + # L CDS end
  
  scale_fill_manual(values = c('pos' = '#B83326FF','neg' = '#1D4497FF', 'unchanged' = 'black')) + 
  geom_smooth(method = loess, se = FALSE, span = 0.05) +
  xlab('Nucleotide Position') + 
  ylab('Differential m6A modification ratio (stm KD - Control)') +
  ylim(-0.30, 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank())
stm_rdmr_transcript + labs(title = 'stm KD CedV L mRNA') 

# 5 - Export pdf plot for each transcript
ggsave(stm_rdmr_transcript, filename = "./NanoCount/plot_output/individual_transcripts/stm/cedv_stm_rdmr_curve_L.pdf", #EDIT FILENAME!
       width = 6.18, height = 2.5)
