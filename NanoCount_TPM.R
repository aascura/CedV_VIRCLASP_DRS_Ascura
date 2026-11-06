# Load in packages
library(dplyr)
library(ggplot2)
library(readr)

# Title: NanoCount_TPM.R (RDMR ModRate Normalized by TPM,Transcript per Million)
# Purpose: Format NanoCount output, calculate average TPM ratio, and normalize calculated RDMR for SPEN and STM samples by TPM calculated by Nanocount, related to Figure 6D-E
# NanoCount reference (associated script): https://github.com/a-slide/NanoCount (NanoCount.sh)
# Author: Andrea Ascura
# Last Update: 18 June 2026

#### SPEN NanoCount ####
# 1 - Load in NanoCount output and add sample column to identify experiment/condition
format_nanocount <- function(path_name, 
                             sample_name){
  nanocount_output <- read_tsv(path_name)
  nanocount_output$sample <- ""
  nanocount_output$sample <- sample_name
  return(nanocount_output)
}

control_1_nanocount <- format_nanocount("./NanoCount/SPEN/transcript_counts_cedv_control_1.tsv", "control_1")
control_2_nanocount <- format_nanocount("./NanoCount/SPEN/transcript_counts_cedv_control_2.tsv", "control_2")
spen_1_nanocount <- format_nanocount("./NanoCount/SPEN/transcript_counts_cedv_spen_1.tsv", "spen_1")
spen_2_nanocount <- format_nanocount("./NanoCount/SPEN/transcript_counts_cedv_spen_2.tsv", "spen_2")

spen_control_combined <- rbind(control_1_nanocount,control_2_nanocount,spen_1_nanocount,spen_2_nanocount)

# 2 - Identify transcript coordinates with viral transcript name
identify_transcript <- function(combined_nanocount_dataset){
  combined_nanocount_dataset$transcript <- ''
  combined_nanocount_dataset$transcript <- ifelse(combined_nanocount_dataset$transcript_name == 'NC_025351.1:56-2010', 'N',
                                                  ifelse(combined_nanocount_dataset$transcript_name == 'NC_025351.1:2014-4517', 'P',
                                                         ifelse(combined_nanocount_dataset$transcript_name == 'NC_025351.1:4521-6125', 'M',
                                                                ifelse(combined_nanocount_dataset$transcript_name == 'NC_025351.1:6129-8166', 'F',
                                                                       ifelse(combined_nanocount_dataset$transcript_name == 'NC_025351.1:8170-10275', 'G',
                                                                              ifelse(combined_nanocount_dataset$transcript_name == 'NC_025351.1:10279-18140','L',
                                                                                     'UNKNOWN'))))))
  return(combined_nanocount_dataset)
}
spen_control_nanocount <- identify_transcript(spen_control_combined)

# 3 - Calculate average TPM between sample replicates (Control_1 & Control_2 / SPEN_1 & SPEN_2) for each transcript
calculate_average_TPM <- function(identified_nanocount_dataset, sample_rep1, sample_rep2, gene){
  split_df <- identified_nanocount_dataset %>%
    dplyr::filter(sample %in% c(sample_rep1,sample_rep2) & transcript == gene) %>%
    dplyr::mutate(average_tpm = mean(tpm)) %>%
    dplyr::relocate(average_tpm, .before = sample) %>%
    dplyr::mutate(TPM_ratio = average_tpm / 1000000) %>%
    dplyr::relocate(TPM_ratio, .after = average_tpm)
  return(split_df)
}
control_1_2_N <- calculate_average_TPM(spen_control_nanocount, 'control_1', 'control_2', 'N')
control_1_2_P <- calculate_average_TPM(spen_control_nanocount, 'control_1', 'control_2', 'P')
control_1_2_M <- calculate_average_TPM(spen_control_nanocount, 'control_1', 'control_2', 'M')
control_1_2_F <- calculate_average_TPM(spen_control_nanocount, 'control_1', 'control_2', 'F')
control_1_2_G <- calculate_average_TPM(spen_control_nanocount, 'control_1', 'control_2', 'G')
control_1_2_L <- calculate_average_TPM(spen_control_nanocount, 'control_1', 'control_2', 'L')

control_1_2_combined <- rbind(control_1_2_N,control_1_2_P,control_1_2_M,control_1_2_F,control_1_2_G,control_1_2_L)

spen_N <- calculate_average_TPM(spen_control_nanocount, 'spen_1', 'spen_2', 'N')
spen_P <- calculate_average_TPM(spen_control_nanocount, 'spen_1', 'spen_2', 'P')
spen_M <- calculate_average_TPM(spen_control_nanocount, 'spen_1', 'spen_2', 'M')
spen_F <- calculate_average_TPM(spen_control_nanocount, 'spen_1', 'spen_2', 'F')
spen_G <- calculate_average_TPM(spen_control_nanocount, 'spen_1', 'spen_2', 'G')
spen_L <- calculate_average_TPM(spen_control_nanocount, 'spen_1', 'spen_2', 'L')

spen_combined <- rbind(spen_N, spen_P, spen_M, spen_F, spen_G, spen_L)

# 4 - Read in RDMR file with transcript/gene identification
dmr_spen <- read.csv('./Nanopore_DRS/cedv_spen_control_rdmr_annotated_regions.csv')

# Split RDMR and add average TPM for corresponding genes from both datasets (RDMR & TPM outputs)
# SPEN Controls (Control_1 & Control_2), SPEN_1 & SPEN_2
combine_modA_tpm <- function(dmr_dataset,
                             tpm_combined_dataset,
                             count_column,
                             total_column, 
                             modA_column) {
  combined_df <- dmr_dataset %>%
    dplyr::select(c(start,end,count_column,total_column, modA_column,gene))
  combined_df$average_tpm <- ''
  combined_df$average_tpm <- tpm_combined_dataset$average_tpm[match(combined_df$gene, tpm_combined_dataset$transcript)]
  combined_df$TPM_ratio <- ''
  combined_df$TPM_ratio <- tpm_combined_dataset$TPM_ratio[match(combined_df$gene, tpm_combined_dataset$transcript)]
  return(combined_df)
}

# sample_a = Control
# sample_b = SPEN/STM
control_1_2_modA_tpm <- combine_modA_tpm(dmr_spen,
                                         control_1_2_combined,
                                         'sample_a_counts', 
                                         'sample_a_total', 
                                         'modA_ratio_control')
                                         
spen_modA_tpm <- combine_modA_tpm(dmr_spen, 
                                  spen_combined,
                                  'sample_b_counts', 
                                  'sample_b_total', 
                                  'modA_ratio_spen')

# 5 - Compute normalized modA_ratio by multiplying by TPM ratio
control_1_2_normalization <- control_1_2_modA_tpm %>%
  dplyr::mutate(normalized_modA = (modA_ratio_control * TPM_ratio)) %>%
  dplyr::mutate(log_normalized_modA = log2(normalized_modA))

spen_normalization <- spen_modA_tpm %>%
  dplyr::mutate(normalized_modA = (modA_ratio_spen * TPM_ratio)) %>%
  dplyr::mutate(log_normalized_modA = log2(normalized_modA))

# 6 - Combined TPM-normalized Control_1/Control_2 dataset with TPM-normalized SPEN_1/SPEN_2 dataset
spen_control_tpm_normalized <- inner_join(control_1_2_normalization, spen_normalization,
                                          join_by(start == start, end == end),
                                          suffix = c("_control","_spen"))

# 7 - Export dataset as csv for calculation of RDMR and statistics in following script: 
write.csv(spen_control_tpm_normalized, './NanoCount/SPEN/spen_control_tpm_normalized.csv', row.names = FALSE)

# END - Proceed to TPM_Normalized_RDMR_Add_Stats.R using 'spen_control_tpm_normalized.csv' as the input for SPEN samples

#### STM NanoCount ####

# 1 - Load in NanoCount output and add sample column to identify experiment/condition
control_3_nanocount <- format_nanocount("./NanoCount/STM/transcript_counts_cedv_control_3.tsv", "control_3")
control_4_nanocount <- format_nanocount("./NanoCount/STM/transcript_counts_cedv_control_4.tsv", "control_4")
stm_1_nanocount <- format_nanocount("./NanoCount/STM/transcript_counts_cedv_stm_1.tsv", "stm_1")
stm_2_nanocount <- format_nanocount("./NanoCount/STM/transcript_counts_cedv_stm_2.tsv", "stm_2")

stm_control_combined <- rbind(control_3_nanocount,control_4_nanocount,stm_1_nanocount,stm_2_nanocount)

# 2 - Identify transcript coordinates with viral transcript name
stm_control_nanocount <- identify_transcript(stm_control_combined)

# 3 - Calculate average TPM between sample replicates (Control_1 & Control_2 / SPEN_1 & SPEN_2) for each transcript
control_3_4_N <- calculate_average_TPM(stm_control_nanocount, 'control_3', 'control_4', 'N')
control_3_4_P <- calculate_average_TPM(stm_control_nanocount, 'control_3', 'control_4', 'P')
control_3_4_M <- calculate_average_TPM(stm_control_nanocount, 'control_3', 'control_4', 'M')
control_3_4_F <- calculate_average_TPM(stm_control_nanocount, 'control_3', 'control_4', 'F')
control_3_4_G <- calculate_average_TPM(stm_control_nanocount, 'control_3', 'control_4', 'G')
control_3_4_L <- calculate_average_TPM(stm_control_nanocount, 'control_3', 'control_4', 'L')

control_3_4_combined <- rbind(control_3_4_N,control_3_4_P,control_3_4_M,control_3_4_F,control_3_4_G,control_3_4_L)

stm_N <- calculate_average_TPM(stm_control_nanocount, 'stm_1', 'stm_2', 'N')
stm_P <- calculate_average_TPM(stm_control_nanocount, 'stm_1', 'stm_2', 'P')
stm_M <- calculate_average_TPM(stm_control_nanocount, 'stm_1', 'stm_2', 'M')
stm_F <- calculate_average_TPM(stm_control_nanocount, 'stm_1', 'stm_2', 'F')
stm_G <- calculate_average_TPM(stm_control_nanocount, 'stm_1', 'stm_2', 'G')
stm_L <- calculate_average_TPM(stm_control_nanocount, 'stm_1', 'stm_2', 'L')

stm_combined <- rbind(stm_N, stm_P, stm_M, stm_F, stm_G, stm_L)

# 4 - Read in RDMR file with transcript/gene identification
dmr_stm <- read.csv('./Nanopore_DRS/cedv_stm_control_rdmr_annotated_regions.csv')

# Split RDMR and add average TPM for corresponding genes from both datasets (RDMR & TPM outputs)
# STM Controls (Control_3 & Control_4), STM_1 & STM_2

# sample_a = Control
# sample_b = STM
control_3_4_modA_tpm <- combine_modA_tpm(dmr_stm,
                                         control_3_4_combined,
                                         'sample_a_counts', 
                                         'sample_a_total', 
                                         'modA_ratio_control')

stm_modA_tpm <- combine_modA_tpm(dmr_stm, 
                                  stm_combined,
                                  'sample_b_counts', 
                                  'sample_b_total', 
                                  'modA_ratio_stm')

# 5 - Compute normalized modA_ratio by multiplying by TPM ra
control_3_4_normalization <- control_3_4_modA_tpm %>%
  dplyr::mutate(normalized_modA = (modA_ratio_control * TPM_ratio)) %>%
  dplyr::mutate(log_normalized_modA = log2(normalized_modA))

stm_normalization <- stm_modA_tpm %>%
  dplyr::mutate(normalized_modA = (modA_ratio_stm * TPM_ratio)) %>%
  dplyr::mutate(log_normalized_modA = log2(normalized_modA))

# 6 - Combined TPM-normalized Control_3/Control_4 dataset with TPM-normalized STM_1/STM_2 dataset
stm_control_tpm_normalized <- inner_join(control_3_4_normalization, stm_normalization,
                                          join_by(start == start, end == end),
                                          suffix = c("_control","_stm"))

# 7 - Export dataset as csv for calculation of RDMR and statistics in following script: 
write.csv(stm_control_tpm_normalized, './NanoCount/STM/stm_control_tpm_normalized.csv', row.names = FALSE)

# END - Proceed to TPM_Normalized_RDMR_Add_Stats.R using 'stm_control_tpm_normalized.csv' as the input for STM samples

all_tpm_spen <- rbind(control_1_2_combined, spen_combined) %>%
  dplyr::relocate('sample', .before = 'tpm')

all_tpm_stm <- rbind(control_3_4_combined, stm_combined) %>%
  dplyr::relocate('sample', .before = 'tpm')

