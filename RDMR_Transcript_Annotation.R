# Load in packages
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(stringr)

# Title: RDMR_Transcript_Annotation.R
# Purpose: Format output from Modkit_Differential_Methylation_Scoring.sh (must be already converted from BED to CSV/XLSX), compute unweighted RDMR for downstream TPM normalization, assign transcript position to start sites, and apply coverage filter.
# NanoCount reference (associated script): https://github.com/a-slide/NanoCount (NanoCount.sh)
# Author: Andrea Ascura
# Last Update: 05 November 2026

#### SPEN Unweighted RDMR Calculation and Transcript Annotations ####

# sample_a = Control
# sample_b = SPEN

# 1 - Load in differential analysis (dmr_output) from Modkit_Differential_Methylation_Scoring.sh and remove unnecessary strings from numeric columns
# Modkit descriptions to dmr single base analysis column names: https://nanoporetech.github.io/modkit/intro_dmr.html
dmr_raw <- read_excel('./Nanopore_DRS/cedv_bmt0.7_dmr_control_spen.xlsx')

dmr <- dmr_raw %>%
  mutate(sample_a_percents = str_remove(sample_a_percents,'a:')) %>%
  mutate(sample_a_counts = str_remove(sample_a_counts,'a:')) %>%
  mutate(sample_b_percents = str_remove(sample_b_percents,'a:')) %>%
  mutate(sample_b_counts = str_remove(sample_b_counts,'a:'))

# 2 - Convert column classes from character to numeric and compute percent difference
dmr$sample_a_percents <- as.numeric(as.character(dmr$sample_a_percents))
dmr$sample_b_percents <- as.numeric(as.character(dmr$sample_b_percents))
dmr$sample_a_counts <- as.numeric(as.character(dmr$sample_a_counts))
dmr$sample_b_counts <- as.numeric(as.character(dmr$sample_b_counts))

# 3 - Calculate modification ratio difference (mod_diff): modA_ratio_spen - modA_ratio_control
dmr_full_spen <- dmr %>%
  mutate(modA_ratio_control = sample_a_counts/sample_a_total) %>% # p1; sample_a = control
  mutate(modA_ratio_spen = sample_b_counts/sample_b_total) %>% # p2; sample_b = spen
  mutate(mod_diff = modA_ratio_spen - modA_ratio_control) %>%
  dplyr::select(start, end, sample_a_counts, sample_a_total, sample_b_counts, sample_b_total, modA_ratio_control, modA_ratio_spen, mod_diff)

# 4 - Annotate start positions with transcript location
#  - Define CedV genomic boundaries (gene start + end) (nucleotide positions) to subset data into regions and annotate full DMR
# https://www.ncbi.nlm.nih.gov/nuccore/JQ001776
# N: 56...2010
# P: 2014...4517
# M: 4521...6125
# F: 6129...8166
# G: 8170...10275
# L: 10279...18140

dmr_full_spen$gene <- ''
dmr_full_spen$gene <- ifelse(dmr_full_spen$start >= 56 & dmr_full_spen$end <= 2010, 'N',
                             ifelse(dmr_full_spen$start >= 2014 & dmr_full_spen$end <= 4517, 'P',
                                    ifelse(dmr_full_spen$start >= 4521 & dmr_full_spen$end <= 6125, 'M',
                                           ifelse(dmr_full_spen$start >= 6129 & dmr_full_spen$end <= 8166, 'F',
                                                  ifelse(dmr_full_spen$start >= 8170 & dmr_full_spen$end <= 10275, 'G',
                                                         ifelse(dmr_full_spen$start >= 10279 & dmr_full_spen$end <= 18140, 'L',
                                                                'UNKNOWN'))))))

# 5 - Filter for sites with at least 50 modified or unmodified A's at a single position
dmr_cov50_spen <- dmr_full_spen %>%
  dplyr::filter(sample_a_total >= 50 & sample_b_total >= 50) # 452 sites

# 6 - Export file dmr_full file with gene locations annotated
write.csv(dmr_full_spen, './Nanopore_DRS/cedv_spen_control_rdmr_annotated_regions.csv', row.names = FALSE)

#### STM2457 Unweighted RDMR Calculation and Transcript Annotations ####

# sample_a = Control
# sample_b = STM2457

# 1 - Load in differential analysis (dmr_output) from Modkit_Differential_Methylation_Scoring.sh and remove unnecessary strings from numeric columns
# Modkit descriptions to dmr single base analysis column names: https://nanoporetech.github.io/modkit/intro_dmr.html
dmr_raw_stm <- read_excel('./Nanopore_DRS/cedv_bmt0.7_dmr_control_stm.xlsx')

dmr_stm <- dmr_raw_stm %>%
  mutate(sample_a_percents = str_remove(sample_a_percents,'a:')) %>%
  mutate(sample_a_counts = str_remove(sample_a_counts,'a:')) %>%
  mutate(sample_b_percents = str_remove(sample_b_percents,'a:')) %>%
  mutate(sample_b_counts = str_remove(sample_b_counts,'a:'))

# 2 - Convert column classes from character to numeric and compute percent difference
dmr_stm$sample_a_percents <- as.numeric(as.character(dmr_stm$sample_a_percents))
dmr_stm$sample_b_percents <- as.numeric(as.character(dmr_stm$sample_b_percents))
dmr_stm$sample_a_counts <- as.numeric(as.character(dmr_stm$sample_a_counts))
dmr_stm$sample_b_counts <- as.numeric(as.character(dmr_stm$sample_b_counts))

# 3 - Calculate modification ratio difference (mod_diff): modA_ratio_stm - modA_ratio_control
dmr_full_stm <- dmr_stm %>%
  mutate(modA_ratio_control = sample_a_counts/sample_a_total) %>% # p1; sample_a = control
  mutate(modA_ratio_stm = sample_b_counts/sample_b_total) %>% # p2; sample_b = stm
  mutate(mod_diff = modA_ratio_stm - modA_ratio_control) %>%
  dplyr::select(start, end, sample_a_counts, sample_a_total, sample_b_counts, sample_b_total, modA_ratio_control, modA_ratio_stm, mod_diff)

dmr_full_stm$gene <- ''
dmr_full_stm$gene <- ifelse(dmr_full_stm$start >= 56 & dmr_full_stm$end <= 2010, 'N',
                            ifelse(dmr_full_stm$start >= 2014 & dmr_full_stm$end <= 4517, 'P',
                                   ifelse(dmr_full_stm$start >= 4521 & dmr_full_stm$end <= 6125, 'M',
                                          ifelse(dmr_full_stm$start >= 6129 & dmr_full_stm$end <= 8166, 'F',
                                                 ifelse(dmr_full_stm$start >= 8170 & dmr_full_stm$end <= 10275, 'G',
                                                        ifelse(dmr_full_stm$start >= 10279 & dmr_full_stm$end <= 18140, 'L',
                                                               'UNKNOWN'))))))

# 4 - Filter for sites with at least 50 modified or unmodified A's at a single position
dmr_cov50_stm <- dmr_full_stm %>%
  dplyr::filter(sample_a_total >= 50 & sample_b_total >= 50) # 505

# 5 - Export file dmr_full file with gene locations annotated
write.csv(dmr_full_stm, './Nanopore_DRS/cedv_stm_control_rdmr_annotated_regions.csv', row.names = FALSE)
