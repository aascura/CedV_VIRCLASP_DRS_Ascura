# Load in packages
library(dplyr)
library(data.table)
library(EnhancedVolcano)

# Title: VIRCLASP_Volcano.R
# Author: Andrea Ascura
# Last Update: 18 November 2025

# Description: Plot depleted and enriched CLASP proteins using volcano plot

# 1 - Load in CedV VIRCLASP hits formatted from CedV_VIRCLASP_Interactome.R and remove Inf/-Inf log2FC values
cedv_all_rbp <- fread('./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_sp_tr_AA.csv') %>% # Contains Sp & Tr for all Species
  dplyr::filter(!log2FC == Inf) %>%
  dplyr::filter(!log2FC == -Inf) %>%

# 2 - Filter VIRCLASP hits for only Sp & Tr Human and CedV Viral Proteins and recombine into single dataset
cedv_sp_tr_human <- cedv_all_rbp[grepl('_HUMAN',
                                       cedv_all_rbp$Fasta.headers,
                                       ignore.case = TRUE)]
cedv_viral_pro <- cedv_all_rbp[grepl('OS=Cedar virus', 
                                     cedv_all_rbp$Fasta.headers,
                                     ignore.case = TRUE)]
cedv_rbps <- rbind(cedv_sp_tr_human, cedv_viral_pro)

# 3 - Assign colors based on fold-change and p-value cutoffs
# Keyvals for Human RBPs + CedV Viral Proteins - Complete parameters preserved below
keyvals <- ifelse(cedv_rbps$Primary_gene_name %in% cedv_viral_pro$Primary_gene_name, '#EDB144FF',
                  ifelse(cedv_rbps$log2FC < -1 & cedv_rbps$adj.pvalue < 0.05, '#1D4497FF', # FC and P cutoff
                         ifelse(cedv_rbps$log2FC > 1 & cedv_rbps$adj.pvalue < 0.05, '#B83326FF','#BAB0AC'))) # FC and P cutoff

keyvals[is.na(keyvals)] <- '#BAB0AC'
names(keyvals)[keyvals == '#EDB144FF'] <- 'CedV proteins'
names(keyvals)[keyvals == '#B83326FF'] <- 'Enriched'
names(keyvals)[keyvals == '#1D4497FF'] <- 'Depleted'

# Optional Step: Make modifications for keyvals to explicitly identify specific candidates
# Skip to Step 4 if no candidates will be identified in separate text bubble
siRNA_candidates <- c('ATRX','CHD1','IFI16','POM121C','SPEN','RBM15','RBM15B','WDR33','SRRM2','AFF4','TTF2') # Edit this list for different candidates

cedv_nocandidates <- cedv_rbps %>%
  filter(!Primary_gene_name %in% siRNA_candidates)
cedv_candidates <- cedv_rbps %>%
  filter(Primary_gene_name %in% siRNA_candidates)
cedv_rbps <- rbind(cedv_nocandidates,cedv_candidates) # Rows appended to bottom will graph point last (bring point to front)

keyvals <- ifelse(cedv_rbps$Primary_gene_name %in% siRNA_candidates, '#f39202',
                  ifelse(cedv_rbps$log2FC < -1 & cedv_rbps$adj.pvalue < 0.05, '#1D4497FF', # FC and P cutoff
                         ifelse(cedv_rbps$log2FC > 1 & cedv_rbps$adj.pvalue < 0.05, '#B83326FF','#BAB0AC'))) # FC and P cutoff

keyvals[is.na(keyvals)] <- '#BAB0AC'
names(keyvals)[keyvals == '#f39202'] <- 'siRNA candidates'
names(keyvals)[keyvals == '#B83326FF'] <- 'Enriched'
names(keyvals)[keyvals == '#1D4497FF'] <- 'Depleted'

# 4 - Create volcano plot of enriched and depleted proteins: Sp & Tr Human Proteins + CedV Proteins
cedv_volcano <- EnhancedVolcano(cedv_rbps,
                                lab = cedv_rbps$Primary_gene_name,
                                x = 'log2FC', y = 'adj.pvalue',
                                ylim = c(0, max(-log10(cedv_rbps[["adj.pvalue"]]), na.rm = TRUE) + 1.5),
                                xlab = bquote(~log[2] ~ "fold change (CLASP/Input)"),
                                ylab = bquote(~Adjusted ~ italic(P) ~value ~(-log[10] ~ italic(P))),
                                title = NULL,
                                subtitle = NULL,
                                caption = NULL,
                                selectLab = ' ', # siRNA_candidates 
                                drawConnectors = TRUE,
                                lengthConnectors = unit(5.0, "npc"),
                                directionConnectors = "both",
                                widthConnectors = 1,
                                arrowheads = FALSE,
                                colCustom = keyvals,
                                #colAlpha = 1,
                                pCutoff = 0.05,
                                FCcutoff = 1, # Fold change cutoff
                                cutoffLineType = 'dashed',
                                pointSize = 4,
                                legendPosition = 'right',
                                boxedLabels = TRUE,
                                legendLabels = c('Depleted','Enriched'),
                                gridlines.major = FALSE, gridlines.minor = FALSE)

cedv_volcano
cedv_volcano + theme(legend.position = 'none') # optional viewing without legend

# 5 - Save plot to PDF and specify dimensions
ggsave(
  filename = 'VIRCLASP_Volcano.pdf',
  plot = cedv_volcano, # last_plot()
  path = './VIRCLASP_Proteomics/plot_outputs',
  width = 9,
  height = 8,
  dpi = 300
)
