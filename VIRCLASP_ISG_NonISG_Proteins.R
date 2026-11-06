# Load in packages
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(rstatix)

# Title: VIRCLASP_ISG_NonISG_Proteins.R
# Author: Andrea Ascura
# Last Update: 18 November 2025

# Description: Determine enrichment of proteins annotated as ISGs in all CLASP protein sample and CedV RNA interactome
# ISG annotation database (Interferome): https://interferome.org/interferome/home.jspx

# Optional (Skip to Step 1 if already performed): 
# Load in VIRCLASP list of human proteins and filter for CedV RNA interactome to generate required input for Interferome database query.
all_clasp <- fread('./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_human_AA.csv')
clasp_genes <- paste(all_clasp$Primary_gene_name, sep = " ", collapse = ',')
write.table(clasp_genes,'./VIRCLASP_Proteomics/reference_datasets/CedV_All_CLASP_Proteins.txt')

# 1 - Load in CLASP IFNome search results and CedV VIRCLASP list
clasp_ifnome_raw <- fread('./VIRCLASP_Proteomics/reference_datasets/IFNome_CLASP_GeneSearchResults.txt')
clasp_ifnome <- clasp_ifnome_raw %>%
  distinct(`Gene Name`, .keep_all = TRUE)

cedv_human_rbp <- fread('./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_human_AA.csv') %>%
  dplyr::rename('SYMBOL' = 'Primary_gene_name') %>%
  distinct(SYMBOL, .keep_all = TRUE)

# 2 - Retrieve gene symbols from CLASP proteins Uniprot IDs
clasp_proteins <- fread('./VIRCLASP_Proteomics/reference_datasets/CLASP_Proteins.csv') %>%
  distinct(Primary_gene_name, .keep_all = TRUE)

clasp_ids_raw <- bitr(geneID = clasp_proteins$Var1,
                      fromType = 'UNIPROT', 
                      toType = c('SYMBOL'),
                      OrgDb = 'org.Hs.eg.db')
clasp_ids <- clasp_ids_raw %>%
  distinct(UNIPROT, .keep_all = TRUE)

# 3 - Use CLASP IFNome results to designate ISGs and Non-ISGs CLASP proteins
clasp_isgs_list <- clasp_ids %>%
  filter(SYMBOL %in% clasp_ifnome$`Gene Name`) %>%
  left_join(cedv_human_rbp) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(Gene_names,log2FC,adj.pvalue)
clasp_isgs_list$type <- 'ISG'

clasp_nonisgs_list <- clasp_ids %>%
  filter(!SYMBOL %in% clasp_ifnome$`Gene Name`) %>%
  left_join(cedv_human_rbp) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(Gene_names,log2FC,adj.pvalue)
clasp_nonisgs_list$type <- 'Non-ISG'

# 4 - Use CLASP IFNome results to designate ISGs and Non-ISGs CLASP proteins for CDF preparation
clasp_isgs <- clasp_ids %>%
  distinct(SYMBOL) %>%
  filter(SYMBOL %in% clasp_ifnome$`Gene Name`) %>%
  left_join(cedv_human_rbp) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(log2FC)
clasp_isgs$type <- 'ISG'
clasp_nonisgs <- clasp_ids %>%
  filter(!SYMBOL %in% clasp_ifnome$`Gene Name`) %>%
  left_join(cedv_human_rbp) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(log2FC)
clasp_nonisgs$type <- 'Non-ISG'

# 5 - Create all CLASP protein baseline by combining ISG and Non-ISG
clasp_all <- rbind(clasp_isgs,clasp_nonisgs)
clasp_all$type <- 'all'

# 6 - Combine ISG, Non-ISG, and all CLASP proteins for CDF prep
combined_clasp <- rbind(clasp_all,clasp_isgs,clasp_nonisgs)

# 7 - Calculate significance of groups ('Control' = 'All CLASP Proteins') using Wilcoxon test
pairwise_wilcox_test(combined_clasp,log2FC ~ type,ref.group = "all",p.adjust.method = "BH")

# 8 - Visualize using CDF
p1 <- ggplot(combined_clasp, aes(log2FC, colour = type)) +
  stat_ecdf() +
  scale_color_manual(name = '',
                     values = c('#BAB0AC', 
                                '#C8570DFF',
                                '#5773C0FF'),
                     labels = c(paste0('All CLASP Proteins (',nrow(clasp_all),')'),
                                paste0('ISGs (', nrow(clasp_isgs),')'),
                                paste0('Non-ISGs (', nrow(clasp_nonisgs),')'))) +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) #+ 
#theme(legend.position = 'none')
p1

# 9 - Export plot to PDF and specify dimensions
ggsave(
  filename = 'CLASP_ISGs_CDF.pdf',
  plot = p1,
  path = './VIRCLASP_Proteomics/plot_outputs',
  width = 6,
  height = 6,
  dpi = 300
)

# 10 - Load in CedV IFNome results
cedv_ifnome_raw <- fread('./VIRCLASP_Proteomics/reference_datasets/IFNome_CedVInteractome_GeneSearchResults.txt')
cedv_ifnome <- cedv_ifnome_raw %>%
  distinct(`Gene Name`, .keep_all = TRUE)

# 11 - Load in CedV Interactome List
cedv_human_rbp_sig <- cedv_human_rbp %>%
  filter(log2FC > 0 & adj.pvalue <= 0.05) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# 12 - Use CedV Interactome IFNome results to designate ISGs and Non-ISGs CedV proteins
cedv_isgs <- cedv_human_rbp_sig %>%
  filter(SYMBOL %in% cedv_ifnome$`Gene Name`) %>%
  dplyr::select(log2FC)
cedv_isgs$type <- 'ISG'
cedv_nonisgs <- cedv_human_rbp_sig %>%
  filter(!SYMBOL %in% cedv_ifnome$`Gene Name`) %>%
  dplyr::select(log2FC)
cedv_nonisgs$type <- 'Non-ISG'

# 13 - Create all CedV protein baseline by combining ISG and Non-ISG
cedv_all <- rbind(cedv_isgs,cedv_nonisgs)
cedv_all$type <- 'all'

# 14 - Combine ISG, Non-ISG, and all CedV proteins
combined_cedv <- rbind(cedv_all, cedv_isgs, cedv_nonisgs)

# 15 - Calculate significance of groups ('Control' = 'All CedV Proteins') using Wilcoxon test
pairwise_wilcox_test(combined_cedv,log2FC ~ type,ref.group = "all",p.adjust.method = "BH")

# 16 - Visualize using CDF
p2 <-ggplot(combined_cedv, aes(log2FC, colour = type)) +
  stat_ecdf() +
  scale_color_manual(name = '',
                     values = c('#BAB0AC', 
                                '#C8570DFF',
                                '#5773C0FF'),
                     labels = c(paste0('CedV Interactome (',nrow(cedv_all),')'),
                                paste0('ISGs (', nrow(cedv_isgs),')'),
                                paste0('Non-ISGs (', nrow(cedv_nonisgs),')'))) +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = 'none')
p2

# 17 - Export plot to PDF and specify dimensions
ggsave(
  filename = 'CedV_Interactome_ISGs_CDF.pdf',
  plot = p2,
  path = './VIRCLASP_Proteomics/plot_outputs',
  width = 6,
  height = 6,
  dpi = 300
)