# Load in packages
library(dplyr)
library(stringr)
library(UniProt.ws)

# Title: CedV_VIRCLASP_Protein_Filtering.R
# Author: Andrea Ascura
# Last Update: 18 November 2025

# Description: Filter and format raw csv file containing CLASP_4SU/Input_4SU proteins output from CedV_VIRCLASP_MSstats.R

# 1 - Load in MSstats groupComparison output and MaxQuant proteinGroups file
pg <- read.table('./VIRCLASP_Proteomics/maxquant_output/proteinGroups.txt', sep = '\t', header = TRUE)

# 2 - Extract Uniprot ID from Fasta.header column
protein_pg <- within(pg, Protein <- str_split_i(Protein.IDs,';',1)) %>%
  dplyr::select(Protein.names, Gene.names, Fasta.headers, Protein)

# 3 - Load in unfiltered CedV VIRCLASP interactome dataset output from CedV_MSstats.R
cedv_raw <- read.csv('./VIRCLASP_Proteomics/msstats_output/Unfiltered_CedV_VIRCLASP_Hits_AA.csv', header = TRUE)

# 4 - Filter CedV VIRCLASP interactome to remove non-numeric p-values
cedv_filt <- cedv_raw %>%
  dplyr::filter(!is.na(adj.pvalue)) %>%
  arrange(desc(log2FC))

# 5 - Identify VIRCLASP hits using extracted Uniprot ID
# This contains Swiss-prot (sp) and isoforms (tr) from ALL species
cedv_hits_sp_tr <- left_join(cedv_filt,protein_pg, by = 'Protein') %>%
  relocate(Gene.names, Protein.names, .before = log2FC) %>%
  relocate(Protein, .before = Fasta.headers) 

# 6 - Fill in empty cells by passing Uniprot_ID in Uniprot ID mapping service (UniProt.ws)
uni_id <- cedv_hits_sp_tr['Protein']
up <- UniProt.ws(taxId=9606)

res <- UniProt.ws::select(
  x = up,
  keys = uni_id,
  columns = c("gene_primary", 
              "gene_names", 
              "protein_name"),
  keytype = "UniProtKB"
) %>% dplyr::select(-From)

colnames(res) <- c('Uniprot_ID','Primary_gene_name', 'Gene_names','Protein_names')

# 7 - Merge Uniprot mapping results with CedV Hits (Sp and Tr)
CedV_VIRCLASP_hits_all <- merge(cedv_hits_sp_tr, res, by.x = 'Protein', by.y = 'Uniprot_ID') %>%
  dplyr::select(Protein, Primary_gene_name, log2FC, adj.pvalue, Gene_names, Protein_names, Fasta.headers) %>%
  rename('Uniprot_ID' = 'Protein') %>%
  arrange(desc(log2FC))

# 8 - Export filtered CedV VIRCLASP list including Swiss-prot (sp) and Trembl (tr) proteins
# This file contains proteins from all species used in proteome searches: Human (Homo Sapiens), Monkey (Chlorocebus sabaeus), and Cedar virus
write.csv(CedV_VIRCLASP_hits_all,'./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_sp_tr_AA.csv',row.names = FALSE)

# 9 - Filter CedV VIRCLASP list for Human Swiss-Prot ('sp') proteins only
CedV_VIRCLASP_hits_human <- within(CedV_VIRCLASP_hits_all, Form <- substr(Fasta.headers,1,2)) %>%
  filter(Form == 'sp') %>%
  dplyr::select(-c(Form)) %>%
  relocate(Uniprot_ID, .before = Fasta.headers) %>%
  arrange(desc(log2FC))

# 10 - Export filtered CedV VIRCLASP list of Human Swiss-Prot proteins
write.csv(CedV_VIRCLASP_hits_human,'./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_human_sp_AA.csv',row.names = FALSE)

# 11 - Load back in VIRCLASP hits human and remove -Inf/Inf log2FC to produce table of all CLASP proteins used for future analyses
# This is the list of 1436 distinct human proteins identified in the CLASP sample
cedv_human_rbp_raw <- fread('./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_human_sp_AA.csv')
cedv_human_rbp <- cedv_human_rbp_raw %>%
  dplyr::filter(!log2FC == Inf) %>%
  dplyr::filter(!log2FC == -Inf) %>%
  distinct(Primary_gene_name, .keep_all = TRUE)

# 12 - Export list of filtered CLASP proteins
write.csv(cedv_human_rbp,'./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_human_AA.csv',row.names = FALSE)

