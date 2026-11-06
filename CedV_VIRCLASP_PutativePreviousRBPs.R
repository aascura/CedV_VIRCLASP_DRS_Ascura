# Load in packages
library(dplyr)
library(readxl)

# Title: CedV_VIRCLASP_PutativePreviousRBPs.R
# Author: Andrea Ascura
# Last Update: 20 November 2025

# Description: Determine proportion of CedV RNA interactome that are putative RBPs vs previously identified RBPs based on other vRNA interactomes

# 1 - Load in list of CLASP proteins
cedv_human_rbp <- fread('./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_human_AA.csv')

# 2 - Filter CLASP list to obtain CedV RNA interactome candidates
cedv_human_rbp_sig <- cedv_human_rbp %>%
  filter(log2FC > 1 & adj.pvalue <= 0.05) %>%
  dplyr::select('Primary_gene_name','Uniprot_ID','log2FC') 

# 3 - Load in reference list of RBPs identified on vRNA compiled from Castello et al., Mol Cell, 2024
castello_vRNA_Interactome <- read_excel('./VIRCLASP_Proteomics/reference_datasets/Castello_MolCell_2024.xlsx', sheet = 'vRNA interactome')
colnames(castello_vRNA_Interactome) <- castello_vRNA_Interactome[1,]
castello_vRNA_Interactome <- castello_vRNA_Interactome[-1,]
castello_select <- castello_vRNA_Interactome %>%
  dplyr::select(`Gene name`)

# 4 - Create custom function to check for RBPs in vRNPs (Castello Reference Only). This can be adapted to include other datasets.
virclasp_hit_checker_vRNP <- function(virclasp_list,castello_dataset) {
  virclasp_overlaps_vRNP <- virclasp_list
  virclasp_overlaps_vRNP$Previous_RBP <- virclasp_list[[1]] %in% castello_dataset[[1]]
  write.csv(virclasp_overlaps_vRNP,'./VIRCLASP_Proteomics/reference_datasets/CedV_VIRCLASP_PreviousPutativeRBPs.csv', row.names = FALSE)
  return(virclasp_overlaps_vRNP)
}

# 5 - Run function on CedV VIRCLASP Hits
res <- virclasp_hit_checker_vRNP(cedv_human_rbp_sig, castello_select)

# 6 - Determine Number of Identified vs Novel RBPs
cedv_identified <- res %>%
  dplyr::select(Primary_gene_name, Castello_RBP) %>%
  filter(Castello_RBP == TRUE) 

cedv_identified <- cedv_human_rbp_sig %>% # Reassociate primary gene name with Uniprot_ID
  filter(Primary_gene_name %in% cedv_identified$Primary_gene_name)
cedv_novel <- cedv_human_rbp_sig %>%
  filter(!Primary_gene_name %in% cedv_identified$Primary_gene_name)

# 7 - Visualize proportion of previously identified RBPs vs Putative RBPs
status <-  c('Previously identified RBP','Putative RBP')
count <- c(nrow(cedv_identified),nrow(cedv_novel))
identified_df <- data.frame(Status = status, Count = count)
identified_df$Category <- 'CedV Interactome' # Add to create a stacked bar plot

ggplot(identified_df, aes(fill = Status, y = Count, x = Category)) + 
  geom_bar(position = 'fill', stat = 'identity', colour = 'black', width = 0.4) + 
  labs(x = "", y = 'Proportion') +
  scale_fill_manual(values = c('#8CC8BC', '#d2e9e4')) + 
  theme_classic() +
  labs(fill = '') +
  geom_text(aes(label = Count),
            position=position_fill(vjust=0.5), colour="black", size = 5)

# 8 - Save plot to PDF and specify dimensions
ggsave(
  filename = 'Cedv_PutativePreviousRBP_Barplot.pdf',
  plot = last_plot(), # last_plot()
  path = './VIRCLASP_Proteomics/plot_outputs',
  width = 6,
  height = 8,
  dpi = 300
)
