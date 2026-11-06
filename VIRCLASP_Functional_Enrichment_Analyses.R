# Load in packages
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)

# Title: VIRCLASP_Functional_Enrichment_Analyses.R
# Author: Andrea Ascura
# Last Update: 18 November 2025

# Description: Create custom plots using output from various function enrichment analyses (InterPro Features, Reactome)

## InterPro Features Enrchment
# 1 - Load in InterPro enrichment data output from STRING database (https://string-db.org/)
# STRING database was queried using CedV RBPs from VIRCLASP list (log2FC > 0) consisting of 203 proteins
interpro <- fread('./VIRCLASP_Proteomics/reference_datasets/CedV_STRING_InterPro_ProteinDomains.tsv') %>%
  dplyr::arrange(desc(`observed gene count`)) %>%
  dplyr::slice_head(n = 11)

# 2 - Convert FDR column to numeric class
interpro$`false discovery rate` <- as.numeric(as.character(interpro$`false discovery rate`))

# 3 - Plot InterPro enrichment data on a lollipop plot
interpro_plt <- ggplot(interpro, aes(x = reorder(`term description`, `observed gene count`), y = `strength`)) +
  geom_segment(aes(x = reorder(`term description`, `observed gene count`),
                   y = 0, yend = `strength`),
               color = "black", lwd = 1) +
  geom_point(aes(size = `observed gene count`, colour = `false discovery rate`)) +
  scale_size(range = c(2,10)) +
  scale_color_gradient(low = "#9CC184FF", high = "#1E3D14FF") +
  xlab('') +
  ylab("Strength") +
  coord_flip() +
  theme_minimal()
interpro_plt

# 4 - Export plot and specify dimensions
ggsave(
  filename = './VIRCLASP_Proteomics/plot_outputs/CedV_InterPro_Enrichment.pdf',
  plot = interpro_plt,
  width = 7,
  height = 5,
  dpi = 300
)

## Reactome Pathway Enrichment
# 1 - Load in Reactome pathway enrichment data output from STRING database (https://string-db.org/)
reactome <- read_xlsx('./VIRCLASP_Proteomics/reference_datasets/CedV_VIRCLASP_ReactomePA_Summary.xlsx')

# 2 - Convert FDR column to numeric class
reactome$p.adjust <- as.numeric(as.character(reactome$p.adjust))

# 3 - Separate viral- and cellular-related pathways for graphing in separate plots
reactome_viral <- reactome %>% 
  dplyr::filter(Category == 'Viral')
reactome_cellular <- reactome %>% 
  dplyr::filter(Category == 'Cellular')

# 4 - Plot InterPro enrichment data on a lollipop plot for both cellular and viral processes
reactome_plt_viral <- ggplot(reactome_viral, aes(x = reorder(Description, -p.adjust), y = Category)) +
  geom_point(aes(size = Count, colour = p.adjust)) +
  scale_size(range = c(2,8)) +
  scale_color_gradient(low = "#9CC184FF", high = "#1E3D14FF") +
  xlab('') +
  ylab("Strength") +
  coord_flip() +
  theme_minimal()
reactome_plt_viral

reactome_plt_cellular <- ggplot(reactome_cellular, aes(x = reorder(Description, -p.adjust), y = Category)) +
  geom_point(aes(size = Count, colour = p.adjust)) +
  scale_size(range = c(1,4)) +
  scale_color_gradient(low = "#9CC184FF", high = "#1E3D14FF") +
  xlab('') +
  ylab("Strength") +
  coord_flip() +
  theme_minimal()
reactome_plt_cellular

# 5 - Export plots and specify dimensions
ggsave(
  filename = './VIRCLASP_Proteomics/plot_outputs/CedV_Reactome_Viral.pdf',
  plot = reactome_plt_viral,
  width = 6,
  height = 7,
  dpi = 300
)

ggsave(
  filename = './VIRCLASP_Proteomics/plot_outputs/CedV_Reactome_Cellular.pdf',
  plot = reactome_plt_cellular,
  width = 8,
  height = 5,
  dpi = 300
)
