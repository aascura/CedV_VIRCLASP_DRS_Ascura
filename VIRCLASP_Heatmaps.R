# Load in packages
library(dplyr)
library(data.table)
library(pheatmap)
library(tidyverse)
library(dendsort)

# Title: VIRCLASP_Heatmaps.R
# Author: Andrea Ascura
# Last Update: 18 November 2025

# Description: Assess clustering of VIRCLASP CLASP and Input protein intensities

# 1 - Load in CedV VIRCLASP Human Hits output from CedV_VIRCLASP_Interactome.R
cedv_rbps_hsa <- fread('./VIRCLASP_Proteomics/CedV_VIRCLASP_hits_human_sp_AA.csv')

# 2 - Load in MSstats output file, mq_process, from CedV_VIRCLASP_MSstats.R containing log2intensities
mq_processed <- fread('./VIRCLASP_Proteomics/msstats_output/mq_processed_proteinleveldata.csv')

# 3 - Extract log intensities from mq_processed based on inclusion in CedV Human Hits
cedv_hsa <- mq_processed %>%
  dplyr::filter(Protein %in% cedv_rbps_hsa$Uniprot_ID)

# 4 - Create custom function to extract protein intensities from specific condition and run
get_run <- function(Run_No, Group_Name, Subject_No){
  df <- cedv_hsa %>%
    dplyr::filter(Run_No == RUN & Group_Name == GROUP & Subject_No == SUBJECT) %>%
    dplyr::select(Protein, LogIntensities)
}

### Heatmap Clustering and Visualization: Input Only
# 5 - Extract log intensities of Input_4SU replicates and combine all Input replicates
Input_list <- list(get_run(9,'Input_4SU','R1'), # Subject = R1, Run = 9,   Rep 1
                   get_run(10,'Input_4SU','R1'), # Subject = R1, Run = 10, Rep 2
                   get_run(11,'Input_4SU','R2'), # Subject = R2, Run = 11, Rep 3
                   get_run(12,'Input_4SU','R2') # Subject = R2, Run = 12 # Rep 4
)

Input_AllReps <- Input_list %>%
  reduce(inner_join, by = 'Protein')

Input_combined <- as.data.frame(do.call(rbind, Input_list)) %>%
  dplyr::select(Protein)
Input_count <- table(Input_combined$Protein) %>% as.data.frame()
Input_proteins <- merge(Input_count, cedv_rbps_hsa, by.x = 'Var1', by.y = 'Uniprot_ID') %>%
  dplyr::select(Primary_gene_name, log2FC, Var1)

# 6 - Reformat combined Input data frame 
rownames(Input_AllReps) <- Input_AllReps$Protein
Input_AllReps <- Input_AllReps %>%
  dplyr::select(-Protein) %>%
  dplyr::rename('Replicate 1' = 'LogIntensities.x',
                'Replicate 2' = 'LogIntensities.y',
                'Replicate 3' = 'LogIntensities.x.x',
                'Replicate 4' = 'LogIntensities.y.y')

# 7 - Split columns into 4 groups and add column annotations
col_groups_input <- substr(colnames(Input_AllReps), 11, 11)
table(col_groups_input)

input_col <- data.frame(Replicate = col_groups_input)
rownames(input_col) <- colnames(Input_AllReps)

# 8 - Sort the dendrogram
sort_hclust <- function(x) as.hclust(dendsort(as.dendrogram(x)))
input_cluster_rows <- sort_hclust(hclust(dist(Input_AllReps)))
input_cluster_cols <- sort_hclust(hclust(dist(t(Input_AllReps))))

# 9 - Scale and visualize Input (All Replicates) using pheatmap functions
scaled_data_input_allreps <- scale(Input_AllReps)
clustered_data_input_allreps <- pheatmap(scaled_data_input_allreps,
                                         cluster_rows = input_cluster_rows,
                                         cluster_cols = input_cluster_cols,
                                         annotation_col = input_col, # Turn on for spacing purposes. Color bars will be manually replaced with text labels
                                         cutree_rows = 5,
                                         cellwidth = 35,
                                         cellheight = NA,
                                         color = colorRampPalette(c("#1D4497FF",'#f1f1f1','#B83326FF')) (99), # 100 color divisions results in white breaks 
                                         show_colnames = TRUE,
                                         breaks = seq(-4,4, length.out = 100),
                                         angle_col = '0',
                                         main = 'Input')

### Heatmap Clustering and Visualization: CLASP Only
# 10 -  Extract logIntensities of CLASP_4SU Replicates into individual data frames and combine all CLASP replicates
CLASP_list <- list(get_run(1,'CLASP_4SU','R1'), # Subject = R1, Run = 1, Rep 1
                   get_run(2,'CLASP_4SU','R1'), # Subject = R1, Run = 2, Rep 2
                   get_run(3,'CLASP_4SU','R2'), # Subject = R2, Run = 3, Rep 3
                   get_run(4,'CLASP_4SU','R2') # Subject = R2, Run = 4, Rep 4
) 

CLASP_AllReps <- CLASP_list %>%
  reduce(inner_join, by = 'Protein')

# 11 - Reformat combined CLASP data frame 
rownames(CLASP_AllReps) <- CLASP_AllReps$Protein
CLASP_AllReps <- CLASP_AllReps %>%
  dplyr::select(-Protein) %>%
  dplyr::rename('Replicate 1' = 'LogIntensities.x',
                'Replicate 2' = 'LogIntensities.y',
                'Replicate 3' = 'LogIntensities.x.x',
                'Replicate 4' = 'LogIntensities.y.y')

# 12 - Split columns into 4 groups and add column annotations
col_groups_clasp <- substr(colnames(CLASP_AllReps), 11, 11)
table(col_groups_clasp)

clasp_col <- data.frame(Replicate = col_groups_clasp)
rownames(clasp_col) <- colnames(CLASP_AllReps)

# Write out CLASP proteins
#write.csv(CLASP_AllReps,'CLASP_Proteins.csv')

# 13 - Sorting the dendrogram
clasp_cluster_rows <- sort_hclust(hclust(dist(CLASP_AllReps)))
clasp_cluster_cols <- sort_hclust(hclust(dist(t(CLASP_AllReps))))
clasp_cluster_cols$order <- c(1,2,3,4) # Override ordering to maintain 1-4 order

# 14 - Scale and visualize Input (All Replicates) using pheatmap functions
scaled_data_CLASP_allreps <- scale(CLASP_AllReps)

clustered_data_CLASP_allreps <- pheatmap(scaled_data_CLASP_allreps,
                                         cluster_rows = clasp_cluster_rows,
                                         cluster_cols = clasp_cluster_cols, # Dendrogram will be manually added 
                                         annotation_col = clasp_col, # Turn on for spacing purposes. Color bars will be replaced with text labels
                                         cutree_rows = 5,
                                         cellwidth = 35,
                                         cellheight = NA,
                                         color = colorRampPalette(c("#1D4497FF",'#f1f1f1','#B83326FF')) (99),# 100 color divisions results in white breaks
                                         show_colnames = TRUE,
                                         breaks = seq(-4,4, length.out = 100),
                                         angle_col = '0',
                                         main = 'CLASP')

### Heatmap Clustering and Visualization using log2 intensities AVERAGE for all replicates in same sample type
# 15 - Calculate average for each protein across replicates and reformat data frame
Input_AllReps$Replicate_Average <- rowMeans(Input_AllReps, na.rm = TRUE)
Input_AllReps_Avg <- Input_AllReps %>%
  tibble::rownames_to_column('Protein') %>%
  dplyr::select(Protein,Replicate_Average)

CLASP_AllReps$Replicate_Average <- rowMeans(CLASP_AllReps, na.rm = TRUE)
CLASP_AllReps_Avg <- CLASP_AllReps %>%
  tibble::rownames_to_column('Protein') %>%
  dplyr::select(Protein,Replicate_Average)

# 16 - Combine Input and CLASP data frame and reformat
Input_CLASP_Averages <- inner_join(Input_AllReps_Avg, # x = Input
                                   CLASP_AllReps_Avg, # y = CLASP
                                   by = 'Protein')
rownames(Input_CLASP_Averages) <- Input_CLASP_Averages$Protein
Input_CLASP_Averages <- Input_CLASP_Averages %>%
  dplyr::select(-Protein) %>%
  dplyr::rename('Input' = 'Replicate_Average.x',
                'CLASP' = 'Replicate_Average.y')

# 17 - Split columns into 4 groups and add column annotations
col_groups_Input_CLASP <- substr(colnames(Input_CLASP_Averages), 1, 5)
table(col_groups_Input_CLASP)

Input_CLASP_col <- data.frame(Sample = col_groups_Input_CLASP)
rownames(Input_CLASP_col) <- colnames(Input_CLASP_Averages)

# 18 - Sorting the dendrogram
Input_CLASP_cluster_rows <- sort_hclust(hclust(dist(Input_CLASP_Averages)))

# 19 - Scale() centers the values in a vector by subtracting the mean and dividing by the standard deviation
scaled_Input_CLASP_Average <- scale(Input_CLASP_Averages)

# 20 - Heatmap visualization
cedv_matrix <- as.matrix(scaled_Input_CLASP_Average)
Input_CLASP_Avg_Phmap <- pheatmap(cedv_matrix, 
                                  cluster_cols = F, # Turn off clustering by columns
                                  cluster_rows = Input_CLASP_cluster_rows,
                                  annotation_col = Input_CLASP_col,
                                  cellwidth = 55,
                                  cellheight = NA,
                                  cutree_rows = 5,
                                  color = colorRampPalette(c("#1D4497FF",'#f1f1f1','#B83326FF')) (100),
                                  cutree_cols = 2,
                                  show_colnames = TRUE,
                                  breaks = seq(-4,4, length.out = 101),
                                  angle_col = '0',
                                  main = 'Input and CLASP Replicate Average')

