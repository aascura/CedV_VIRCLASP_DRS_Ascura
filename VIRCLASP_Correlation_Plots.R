# Load in packages
library(ggplot2)
library(ggpubr)
library(dplyr)
library(corrplot)

# Title: VIRCLASP_Correlation_Plots.R
# Author: Andrea Ascura
# Last Update: 18 November 2025

# Description: Assess correlation of VIRCLASP CLASP and Input samples at peptide intensity level

# 1 - Load in mq_format.csv (peptide intensity level data) output from CedV_VIRCLASP_MSstats.R 
mq_format <- read.csv('./VIRCLASP_Proteomics/msstats_output/mq_format.csv', header = TRUE)

# 2 - Read in Annotation.csv file used in CedV_VIRCLASP_MSstats.R
annot <- read.csv("./VIRCLASP_Proteomics/maxquant_output/Annotation.csv", header=TRUE)

# 3 - Filter and format dataframe containing peptide intensities
raw_data <- as.data.frame(mq_format) %>%  
  dplyr::rename('Log2Intensity' = 'Intensity', 'Raw.file' = 'Run') %>% 
  inner_join(annot, by = 'Raw.file') %>%
  dplyr::select(ProteinName,Condition.x,BioReplicate.x,Raw.file,Log2Intensity,Run) %>%
  dplyr::rename('Condition' = 'Condition.x', 'BioReplicate' = 'BioReplicate.x') %>%
  relocate(Run, .after = BioReplicate)

# 4 - Calculate Log2 Intensities 
raw_data$Log2Intensity <- log2(raw_data$Log2Intensity)

# 5 - Create custom function to extract intensities from specific condition and run
rep_filter <- function(Dataset, Condition_Type, BioRep, Run_No){
  Dataset %>%
    dplyr::filter(Condition == Condition_Type & BioReplicate == BioRep & Run == Run_No) %>%
    dplyr::select(Log2Intensity) 
}
# 6 - Extract intensities from Input_4SU replicates: 
Input_4SU_R1_5 <- rep_filter(raw_data,'Input_4SU','R1',5) # Rep 1
Input_4SU_R1_6 <- rep_filter(raw_data,'Input_4SU','R1',6) # Rep 2
Input_4SU_R2_7 <- rep_filter(raw_data,'Input_4SU','R2',7) # Rep 3
Input_4SU_R2_8 <- rep_filter(raw_data,'Input_4SU','R2',8) # Rep 4

# 7 - Combine Input_4SU replicate intensities for two example replicates
Input_Reps_5_7 <- cbind(Input_4SU_R1_5,Input_4SU_R2_7) %>%
  na.omit()
colnames(Input_Reps_5_7) <- c('Log2Intensity_Input4SU_R1_5','Log2Intensity_Input4SU_R2_7')

# 8 - Visualize scatterplot with correlation coefficients for Input_4SU replicates
Input_Reps_Corplt <- ggscatter(Input_Reps_5_7, x = "Log2Intensity_Input4SU_R1_5", y = "Log2Intensity_Input4SU_R2_7",
                               add = "reg.line", 
                               add.params = list(color = '#1D4497FF'),
                               palette = 'npg',
                               xlab = 'Log2 Intensity Replicate 1 (+4SU Input)', ylab = 'Log2 Intensity Replicate 3 (+4SU Input)',
                               font.x = c(12,'bold'), font.y = c(12,'bold'),
                               xlim = c(15,35), ylim = c(15,35))+ stat_cor(method = "pearson", label.x = 30, label.y = 15)

Input_Pt_Num <- paste('N = ', nrow(Input_Reps_5_7))
ggpar(Input_Reps_Corplt, 
      caption = Input_Pt_Num)

# 9 - Extract intensities from CLASP_4SU replicates:
CLASP_4SU_R1_13 <- rep_filter(raw_data,'CLASP_4SU','R1',13) # Rep 1
CLASP_4SU_R1_14 <- rep_filter(raw_data,'CLASP_4SU','R1',14) # Rep 2
CLASP_4SU_R2_15 <- rep_filter(raw_data,'CLASP_4SU','R2',15) # Rep 3
CLASP_4SU_R2_16 <- rep_filter(raw_data,'CLASP_4SU','R2',16) # Rep 4

# 10 - Combine CLASP_4SU replicate intensities for two example replicates
CLASP_Reps_13_15 <- cbind(CLASP_4SU_R1_13,CLASP_4SU_R2_15) %>%
  na.omit()
colnames(CLASP_Reps_13_15) <- c('Log2Intensity_CLASP4SU_R1_13','Log2Intensity_CLASP4SU_R2_15')

# 11 - Visualize scatterplot with correlation coefficients for CLASP_4SU replicates
CLASP_Reps_Corplt <- ggscatter(CLASP_Reps_13_15, x = "Log2Intensity_CLASP4SU_R1_13", y = "Log2Intensity_CLASP4SU_R2_15",
                               add = "reg.line", 
                               add.params = list(color = '#1D4497FF'),
                               palette = 'npg',
                               xlab = 'Log2 Intensity Replicate 1 (+4SU CLASP)', ylab = 'Log2 Intensity Replicate 3 (+4SU CLASP)',
                               font.x = c(12,'bold'), font.y = c(12,'bold'),
                               xlim = c(15,30), ylim = c(15,30)
)+ stat_cor(method = "pearson", label.x = 27, label.y = 15)

CLASP_Reps_Corplt
CLASP_Pt_Num <- paste('N = ', nrow(CLASP_Reps_13_15))
ggpar(CLASP_Reps_Corplt, 
      caption = CLASP_Pt_Num)

# 12 - Setup up combined dataset to evaluate correlation for ALL Input and CLASP replicates 1-4
Input_CLASP_Reps_All <- cbind(Input_4SU_R1_5,
                              Input_4SU_R1_6,
                              Input_4SU_R2_7,
                              Input_4SU_R2_8,
                              CLASP_4SU_R1_13,
                              CLASP_4SU_R1_14,
                              CLASP_4SU_R2_15,
                              CLASP_4SU_R2_16) %>% na.omit()
colnames(Input_CLASP_Reps_All) <- c('Input Rep. 1',
                                    'Rep. 2',
                                    'Rep. 3',
                                    'Rep. 4',
                                    'CLASP Rep. 1',
                                    'Rep. 2',
                                    'Rep. 3',
                                    'Rep. 4')

# 13 - Visualize full correlation matrix of Input and CLASP replicates using corrplot
custom_scale <- colorRampPalette(c('#eff1ff','#1D4497FF'))
custom_range <- custom_scale(600)

corrplot(cor(Input_CLASP_Reps_All), # Pearson correlation is used
         method = 'shade',
         type = 'lower',
         cl.pos = 'r',
         is.corr = FALSE,
         col = custom_range,
         col.lim = c(0.50,1),
         bg = 'white',
         outline = TRUE,
         tl.col = 'black',
         addgrid.col = 'black',
         diag = TRUE,
         tl.srt = 0,
         addCoef.col = 'white')