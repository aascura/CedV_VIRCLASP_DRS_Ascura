# Load in packages
library(MSstats)
library(dplyr)

# Title: CedV_VIRCLASP_MSstats.R
# Author: Andrea Ascura
# Last Update: 18 November 2025

# Description: Statistical analysis of VIRCLASP proteins from MaxQuant identification 

# 1 - Read in proteinGroups.txt file (MaxQuant output file)
proteinGroups <- read.table("./VIRCLASP_Proteomics/maxquant_output/proteinGroups.txt", 
                            sep="\t", header=TRUE) # Check original proteinGroups output file - changed from sep="\t" to sep=""

# 2 - Read in evidence.txt file (MaxQuant output file)
evi <- read.table("./VIRCLASP_Proteomics/maxquant_output/evidence.txt", sep="\t", header=TRUE)

# 3 - Read in Annotation.csv (User-created file including raw file names, condition, biological replicates, and labeling type)
# Note: Refer to file for required formatting
annot <- read.csv("./VIRCLASP_Proteomics/maxquant_output/Annotation.csv", header=TRUE)

# 4 - Remove keratin from proteinGroups
mq_proteinGroups <- proteinGroups[!grepl("keratin",
                                         proteinGroups$Fasta.headers,
                                         ignore.case = TRUE), ] 

# 5 - Keep all potential contaminants (may be RBPs)
for (i in 1:nrow(mq_proteinGroups)) {
  if (mq_proteinGroups$Potential.contaminant[i] == "+" || is.na(mq_proteinGroups$Potential.contaminant[i])) {
    mq_proteinGroups$Potential.contaminant[i] <- ""
  }
}

# 6 - Convert MaxQuant output into MSstats format - retrieves peptide ion intensity from evidence.txt
mq_format <- MaxQtoMSstatsFormat(evidence = evi,
                                 annotation = annot,
                                 proteinGroups = mq_proteinGroups,
                                 proteinID = "Leading.razor.protein",
                                 useUniquePeptide = TRUE,
                                 summaryforMultipleRows = max,
                                 removeFewMeasurements = TRUE,
                                 removeMpeptides = FALSE,
                                 removeOxidationMpeptides = FALSE,
                                 removeProtein_with1Peptide = FALSE)

# 7 - Process data to performm normalization and summarization
mq_processed <- dataProcess(mq_format,
                            logTrans = 2,
                            normalization = 'QUANTILE',
                            featureSubset = "all",
                            summaryMethod = "TMP",
                            remove50missing = FALSE,
                            censoredInt = "NA", # IMPORTANT FOR MAXQUANT!!
                            MBimpute = TRUE,
                            fix_missing = FALSE,
                            maxQuantileforCensored = 0.999)

# 8 - Visualize effects of normalization between all proteins - should be stable
dataProcessPlots(data = mq_processed,
                 which.Protein = "allonly",
                 type="QCplot",
                 ylimUp=35,
                 width=5,
                 height=5)

# 9 - Group Comparison: Identify differential enriched proteins
# 9A - Previous condition order for matrix: CLASP_4SU, CLASP_no4SU, Input_4SU, Input_no4SU
levels(mq_processed$ProteinLevelData$GROUP)

# 9B - Setup comparison for CLASP_4SU / Input_4SU based on order returned from Step 9A.
contrast.matrix <- NULL
comparison <- matrix(c(1, 0, -1, 0), nrow = 1)
contrast.matrix <- rbind(contrast.matrix, comparison)
row.names(contrast.matrix) <- c('CLASP_4SU vs. Input_4SU')
colnames(contrast.matrix) <- c('CLASP_4SU',
                               'CLASP_no4SU',
                               'Input_4SU',
                               'Input_no4SU')

mq_comparison <- groupComparison(contrast.matrix, data = mq_processed)

# 10 - Visually confirm presence of differential abundant proteins: CLASP_4SU/Input_4SU
groupComparisonPlots(data = mq_comparison$ComparisonResult,
                     type = 'VolcanoPlot',
                     width = 5,
                     height = 5,
                     address = 'mq_comparison_plot')

# 11 - Export file containing proteins from CLASP_4SU/Input_4SU comparison
# This is the input file required for CedV_VIRCLASP_Interactome.R
# Note: This file is the raw, unfiltered groupComparison of log2(CLASP_4SU/Input_4SU) with nonnumeric pvalues, all species, and isoforms
write.csv(mq_comparison$ComparisonResult,'./VIRCLASP_Proteomics/msstats_output/Unfiltered_CedV_VIRCLASP_Hits_AA.csv', row.names = FALSE)

# 12 - Export protein intensities for all samples 
# This is the input file required for CedV_VIRCLASP_Correlation_Plots.R
write.csv(mq_format, './VIRCLASP_Proteomics/msstats_output/mq_format.csv', row.names = FALSE)

# 13 - Export protein level data for all samples 
# This is the input file required for CedV_VIRCLASP_Heatmap.R
write.csv(mq_processed$ProteinLevelData, './VIRCLASP_Proteomics/msstats_output/mq_processed_proteinleveldata.csv', row.names = FALSE)
