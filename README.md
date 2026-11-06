# CEDV VIR-CLASP Proteomics & DRS Analysis

This is the collection of scripts used to analyze VIR-CLASP MS-based proteomics and Nanopore Direct RNA-sequencing data for the CedV manuscript (Ascura et al., Cell Rep, 2026). For VIR-CLASP Proteomics, raw MS data is processed using MaxQuant, which are then analyzed using indicated custom R scripts in (**R Project: CedV_VIRCLASP_DRS.Rproj**). For Nanopore DRS, raw POD5 files are first processed using indicated shell scripts followed by custom R scripts. Scripts are organized by experiment with the relevant manuscript figure number, script description, required inputs, and expected outputs included. Scripts should be executed in numerical order within each section to produce appropriate outputs that may be required for subsequent analyses. For optimal analysis, ensure entire Rproj and contents are downloaded together, maintaining folder structure.

## CedV VIR-CLASP MS-based Proteomics

##### 1) CedV_VIRCLASP_MSstats.R

*Description*: Uses MaxQuant outputs to perform statistical analysis of VIR-CLASP proteins\
*Required inputs*: MaxQuant outputs (evidence.txt and proteinGroups.txt) and Annotation.csv

- Note: Annotation.csv is a custom, user-defined file of MS raw filename (Column 1), condition (Column 2), replicate (Column 3), and isotype label type (Column 4), where L is used for untagged/'light' normal masses. Refer to example in maxquant_output folder.

*Expected output(s)*: File with proteins from CLASP_4SU/Input_4SU comparison (Unfiltered_CedV_VIRCLASP_Hits_AA.csv)

##### 2) CedV_VIRCLASP_Protein_Filtering.R

*Description*: Filter and format raw csv file containing CLASP_4SU/Input_4SU proteins output from CedV_VIRCLASP_MSstats.R\
*Required inputs*: MSstat output (Unfiltered_CedV_VIRCLASP_Hits_AA.csv) and MaxQuant output (proteinGroups.txt)\
*Expected output(s)*: List of VIRCLASP proteins: Human Sp (CedV_VIRCLASP_hits_human_sp_AA.csv/CedV_VIRCLASP_hits_human_AA.csv), all species Sp+Tr (CedV_VIRCLASP_hits_sp_tr_AA.csv)\
*Additional Notes*: CedV_VIRCLASP_hits_human_AA.csv contains the list of **1436 distinct human proteins** identified in the CLASP sample where the **147 CedV RNA interactome proteins ("CedV VIRCLASP proteins")** are filtered from.

##### 3) VIRCLASP_Correlation_Plots.R

*Manuscript Figure*: Figure 1\
*Description*: Assess correlation of VIRCLASP CLASP and Input samples at peptide intensity level\
*Required inputs*: MSstat output (mq_format.csv) and Annotation.csv\
*Expected output(s)*: CLASP vs Input correlation matrix, individual correlation scatterplots (for CLASP or Input)

##### 4) VIRCLASP_Heatmaps.R

*Manuscript Figure*: Supplementary Figure 1\
*Description*: Assess clustering of VIRCLASP CLASP and Input protein intensities between replicates and sample types\
*Required inputs*: VIRCLASP hits (CedV_VIRCLASP_hits_human_sp_AA.csv) and formatted mq_process.csv output from CedV_VIRCLASP_MSstats.R (mq_processed_proteinleveldata.csv)\
*Expected output(s)*: Clustered heatmaps for CLASP replicates 1-4, Input replicates 1-4, and CLASP vs Input sample comparison

##### 5) VIRCLASP_Volcano.R

*Manuscript Figure*: Figure 1\
*Description*: Plot depleted and enriched CLASP proteins using volcano plot\
*Required inputs*: VIRCLASP hits (CedV_VIRCLASP_hits_sp_tr_AA.csv)\
*Expected output(s)*: Volcano plot of depleted and enriched CLASP proteins (VIRCLASP_Volcano.pdf)

##### 6) CedV_VIRCLASP_PutativePreviousRBPs.R

*Manuscript Figure*: Supplementary Figure 1\
*Description*: Determine proportion of CedV RNA interactome that are putative RBPs vs previously identified RBPs based on other vRNA interactomes\
*Required inputs*: VIRCLASP list (CedV_VIRCLASP_hits_human_AA) and compiled dataset reference (Castello_MolCell_2024.xlsx)\
*Expected output(s)*: Barplot of previously identified vs putative RBPs in CedV RNA interactome (Cedv_PutativePreviousRBP_Barplot.pdf)

##### 7) VIRCLASP_ISG_NonISG_Proteins.R

*Manuscript Figure*: Figure 2\
*Description*: Determine enrichment of proteins annotated as ISGs in all CLASP protein sample and CedV RNA interactome\
*Required inputs*: VIRCLASP list (CedV_VIRCLASP_hits_human_AA), IFNome query results (IFNome_CLASP_GeneSearchResults.txt,IFNome_CedVInteractome_GeneSearchResults.txt), and CLASP protein names (CLASP_Proteins.csv)\
*Expected output(s)*: CDF of ISGs vs NonISG candidates identified in all CLASP proteins (CLASP_ISGs_CDF.pdf) or CedV RNA interactome (CedV_Interactome_ISGs_CDF.pdf)

##### 8) VIRCLASP_Functional_Enrichment_Analyses.R

*Manuscript Figure*: Figure 2\
*Description*: Create custom plots using output from various function enrichment analyses (InterPro Features, Reactome)\
*Required inputs*:InterPro Enrichment results (CedV_STRING_InterPro_ProteinDomains.tsv) and Reactome enrichment results (CedV_VIRCLASP_ReactomePA_Summary.xlsx) output from STRING database server\
*Expected output(s)*:Lollipop/Dot plots (CedV_InterPro_Enrichment.pdf,CedV_Reactome_Viral.pdf,CedV_Reactome_Cellular.pdf

------------------------------------------------------------------------

## Nanopore DRS Analysis in SPEN KD/METTL3i (STM2457) vs Control CedV transcripts

*Manuscript Figure*: Figure 6

##### 1) POD5_Merge.sh (Optional)

*Description*: Inspect DRS raw POD5 files and merge files of the same biological replicate and condition\
*Required inputs*: POD5 files\
*Expected outputs*: POD5 files

##### 2) Dorado_Basecalling_m6ADRACH.sh

*Description*: Perform Dorado basecaller to basecall POD5 file into reads with modified base calling for m6A in DRACH context\
*Required inputs*: POD5 files\
*Expected outputs*: BAM files (use only those in bam pass folder)

##### 3) Dorado_Minimap2_Alignment.sh

*Description*: Use Minimap2 to align basecalled reads to CedV and Hg38 reference sequences.\
*Required inputs*: Basecalled BAM files and reference genomes with GTF files\
*Expected outputs*: Aligned BAM files

##### 4) Sort_Index_Aligned_Reads.sh

*Description*: Use Samtools to sort and index aligned reads\
*Required inputs*: Aligned BAM files\
*Expected outputs*: BAM and BAI files

##### 5) NanoCount.sh

*Description*: Compute transcript abundance and TPM (transcripts per million) for each CedV transcript in Control, SPEN, STM2457 samples to generate normalization factor for downstream m6A modification ratios\
*Required inputs*: Sorted and indexed (aligned) BAM files from samtools\
*Expected outputs*: TSV file of TPMs/transcript abundance per transcript in each sample condition

##### 6) Modkit_Modbam_Pileup.sh

*Description*: Tabulate counts of base modifications using Modkit from every sequencing read over each aligned reference  genomic position and apply modification threshold\
*Required inputs*: Sorted and indexed BAM and BAI files\
*Expected outputs*: bedMethyl files

##### 7) Compress_Tabix_Index.sh

*Description*: Compress and tabix index bedMethyl files to prep for differential methylaton scoring\
*Required inputs*: bedMethyl files\
*Expected outputs*: BED and BED.GZ files

##### 8) Modkit_Differential_Methylation_Scoring.sh

*Description*: Perform differential methylation analysis (WT vs SPEN KD/METTL3i) at single base positions with basecalled A (without region segmentation)\
*Required inputs*: Compressed and indexed BED.GZ files for group comparisons (SPEN KD vs Control) and reference genome CedV FASTA\
*Expected outputs*: BED files

##### 9) RDMR_Transcript_Annotation.R

*Description*: Format output from Modkit_Differential_Methylation_Scoring.sh (converted to CSV/XLSX), compute unweighted RDMR for downstream TPM normalization, assign transcript position to start sites, and apply coverage filter. *Required inputs*: Methylation scoring metrics of modified positions in CSV/XLSX format (converted from BED) (cedv_bmt0.7_dmr_control_spen.xlsx/cedv_bmt0.7_dmr_control_stm.xlsx)\
*Expected outputs*: Filtered and transcript-annotated modified positions (.CSV)

##### 10) NanoCount_TPM.R

*Description*: NanoCount output, calculate average TPM ratio, and normalize calculated RDMR for SPEN and STM samples by TPM calculated by Nanocount\
*Required inputs*: Filtered and transcript-annotated modified positions (.CSV)\
*Expected outputs*: TPM-weighted RDMR at modified positions (.CSV)

##### 11) TPM_Normalized_RDMR_Add_Stats.R

*Description*: Calculate TPM-weighted RDMR (aka TPM-weighted signal) and compute FDR for each differentially modified region (single nucleotide position) to prioritize significant sites\
*Required inputs*: Filtered and transcript-annotated modified positions (.CSV) from Step 10\
*Expected outputs*: RDMR and statistics (FDR,etc.) for each nucleotide position (.CSV)

##### 12) TPM_Normalized_RDMR_ModRate_Differences_Individual_Transcripts.R

*Description*: Visualize RDMR with conditional means on individual transcripts\
*Required inputs*: RDMR and statistics (FDR,etc.) for each nucleotide position (.CSV)\
*Expected outputs*: Plot of significant RDMRs at individual nucleotide positions on CedV transcripts

##### 13) TPM_Normalized_RDMR_m6A_Location.R

*Description*: Evaluate % of TPM-normalized RDMR sites that are located in each transcript with RDMR thresholds, related to Table S4\
*Required inputs*: RDMR and statistics (FDR,etc.) for each nucleotide position from Step 11 (.CSV)\
*Expected outputs*: Table of RDMRs (% and Raw Number) located per transcript with RDMR thresholds (.CSV)

##### 14) M6A_A_transcript_ratio.R

*Description*: Calculate m6A-to-A conversion rate per transcript (Horizontal m6A rate per transcript A's), related to Table S4\
*Required inputs*: RDMR and statistics (FDR,etc.) for each nucleotide position from Step 11 (.CSV)\
*Expected outputs*: Table of m6A-to-A conversion rate (% and Raw Number) located per transcript with RDMRs (.CSV)

##### 15) Absolute_ModRate_Metagene.R

*Description*: Visualize absolute modification counts (no TPM_normalization) in Control and SPEN KD/METTL3i CedV transcripts as a metagene-like plot, related to Figure S6H\
*Required inputs*: RDMR and statistics (FDR,etc.) for each nucleotide position from Step 11 (.CSV) - this script will use unnormalized metrics\
*Expected outputs*: Metagene-like plot of absolute modification counts for all conditions(.CSV)

##### 16) Modkit_Differential_Methylation_Scoring_IndividualReplicates.sh

*Description*: Perform differential methylation analysis (WT vs SPEN KD) between biological replicates and sample types at single base positions with basecalled A (without region segmentation) for downstream correlation\
*Required inputs*: Compressed and indexed BED.GZ files for group comparisons (SPEN KD/METTL3i vs Control) and reference genome CedV FASTA\
*Expected outputs*: BED files

##### 17) ModRate_Correlation_Plots.R

*Description*: Assess and plot correlation for SPEN KD/METTL3i and Control Replicates for RDMR from Modkit\
*Required inputs*: CSV/XLSX converted from Step 16 BED files\
*Expected outputs*: Correlation scatterplots
