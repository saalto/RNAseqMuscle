# load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
# load libraries
library("DESeq2")
library("ggplot2")
library("factoextra")

# ---------------------------------------------------------------------------------------
# --- Tibialis Anterior (TA) Control versus TA Mutant Tissue Comparison ---

# --- Use 01_combine-data-files.R script to generate the combined text/Data variables ---
  # location holding the count files - each count file is assumed to be one gene per line,
  #   with a single count value (one integer), tab separated, after the gene symbol name
  data_path <- 'C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior'
  # list of files (filenames) to mush together - the final order depends on the order you
  #   list them in this vector
  files_to_agg <- c('sorted.TF1_Ta_Flz16.txt','sorted.TF2_Ta_M11.txt','sorted.TF3_Ta_Flz16.txt',
                    'sorted.TM1_Ta_M2.txt','sorted.TM2_Ta_M13.txt','sorted.TM3_Ta_Flz13.txt',
                    'sorted.TM4_Ta_M12.txt')
  # number of quality control (non-gene) lines at the end of the file, in case this changes
  n_qc <- 5
  # name for the combined file (it is saved to data_path)
  comb_file <- 'tibialis-data-combined.txt'

# -- (2) Use 01_combine-data-files.R script to generate the combined text/data variables --

# -- (3) Load these variables prior to running 02_edgeR.DE.R --# All 7 samples are present in the following file
  data_file = 'C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior/tibialis-data-combined.txt'
  # set up the two groups (in this case, I have 3 wt replicates followed by 4 mut replicates)
  groups <- factor(c(1,1,1,2,2,2,2))
  # threshold for significance
  adjp <- 0.05
  # drop genes with very low total counts? (recommended)
  drop_low = TRUE
  # changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
  qlt_test = FALSE

# -- (4) Load 02_edgeR_DE.R function

# -- (5) Run the following line to generate differential expression (DE) analysis --
  deedgTA_flt <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

# -- (6) Examine the DE analysis using DESeq2 and PCA --
  # set working directory
  setwd("C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior")
  directory <- "C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior"
  # set up DESeq2 data, based on file names of HTSeq counts in working directory
  sampleFiles <- dir(pattern = 'sorted')
  identifierNames <- c("C1", "C2", "C3", "M1", "M2", "M3", "M4")
  ConditionMatch <- regexpr(pattern = '[A-Z]+', sampleFiles)
  sampleConditions <- regmatches(sampleFiles, ConditionMatch)
  sampleTable <- data.frame(sampleName = identifierNames, fileName = sampleFiles, 
                            condition = sampleConditions)
  # calculate DESeq2 from HTSeq count tables and filter out zero counts
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                         directory = directory, design = ~ condition)
  ddsHTSeqFiltered <- ddsHTSeq [ rowSums(counts(ddsHTSeq)) > 0, ] 
  ddsHTSeqFiltered <- DESeq(ddsHTSeqFiltered)
  # visualize Data Transformations by generating log2 normalized count matrix
  rld <- rlog(ddsHTSeqFiltered, blind = FALSE)
  logTransCounts <- assay(rld)
  # examine principle component analysis (PCA) based on log2 normalized count matrix
  OGPCAN <-prcomp(logTransCounts, center = T, scale = F, tol = 0)
  # Eigenvalues
  eig.val <- get_eigenvalue(OGPCAN)
  OGPCAN_matrix <- as.data.frame(OGPCAN$rotation)
  OGPCAN_matrix$Condition <- c("Control","Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant")
  ggplot(OGPCAN_matrix, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) + geom_text(aes(label=row.names(OGPCAN_matrix)),vjust=0, hjust=0)+
    theme(axis.text.x = element_text(size = 14, color = "black"),
          axis.title.x = element_text(size  = 16, face = "bold"),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 16, face = 'bold'),
          legend.text = element_text(size = 14)) +
    scale_color_discrete(name = "Sample") +
    xlab(paste0("PC1: ", sprintf("%.3f", eig.val$variance.percent[1]), "% variance")) +
    ylab(paste0("PC2: ", sprintf("%.3f", eig.val$variance.percent[2]), "% variance"))

# -- (6) Examine the DE analysis using heatmap of z-scores --
  # Run 07_zscoreGeneration.R, Section 3: Heat map of Tibialis Anterior data 

# -- (7) Continue to section below to examine TA dataset with renaming Control 2 to Mutant 1 --

#---------------------------------------------------------------------------------------------------------------------------
# --- Tibialis Anterior (TA) Control versus TA Mutant Tissue Comparison WITH RENAMED FILES---

# -- (1) Load these variables prior to running 01_combine-data-files.R --
  # location holding the count files - each count file is assumed to be one gene per line,
  #   with a single count value (one integer), tab separated, after the gene symbol name
  data_path <- 'C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior/renamed_files'
  # list of files (filenames) to mush together - the final order depends on the order you
  #   list them in this vector
  files_to_agg <- c('sorted.TF1_Ta_Flz16.txt','sorted.TF2_Ta_Flz16.txt','sorted.TM1_Ta_M11.txt',
                    'sorted.TM2_Ta_M2.txt','sorted.TM3_Ta_M13.txt','sorted.TM4_Ta_Flz13.txt','sorted.TM5_Ta_M12.txt')
  # number of quality control (non-gene) lines at the end of the file, in case this changes
  n_qc <- 5
  # name for the combined file (it is saved to data_path)
  comb_file <- 'AllTibialis_renamed_ReadCounts.txt'

# -- (2) Use 01_combine-data-files.R script to generate the combined text/Data variables --

# -- (3) Load these variables prior to running 02_edgeR.DE.R --
  # data_file - the text file of the combined read counts
  data_file <- 'C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior/renamed_files/AllTibialis_renamed_ReadCounts.txt'
  # set up the two groups (in this case, I have 2 wt replicates followed by 4 mut replicates)
  groups <- factor(c(1,1,2,2,2,2,2))
  # threshold for significance
  adjp <- 0.05
  # drop genes with very low total counts? (recommended)
  drop_low = TRUE
  # changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
  qlt_test = FALSE

# -- (4) Use 02_edgeR_DE.R script to generate differential expression (DE) analysis --
  deedgTA_flt_renamed <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

# -- (5) Examine the DE analysis using DESeq2 and PCA --
  # set working directory
  setwd("C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior/renamed_files")
  directory <- "C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis anterior/renamed_files"
  # set up DESeq2 data, based on file names of HTSeq counts in working directory
  sampleFiles <- dir(pattern = 'sorted')
  identifierNames <- c("C1", "C2", "M1", "M2", "M3", "M4", "M5")
  ConditionMatch <- regexpr(pattern = '[A-Z]+', sampleFiles)
  sampleConditions <- regmatches(sampleFiles, ConditionMatch)
  sampleTable <- data.frame(sampleName = identifierNames, fileName = sampleFiles, 
                            condition = sampleConditions)
  # calculate DESeq2 from HTSeq count tables and filter out zero counts
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                         directory = directory, design = ~ condition)
  ddsHTSeqFiltered <- ddsHTSeq [ rowSums(counts(ddsHTSeq)) > 0, ] 
  ddsHTSeqFiltered <- DESeq(ddsHTSeqFiltered)
  # visualize Data Transformations by generating log2 normalized count matrix
  rld <- rlog(ddsHTSeqFiltered, blind = FALSE)
  logTransCounts <- assay(rld)
  # examine principle component analysis (PCA) based on log2 normalized count matrix
  OGPCAN <-prcomp(logTransCounts, center = T, scale = F, tol = 0)
  # Eigenvalues
  eig.val <- get_eigenvalue(OGPCAN)
  OGPCAN_matrix <- as.data.frame(OGPCAN$rotation)
  OGPCAN_matrix$Condition <- c("Control","Control","Mutant", "Mutant", "Mutant", "Mutant", "Mutant")
  ggplot(OGPCAN_matrix, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) + geom_text(aes(label=row.names(OGPCAN_matrix)),vjust=0, hjust=0)+
    theme(axis.text.x = element_text(size = 14, color = "black"),
          axis.title.x = element_text(size  = 16, face = "bold"),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 16, face = 'bold'),
          legend.text = element_text(size = 14)) +
    scale_color_discrete(name = "Sample") +
    xlab(paste0("PC1: ", sprintf("%.3f", eig.val$variance.percent[1]), "% variance")) +
    ylab(paste0("PC2: ", sprintf("%.3f", eig.val$variance.percent[2]), "% variance"))

# -- (6) Examine the DE analysis using heatmap of z-scores --
  # Run 07_zscoreGeneration.R, Section 4: Heat map of Tibialis Anterior data with the renamed samples