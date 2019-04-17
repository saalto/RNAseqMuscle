# --- Location of Iteration #2 Processed Data ---
home_loc <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/'

# ------------------------------------------------------------------------------------------------------

# --- Soleus (Sol) Control versus Sol Mutant Tissue Comparison ---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.Flz_Sol_S1.txt','sorted.Flz_Sol_S2.txt','sorted.Flz_Sol_S3.txt',
                  'sorted.Flz_Sol_S4.txt','sorted.Flz_Sol_S5.txt','sorted.Flz_Sol_S6.txt',
                  'sorted.Mck_Sol_S1.txt','sorted.Mck_Sol_S2.txt','sorted.Mck_Sol_S3.txt',
                  'sorted.Mck_Sol_S4.txt','sorted.Mck_Sol_S5.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'soleus-data-combined.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
#data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/soleus-data-combined.txt'
#data_file = 'C:/Users/sarah/Documents/AllSoleusReadCountsWO6.txt'
data_file = 'C:/Users/sarah/OneDrive/Documents/2019/01_2019_Winter/PROJECT/RNASeq_Analysis/2019_03_15/AllSoleus_ReadCounts.txt'

# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','wt 3','wt 4','wt 5','wt 6','mut 1','mut 2','mut 3','mut 4','mut 5')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,1,1,1,2,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and mut = treated here)
condition <- c(rep("untreated",6),rep("treated",5))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 -12

deresSOL <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgSOL_flt <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

# ------------------------------------------------------------------------------------------------------

# --- Soleus (Sol) Control versus Sol Mutant Tissue Comparison WITH SAMPLE C6 REMOVED---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.Flz_Sol_S1.txt','sorted.Flz_Sol_S2.txt','sorted.Flz_Sol_S3.txt',
                  'sorted.Flz_Sol_S4.txt','sorted.Flz_Sol_S5.txt',
                  'sorted.Mck_Sol_S1.txt','sorted.Mck_Sol_S2.txt','sorted.Mck_Sol_S3.txt',
                  'sorted.Mck_Sol_S4.txt','sorted.Mck_Sol_S5.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'AllSoleusWO6_ReadCounts.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
#data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/soleus-data-combined.txt'
#data_file = 'C:/Users/sarah/Documents/AllSoleusReadCountsWO6.txt'
data_file = 'C:/Users/sarah/OneDrive/Documents/2019/01_2019_Winter/PROJECT/RNASeq_Analysis/2019_03_15/AllSoleus_ReadCounts.txt'

# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','wt 3','wt 4','wt 5','mut 1','mut 2','mut 3','mut 4','mut 5')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,1,1,2,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and mut = treated here)
condition <- c(rep("untreated",5),rep("treated",5))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 -12

deresSOLWO6 <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgSOLWO6_flt <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

#--------------------------------------------------------------------------------------------------------------------------

# --- Tibialis Anterior (TA) Control versus TA Mutant Tissue Comparison WITHOUT REMOVAL AND RENAMING FILES---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.TF1_Ta_Flz16.txt','sorted.TF2_Ta_M11.txt','sorted.TF3_Ta_Flz16.txt',
                  'sorted.TM1_Ta_M2.txt','sorted.TM2_Ta_M13.txt','sorted.TM3_Ta_Flz13.txt',
                  'sorted.TM4_Ta_M12.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'AllTibialis_ReadCounts.txt'

# All 7 samples are present in the following file
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/tibialis-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','wt 3', 'mut 1','mut 2','mut 3','mut 4')
# set up the two groups (in this case, I have 2 wt replicates followed by 4 mut replicates)
groups <- factor(c(1,1,1,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and mut = treated here)
condition <- c(rep("untreated",3),rep("treated",4))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 - 8

deresTA <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgTA_flt <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

#---------------------------------------------------------------------------------------------------------------------------

# --- Tibialis Anterior (TA) Control versus TA Mutant Tissue Comparison WITH RENAMED FILES---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/renamed_files'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.TF1_Ta_Flz16.txt','sorted.TF2_Ta_Flz16.txt','sorted.TM1_Ta_M11.txt',
                  'sorted.TM2_Ta_M2.txt','sorted.TM3_Ta_M13.txt','sorted.TM4_Ta_Flz13.txt','sorted.TM5_Ta_M12.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'AllTibialis_renamed_ReadCounts.txt'

# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# data_file - the text file of the combined read counts
data_file <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/renamed_files/tibialis-data-combined1.txt'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','mut 1','mut 2','mut 3','mut 4', 'mut 5')
# set up the two groups (in this case, I have 2 wt replicates followed by 4 mut replicates)
groups <- factor(c(1,1,2,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and mut = treated here)
condition <- c(rep("untreated",2),rep("treated",5))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 - 8

deresTA_renamed <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgTA_flt_renamed <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

#--------------------------------------------------------------------------------------------------------------

# --- Tibialis Anterior (TA) Control versus TA Mutant Tissue Comparison WITH RENAMED AND REMOVED FILES---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
# name for the combined file (it is saved to data_path)
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/renamed_files_2'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.TF1_Ta_Flz16.txt','sorted.TF2_Ta_Flz16.txt','sorted.TM1_Ta_M11.txt',
                  'sorted.TM2_Ta_M2.txt','sorted.TM4_Ta_Flz13.txt','sorted.TM5_Ta_M12.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'AllTibialis_renamedremoved_ReadCounts.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/tibialis-data-combined2.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','mut 1','mut 2','mut 3','mut 4')
# set up the two groups (in this case, I have 2 wt replicates followed by 4 mut replicates)
groups <- factor(c(1,1,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and mut = treated here)
condition <- c(rep("untreated",2),rep("treated",4))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 - 8

deresTA_renamedremoved <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgTA_flt_renamedremoved <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

#---------------------------------------------------------------------------------------------------------------------------------

# --- Gastrocnemius (Ga) Control versus Ga Mutant Tissue Comparison WITHOUT REMOVAL AND RENAMING---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/gastrocnemius'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.GF1_ga_flz13.txt','sorted.GF2_ga_flz14.txt',
                  'sorted.GM1_ga_m11.txt','sorted.GM2_ga_m12.txt',
                  'sorted.GM3_ga_m2.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'AllGastro_ReadCounts.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/gastro-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','mut 1','mut 2','mut 3')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("untreated",2),rep("treated",3))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 -6

deresGA <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgGA_flt <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

#----------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------

# --- Gastrocnemius (Ga) Control versus Ga Mutant Tissue Comparison WITH RENAMING---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/gastrocnemius/renamed_files'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.GF1_ga_flz13.txt','sorted.GF2_ga_flz14.txt',
                  'sorted.GF3_ga_m11.txt','sorted.GM2_ga_m12.txt',
                  'sorted.GM3_ga_m2.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'AllGastro_renamed_ReadCounts.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/gastro-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','mut 1','mut 2','mut 3')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("untreated",3),rep("treated",2))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 -6

deresGA <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgGA_flt_renamed <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

#--------------------------------------------------------------------------------------------------------------
# --- Gastrocnemius (Ga) Control versus Ga Mutant Tissue Comparison WITH REMOVAL AND RENAMING---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/gastrocnemius/renamed_files2'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.GF1_ga_flz13.txt','sorted.GF2_ga_m11.txt',
                  'sorted.GM2_ga_m12.txt','sorted.GM3_ga_m2.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'AllGastro_renamedremoved_ReadCounts.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/gastro-data-combined2.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','mut 1','mut 2')
# set up the two groups (in this case, I have 6 wt  replicates followed by 5 mut replicates)
groups <- factor(c(1,1,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("untreated",2),rep("treated",2))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 - 5

deresGA <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgGA_flt_renamedremoved  <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

# -----------------------------------------------------------------------
# --- Combining non-oxidative myofiber types TA and Gastro ---
# --- Tibialis Anterior (TA) Control versus TA Mutant Tissue Comparison WITH RENAMED FILES---
# --- Gastrocnemius (Ga) Control versus Ga Mutant Tissue Comparison WITH REMOVAL AND RENAMING---

# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/combined'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.TF1_Ta_Flz16.txt','sorted.TF2_Ta_Flz16.txt','sorted.GF1_ga_flz13.txt',
                  'sorted.GF2_ga_m11.txt','sorted.TM1_Ta_M11.txt','sorted.TM2_Ta_M2.txt',
                  'sorted.TM3_Ta_M13.txt','sorted.TM4_Ta_Flz13.txt','sorted.TM5_Ta_M12.txt', 
                  'sorted.GM2_ga_m12.txt','sorted.GM3_ga_m2.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'TibialisGastro_combined_ReadCounts.txt'


# directory to save results - no trailing forward slash!
#save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_13'
# data_file - the text file of the combined read counts
data_file <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/combined/TibialisGastro_combined_ReadCounts.txt'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wtTA 1','wtTA 2','wtGA 3','wtGA 4','mutTA 1','mutTA 2','mutTA 3','mutTA 4', 'mutTA 5', 'mutGA 6', 'mutGA 7')
# set up the two groups (in this case, I have 2 wt replicates followed by 4 mut replicates)
groups <- factor(c(1,1,1,1,2,2,2,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and mut = treated here)
condition <- c(rep("untreated",4),rep("treated",7))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

# DESeq2_DE.R and edgeR_DE.R functions will need to modified to make the correct dimensions for the data.frame ---
# Number of columns 2 - 8
deresTAGA_combined <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedgTAGA_combined <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)
