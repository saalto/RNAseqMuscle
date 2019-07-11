# --- Location of Iteration #2 Processed Data ---
home_loc <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/'

#-------------------------------------------------------------------------------------------------------------------------

# --- Soleus versus TA Control Tissue Comparison ---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/control/TAvsSOL'

# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.Ta_Flz_S1.txt','sorted.Ta_Flz_S2.txt','sorted.Sol_Flz_S1.txt','sorted.Sol_Flz_S2.txt','sorted.Sol_Flz_S3.txt',
                  'sorted.Sol_Flz_S4.txt','sorted.Sol_Flz_S5.txt')

# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5

# name for the combined file (it is saved to data_path)
comb_file <- 'soleusandTAC-data-combined.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/control/TAvsSOL/soleusandTAC-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_14/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('TA1', 'TA2','Sol 1','Sol 2','Sol 3','Sol 4','Sol 5')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,2,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("TA",2), rep("Sol",5))
adjp <- 0.05
drop_low <- TRUE
# changing between quasi-likelihood test (TRUE) to likelihood F test (FALSE)
qlt_test = FALSE

deres_SolvsTA_control <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedg_SOLvsTA_control <- edgeR_DE(data_file, groups, adjp, drop_low, qlt_test)

#----------------------------------------------------------------------------------------------------------------------------------

# --- Soleus versus GA Control Tissue Comparison ---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/control/GAvsSOL'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.Sol_Flz_S1.txt','sorted.Sol_Flz_S2.txt','sorted.Sol_Flz_S3.txt',
                  'sorted.Sol_Flz_S4.txt', 'sorted.Sol_Flz_S5.txt','sorted.Sol_Flz_S6.txt',
                   'sorted.Ga_Flz_S1.txt','sorted.GA_Flz_S2.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'soleusandGAC-data-combined2.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2019/01_2019_Winter/PROJECT/RNAseq_analysis/2019_01_31/soleusandGAC-data-combined2.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2019/01_2019_Winter/PROJECT/RNAseq_analysis/2019_01_31/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('Sol 1','Sol 2','Sol 3','Sol 4','Sol 5','Sol 6', 'GA 1','GA 2')
# set up the two groups (in this case, I have 6 Sol replicates followed by 2 Ga replicates)
groups <- factor(c(1,1,1,1,1,1,2,2))
# condition vector - tells DESeq what contrast to do (I called 1 = untreated OR GROUP 1 and 2 = treated here OR GROUP 2; LOGFC(group1/group2))
condition <- c(rep("SOL",6),rep("GA",2))
adjp <- 0.05
drop_low <- TRUE

deres <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedg <- edgeR_DE(data_file, groups, adjp, drop_low)

#-----------------------------------------------------------------------------------------------------------------------------------

# --- TA versus GA Control Tissue Comparison ---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/control/GAvsTA'
# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.GF1_ga_flz13.txt','sorted.GF2_ga_m11.txt','sorted.TF1_Ta_Flz16.txt',
                  'sorted.TF2_Ta_Flz16.txt')
# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5
# name for the combined file (it is saved to data_path)
comb_file <- 'GAandTAC-data-combined.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_14/GAandTAC-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_14/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('Ga 1','Ga 2','TA 1','TA 2')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("GA",2),rep("TA",2))
adjp <- 0.05
drop_low <- TRUE

deres <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedg <- edgeR_DE(data_file, groups, adjp, drop_low)
