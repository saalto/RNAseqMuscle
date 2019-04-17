# --- Location of Iteration #2 Processed Data ---
home_loc <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/'

# -----------------------------------------------------------------------------------------------------------------------

# --- Soleus versus TA Mutant Tissue Comparison ---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/mutant/TAvsSOL'

# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.Flz_Sol_S1.txt','sorted.Flz_Sol_S2.txt','sorted.Flz_Sol_S3.txt',
                  'sorted.Flz_Sol_S4.txt','sorted.Flz_Sol_S5.txt','sorted.Flz_Sol_S6.txt',
                  'sorted.Mck_Sol_S1.txt','sorted.Mck_Sol_S2.txt','sorted.Mck_Sol_S3.txt',
                  'sorted.Mck_Sol_S4.txt','sorted.Mck_Sol_S5.txt')

# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5

# name for the combined file (it is saved to data_path)
comb_file <- 'soleusandTAM-data-combined.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/soleusandTAM-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('Sol 1','Sol 2','Sol 3','Sol 4','Sol 5','TA 1','TA 2', 'TA 3', 'TA 4')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,1,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("Sol",5),rep("TA",4))


# --- Soleus versus GA Mutant Tissue Comparison ---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/mutant/GAvsSOL'

# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.Flz_Sol_S1.txt','sorted.Flz_Sol_S2.txt','sorted.Flz_Sol_S3.txt',
                  'sorted.Flz_Sol_S4.txt','sorted.Flz_Sol_S5.txt','sorted.Flz_Sol_S6.txt',
                  'sorted.Mck_Sol_S1.txt','sorted.Mck_Sol_S2.txt','sorted.Mck_Sol_S3.txt',
                  'sorted.Mck_Sol_S4.txt','sorted.Mck_Sol_S5.txt')

# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5

# name for the combined file (it is saved to data_path)
comb_file <- 'soleusandGAM-data-combined.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/soleusandGAM-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('Sol 1','Sol 2','Sol 3','Sol 4','Sol 5','GA 1','GA 2')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,1,1,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("Sol",5),rep("GA",2))
adjp <- 0.05
drop_low <- TRUE

# --- TA versus GA Mutant Tissue Comparison ---
# --- Use combine-data-files.R script to generate the combined text/Data variables ---
# location holding the count files - each count file is assumed to be one gene per line,
#   with a single count value (one integer), tab separated, after the gene symbol name
data_path <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/mutant/GAvsTA'

# list of files (filenames) to mush together - the final order depends on the order you
#   list them in this vector
files_to_agg <- c('sorted.Flz_Sol_S1.txt','sorted.Flz_Sol_S2.txt','sorted.Flz_Sol_S3.txt',
                  'sorted.Flz_Sol_S4.txt','sorted.Flz_Sol_S5.txt','sorted.Flz_Sol_S6.txt',
                  'sorted.Mck_Sol_S1.txt','sorted.Mck_Sol_S2.txt','sorted.Mck_Sol_S3.txt',
                  'sorted.Mck_Sol_S4.txt','sorted.Mck_Sol_S5.txt')

# number of quality control (non-gene) lines at the end of the file, in case this changes
n_qc <- 5

# name for the combined file (it is saved to data_path)
comb_file <- 'GAandTAM-data-combined.txt'

# --- Use DESeq2_DE.R and edgeR_DE.R scripts to generate Data variables to conduct counting analysis ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/GAandTAM-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('TA 1','TA 2','TA 3','TA 4','GA 1','GA 2')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,1,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("TA",4),rep("GA",2))
adjp <- 0.05
drop_low <- TRUE