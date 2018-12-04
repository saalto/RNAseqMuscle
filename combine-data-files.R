# Concatenates several count files into one large file with multiple columns; final file will also
#   lack the QC lines (last 5 or so)
library(tidyr)

# ---- change stuff here ----
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
# ---------------------------
# start with the first file
master <- read.csv(paste(data_path,files_to_agg[1],sep='/'),sep='\t',header=FALSE)
new_colnames <- c("V1","V2")

# bind each count column to the first file
for(i in 2:length(files_to_agg)){
  f <- read.csv(paste(data_path,files_to_agg[i],sep='/'),sep='\t',header=FALSE)
  master <- cbind(master,f$V2)
  new_colnames[i+1] <- paste("V",i+1,sep="")
}

# fix the column names
colnames(master) <- new_colnames

# remove the QC lines
n <- dim(master)[1]
master <- master[1:(n-n_qc),]

# save to a csv file
write.table(master,file=paste(data_path,comb_file,sep='/'),sep=",",
            quote=FALSE,row.names=FALSE,col.names=FALSE)
