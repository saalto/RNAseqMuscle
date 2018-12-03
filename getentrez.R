# Queries the mouse bioconductor DB with gene names and fetches corresponding
#   ENTREZ IDs.  Currently, does not deal gracefully with situations in which
#   the mapping is one-to-many.
# Proper way to do that is:
#   1. Fetch each ID one-at-a-time
#   2. If only one match is returned, add that symbol:ID pair to the dataframe
#   3. If more than one is returned, do reverse (ID->symbol) matches to figure
#       out which ID to use.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler)

# --- make changes here ---
# directory (no trailing forward slash!) where the symbol file lives
data_dir <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus'

# location to write the symbol-to-entrez key
save_dir <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus'

# # this should be a file with a list of gene symbols, one per line
# name_file <- load(textfile.name)

# in case you need to analyze a different organism
db_name <- "org.Mm.eg.db"
# -------------------------

# --- for-loop to go through each text file and generate a specific rds file ---
files <- dir(data_dir, pattern ='*.txt')

for (i in 1:length(files)){
  # vector of gene names
  textfile <- read.delim(paste(data_dir,files[i],sep='/'), header=FALSE)
  # t <- read.delim(paste(data_dir,name_file,sep='/'),header=FALSE)
  textfile <- as.vector(textfile$V1)
  
  # symbol-to-entrez map
  et <- bitr(textfile,fromType="SYMBOL",toType="ENTREZID",OrgDb=db_name)
  
  # save
  outfile <- paste(i,"symbol-to-entrez.rds",sep = '_')
  saveRDS(et,file=paste(data_dir,outfile,sep='/'))
}
