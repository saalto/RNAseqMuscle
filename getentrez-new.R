# Queries the mouse bioconductor DB with gene names and fetches corresponding
#   ENTREZ IDs.  If there is a 1-many mapping, a list of all IDS that hit are
#   returned.  If no match is found, NA is returned.

library('org.Mm.eg.db')

# --- make changes here ---
# directory (no trailing forward slash!) where the symbol file lives
data_dir <- '../data/2018.11.15/soleus-data'
# this should be a file with a list of gene symbols, one per line
name_file <- 'gene-symbols.txt'
# location to write the symbol-to-entrez key
save_dir <- '../data/2018.11.15/soleus-data'
# -------------------------

# vector of gene names
genes <- read.csv(paste(data_dir,name_file,sep='/'),header=FALSE)
genes <- as.vector(t$V1)

# symbol-to-entrez map
x <- mapIds(org.Mm.eg.db,t,'ENTREZID','SYMBOL',multiVals='list')

# save
saveRDS(et,file=paste(data_dir,'symbol-to-entrez.rds',sep='/'))
