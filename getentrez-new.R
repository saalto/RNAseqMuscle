# Queries the mouse bioconductor DB with gene names and fetches corresponding
#   ENTREZ IDs.  If there is a 1-many mapping, a list of all IDS that hit are
#   returned.  If no match is found, NA is returned.

library('org.Mm.eg.db')

# --- make changes here ---
# directory (no trailing forward slash!) where the symbol file lives
data_dir <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus'

# this should be a file with a list of gene symbols, one per line
#name_file <- 'gene-symbols.txt'

# location to write the symbol-to-entrez key
save_dir <- 'C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus'
# -------------------------

# --- for-loop to go through each text file and generate a specific rds file ---
files <- dir(data_dir, pattern ='*.txt')

for (i in 1:length(files)){
  # vector of gene names
  textfile <- read.delim(paste(data_dir,files[i],sep='/'), header=FALSE)
  textfile <- as.vector(textfile$V1)

  # selecting the genes with <NA>, which stands for multiple IDs
  t <- select(org.Mm.eg.db,textfile, 'ENTREZID', 'SYMBOL')
  t$V3 <- is.na(t$ENTREZID)
  #length(t$V3[t$V3 == TRUE])
  multiples <-subset(t, V3 == TRUE)
  multiples <- multiples$SYMBOL
  x <- mapIds(org.Mm.eg.db,multiples,'ENTREZID','SYMBOL',multiVals='list')
  # keys entered are not valid keys for 'SYMBOL'
  
  # save, 
  outfile <- paste(i,"symbol-to-Multipleentrez.rds",sep = '_')
  saveRDS(x,file=paste(data_dir,outfile,sep='/'))
  
  # symbol-to-Oneentrez map
  x <- mapIds(org.Mm.eg.db,textfile,'ENTREZID','SYMBOL',multiVals='first')
    # save, 
  outfile <- paste(i,"symbol-to-Oneentrez.rds",sep = '_')
  saveRDS(x,file=paste(data_dir,outfile,sep='/'))

}

b <-readRDS('1_symbol-to-entrez.rds')
View(b)

