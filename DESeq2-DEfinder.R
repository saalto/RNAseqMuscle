# Computing significantly overexpressed genes using DESeq2
library("DESeq2")

# ----- make changes here ----
# a very specific format is assumed for the RNASeq count data
#   +comma separated value file
#   +no headers
#   +first column is gene symbols
#   +all conditions (say, wt and mut for a single tissue) are concatenated
#   +DO NOT open this file in Excel before doing the R read, since then
#     Excel interprets some gene names (March1 - March10) as dates

data_file = '../../data/2018.11.15/soleus-data/soleus-all-data.txt'
# directory to save results - no trailing forward slash!
save_path = '..'
# column names; these are basically arbitrary and do not include the column of gene names
colnames <- c('wt 1','wt 2','wt 3','wt 4','wt 5','wt 6','mut 1','mut 2','mut 3','mut 4','mut 5')
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("untreated",6),rep("treated",5))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE

# -------------------------

# read in count data
cts <- read.csv(data_file,header=FALSE)

# separate into gene names and just counts
genes <- cts[,c(1)]
cts <- cts[,c(2:12)]

# change row/col names (script will fail here if Excel was used to open and save the data file)
rownames(cts) <- genes
colnames(cts) <- colnames

coldata <- data.frame(colnames(cts), "condition" = condition)

dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design = ~ condition)

# pre-filtering of stuff with very few counts
if(drop_low) {
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
}

# run the analysis, with the requested p-value cutoff
dds <- DESeq(dds)
res <- results(dds,alpha=adjp)

# shinkage (useful for visualization and ranking)
#resLFC <- lfcShrink(dds,coef="condition_untreated_vs_treated",type="apeglm")
# order results table by smallest p-value
#resOrdered <- res[order(res$padj)]

# save out universe and gene lists for OA analysis
deindex <- which(res$padj < adjp, arr.ind=TRUE)
genelist <- row.names(res)[deindex]
universe <- row.names(res)
saveRDS(genelist,file=paste(save_path,'de-deseq-genes.rds',sep='/'))
saveRDS(universe,file=paste(save_path,'de-deseq-universe.rds',sep='/'))
# save full results table
saveRDS(res,file=paste(save_path,'de-deseq-results.rds',sep='/'))
