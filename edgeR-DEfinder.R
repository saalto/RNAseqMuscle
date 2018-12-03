# Computing differentially expressed genes with edgeR
library("edgeR")

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
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,1,1,1,2,2,2,2,2))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE
# -------------------------
# read in count data
# NOTES:
#   +The last 6 or so quality control lines (__ambiguous, __alignmen_not_unique, etc) have been
#     removed from this file
#   +This is both wt and mut conditions, column concatenated (wt first, then mutant)
#   +DO NOT open this in Excel, as it interprets some gene names (March1 - March10) as dates

cts <- read.csv(data_file,header=FALSE)

# separate into gene names and just counts
genes <- cts[,c(1)]
cts <- cts[,c(2:12)]

# set the rownames as the gene symbols
rownames(cts) <- genes

# pre-filtering of stuff with very few counts
if(drop_low){
  keep <- rowSums(cts) >= 10
  cts <- cts[keep,]
}
# save the gene universe
universe <- rownames(cts)

# now do the calculations
y <- DGEList(counts=cts,group=groups)
y <- calcNormFactors(y)
design <- model.matrix(~groups)
y <- estimateDisp(y,design)

# here are the tests for significance
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)

# now extract the genes we want
tagtable <- topTags(qlf,n=length(qlf$df.total))

# save out DE genes, universe, and results table
deindex <- which(tagtable$table$FDR < adjp, arr.ind=TRUE)
genelist <- row.names(tagtable$table)[deindex]
saveRDS(genelist,file=paste(save_path,'de-edgeR-genes.rds',sep='/'))
saveRDS(universe,file=paste(save_path,'de-deseq-universe.rds',sep='/'))
saveRDS(tagtable$table,file=paste(save_path,'de-edgeR-results.rds',sep='/'))
