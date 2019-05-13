# Need to load the jaccard.R and logseq.R functions before running this function.

# --- Change these conditions ---
data_file = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/soleus-data-combined.txt'
# directory to save results - no trailing forward slash!
save_path = 'C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_12/'
# column names; these are basically arbitrary and do not include the column of gene names
col_names <- c('wt 1','wt 2','wt 3','wt 4','wt 5','wt 6','mut 1','mut 2','mut 3','mut 4','mut 5')
# set up the two groups (in this case, I have 6 wt replicates followed by 5 mut replicates)
groups <- factor(c(1,1,1,1,1,1,2,2,2,2,2))
# condition vector - tells DESeq what contrast to do (I called wt = untreated and wt = treated here)
condition <- c(rep("untreated",6),rep("treated",5))
# threshold for significance
adjp <- 0.05
# drop genes with very low total counts? (recommended)
drop_low = TRUE

deres <- DESeq2_DE(data_file, col_names, condition, adjp, drop_low)
deedg <- edgeR_DE(data_file, groups, adjp, drop_low)
# --------------------------------------------------------------------------------------

s = logseq(0.1, 10^-5, 1000)

table1 <- data.frame("pvalue"=integer(),"condition"= integer(),"N"=integer(), "MedianLogFold" = integer())
table2 <- data.frame("pvalue"=integer(),"condition"= integer(),"N"=integer(), "MedianLogFold" = integer())
table4 <- data.frame("pvalue"=integer(), "JaccardIndex" = integer())

# --- Creating a data.frame of the DESeq and edgeR data ---
count= 0
for(i in s){
  # --- Creating the subset of data results ---
  DESE <- deres$table[which(deres$table$padj< i),]
  EDR <- deedg$table[which(deedg$table$PValue< i),]

  # --- Number of genes at specific p-values for each method---
  a <- nrow(DESE)
  d <-c("DESeq")
  b <- nrow(EDR)
  e <-c("edgeR")
  
  # --- Calculating the median of the log-fold change for each method ---
  deseq <- median(abs(DESE$log2FoldChange))
  edge <- median(abs(EDR$logFC))
  
  # --- Isolate the gene symbols for each method---
  xrow <- row.names(DESE)
  yrow <- row.names(EDR)
  
  # --- Filling in the tables ---
  table1[nrow(table1)+1,] <- c(i,d,a,deseq)
  table2[nrow(table2)+1,] <- c(i,e,b,edge)
  
  # --- Jaccard Index value at specific p-values ---
  f <- jaccard(xrow,yrow)
  table4[nrow(table4)+1,] <- c(i,f)
  
  count = count + 1
}

# --- Removing the NA at the end of the tables ---
table1 <- table1[1:count-1,]
table2 <- table2[1:count-1,]
table4 <- table4[1:count-1,]

# --- Combine DESeq and edgeR tables ---
table3 <- rbind(table1,table2)

# --- Converting characters into numericals ---
table3$pvalue <-as.numeric(table3$pvalue)
table3$MedianLogFold <- as.numeric(table3$MedianLogFold)
table3$N <- as.numeric(table3$N)

table4$pvalue <-as.numeric(table4$pvalue)
table4$JaccardIndex <- as.numeric(table4$JaccardIndex)

# --- Data Visualization ---
library("ggplot2")
library("gridExtra")
library(ggpubr)
library(RColorBrewer)
p1 <- ggplot(table3, aes(x= -log10(pvalue), y=N, color=condition)) + geom_point() +scale_x_log10()+theme(legend.position="none", axis.title.x = element_blank())+scale_color_manual(values = c("red", "blue")) + geom_vline(xintercept = -log10(0.01))
p2 <- ggplot(table3, aes(x= -log10(pvalue), y=MedianLogFold, color=condition)) + geom_point() + scale_x_log10()+theme(legend.position="none") + scale_color_manual(values = c("red", "blue")) +geom_vline(xintercept = -log10(0.01))
p3 <- ggplot(table4, aes(x= -log10(pvalue), y=JaccardIndex)) + geom_point() + scale_x_log10()+theme(axis.title.x = element_blank()) + geom_vline(xintercept = -log10(0.01))

ggarrange(p1, p3, p2, ncol=1, nrow = 3, common.legend = TRUE, legend = "bottom")

