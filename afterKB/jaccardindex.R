setwd("C:/Users/sarah/OneDrive/Documents/2018/04_2018_Fall/RNAseq_analysis/2018_12_05")
x <- readRDS('de-deseq-results.rds')
y <- readRDS('de-edgeR-results.rds')

table4 <- data.frame("pvalue"=integer(),"condition"= integer(), "jaccard" = integer())

# --- Creating a data.frame of the DESeq data ---
count= 0

for(i in seq(0.1, 0.0001, -0.0001)){
  # --- Creating the subset of data results ---
  x <- x[which(x$padj< i),]
  y <- y[which(y$PValue< i),]
  
  # --- Isolate the gene symbols ---
  xrow <- row.names(x)
  yrow <- row.names(y)
  
  # --- Jaccard Index value at specific p-values ---
  a <- jaccard(xrow,yrow)
  d <-c("DESeq")

  count = count + 1
  table4[nrow(table4)+1,] <- c(i,d,a)
  
}

table4 <- table4[1:count-1,]
#table4

# --- Combine DESeq and edgeR tables ---
table4$pvalue <-as.numeric(table4$pvalue)
table4$jaccard <- as.numeric(table4$jaccard)


library("ggplot2")
library(scales)
p4 <- ggplot(table4, aes(x=(-log(pvalue)), y=jaccard)) + geom_point() + scale_x_log10()+ geom_vline(xintercept = -log(0.05))
