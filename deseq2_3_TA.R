#---Clear data and load packages-----------------------
rm(list = ls())
dev.off()
source('https://bioconductor.org/biocLite.R')

# Note: MASS, stringi, and survival packages had to be 'updated' after loading each library. Possible bug with new version of R.
biocLite('DESeq2')
install.packages("stringi",type="win.binary")
library("DESeq2")
library("ggplot2")
biocLite("ComplexHeatmap")
library("ComplexHeatmap")
biocLite("pamr")
library("pamr")
biocLite("MCL")
library("MCL")

#---Set working directory-----------------------------
setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/RNAseq_TAtextfiles/")
directory <- "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/RNAseq_TAtextfiles/"

#---Set up DESeq2 data, based on names of HTSeq counts in working directory-----
sampleFiles <- dir(pattern = 'sorted')
print(sampleFiles)
#sample group set up
ConditionMatch <- regexpr(pattern = '[A-Z]+', dir(pattern = '.txt'))
print(ConditionMatch)
sampleConditions <- regmatches(dir(pattern = '*.txt'), ConditionMatch)
print(sampleConditions)
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleConditions)
print(sampleTable)

#---Calculate DESeq2 from HTSeq count tables-------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)
#print(ddsHTSeq)

#-------Filter out genes with zero counts--------------------------
ddsHTSeqFiltered <- ddsHTSeq [ rowSums(counts(ddsHTSeq)) > 0, ]
ddsHTSeqFiltered <- DESeq(ddsHTSeqFiltered)
#print(ddsHTSeqFiltered)

#---Visualizations of Data Transformations--------------------------------
#Generate log2 normalized count matrix and variance stabilizing matrix
rld <- rlog(ddsHTSeqFiltered, blind = FALSE)
vsd <- varianceStabilizingTransformation(ddsHTSeqFiltered, blind = FALSE)
class(vsd)
head(colData(vsd),3)
logTransCounts <- assay(rld)
head(logTransCounts)

#---Comparing heatmaps of the three visualization methods------------------
library("RColorBrewer")
install.packages("gplots")
library("gplots")

select <- order(rowMeans(counts(ddsHTSeqFiltered,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap.2(counts(ddsHTSeqFiltered,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

#---Comparing samples to each other for correlation patterning-----------
heatmap.2(mat, trace = "none", col = rev(hmcol), margins = c(13,13))
dev.off()

#---The SHORT way to generate PCA plots---------
plotPCA(vsd, "condition")
plotPCA(rld, "condition")

#---MDS plot using Euclidean distances----
distsRL <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(distsRL)
mdsData <- data.frame(cmdscale(DistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3)

#---Poisson Distance Plot------------
install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsHTSeqFiltered)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(ddsHTSeqFiltered)))
ggplot(mdsPois, aes(X1,X2,color=condition)) + geom_point(size=3)

#---Singling out specific genes based on log counts---------------------
logTransCounts[grep("Pitx2", rownames(logTransCounts)), ]
logTransCounts[grep("Foxo3", rownames(logTransCounts)), ]
logTransCounts[grep("Arf1", rownames(logTransCounts)), ]
#rld <- rlog (ddsHTSeqFiltered, blind = F, logTransCounts[grep("Pitx2", rownames(logTransCounts)), ])
#logTransCounts[grep("Foxo3", logTransCounts), ]
#mcols(object)$dispFit

#---The LONG way to calculate Principle component analysis based on log2 normalized count matrix---
OGPCAN <-prcomp(logTransCounts, center = T, scale = F, tol = 0)
print(OGPCAN)
OGPCAN_matrix <- as.data.frame(OGPCAN$rotation)
OGPCAN_matrix <- OGPCAN_matrix[, 1:7]
print(OGPCAN_matrix)
#PCA plot replicates set up
OGPCAN_matrix$Condition <- c(rep("TF", 3), rep("TM", 4)) 
#rep("m1", 1), rep("m2", 1), rep("m3", 1))
#print(OGPCAN_matrix)
#Plot PCA
ggplot(OGPCAN_matrix, aes(PC2, PC1, color = Condition)) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size  = 16, face = "bold"),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14)) +
  #scale_y_continuous(limits = c(0.3, 0.4)) +
  #scale_x_continuous(limits = c(-5, 5)) +
  scale_color_discrete(name = "Sample")
#print(ddsHTSeqFiltered)

#---Calculate differentially expressed genes from DESeq2---------------------
res.TA_W_M <- results(ddsHTSeqFiltered, contrast = c("condition", "TM", "TF"))
mcols(res.TA_W_M, use.names=TRUE)
summary(res.TA_W_M)

#Calculate differentially expressed genes from DESeq2 object based on adjusted p-value
res.TA_W_M.05 <- results(ddsHTSeqFiltered, alpha=0.05)
table(res.TA_W_M.05$padj < 0.05)

#Calculate differentially expressed genes from DESeq2 object based on log fold change equal to 0.5 (2^0.5)
res.TA_W_MLFC1 <- results(ddsHTSeqFiltered, lfcThreshold=0.5)
table(res.TA_W_MLFC1$padj < 0.1)

#---MA plot from res---------------------------------------
plotMA.TA_W_M <- plotMA(res.TA_W_M, ylim=c(-2,2))
mcols(res.TA_W_M, use.names = TRUE)

#---filter DE genes using subsets----------------------------------
#Set FDR threshold of p <= 0.05
res.TA_W_M_filtered <- as.data.frame(as.matrix(subset(res.TA_W_M, pvalue <= 0.05)))
#head(res.TA_W_M_filtered)
#dim(res.TA_W_M_filtered)

#Set adjusted p-value to 0.05
res.TA_W_M_filtered <- as.data.frame(as.matrix(subset(res.TA_W_M, padj <= 0.05)))
#head(res.TA_W_M_filtered)
#dim(res.TA_W_M_filtered)

#Set log fold-change and adjusted p-value
res.TA_W_M_filtered$absFC <- abs(res.TA_W_M_filtered$log2FoldChange)
res.TA_W_M_filtered2 <- subset(res.TA_W_M_filtered, absFC > 0.5)
res.TA_W_M_filtered3 <- subset(res.TA_W_M_filtered2, padj < 0.05)
#head(res.TA_W_M_filtered)
#dim(res.TA_W_M_filtered)
#head(res.TA_W_M_filtered2)
#dim(res.TA_W_M_filtered2)
#head(res.TA_W_M_filtered3)
#dim(res.TA_W_M_filtered3)

#---Print out the filtered data as a text file----------------------------
write.table(res.TA_W_M_filtered, "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/RNAseq_analysis/res.TA_W_M_filtered_20180628.txt", sep ="\t")
#write.csv(as.data.frame(res.TA_W_M_filtered), file = "")

#---Heatmap of the 30 most significant fold-change genes--------------
install.packages("pheatmap")
library("pheatmap")
Mat <- assay(vsd)[ head(order(res.TA_W_M$padj),30), ]
Mat <- Mat - rowMeans(Mat)
df <- as.data.frame(colData(vsd)[,c("condition")])
pheatmap(Mat)

#---Heatmap plot of 200 gene expressions and significant Log2fold changes-----
var_genes <- apply(logTransCounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]
highly_variable_lcpm <- logTransCounts[select_var,]
Heatmap(highly_variable_lcpm)
#nrow(res.TA_W_M_filtered)
#ncol(res.TA_W_M_filtered) 
#head(res.TA_W_M_filtered)
#head(var_genes)
#head(select_var)
#dim(highly_variable_lcpm)
#head(highly_variable_lcpm)

#---Heatmap of most significant genes from filtered-out dataset of the overall fold-change between the mutant and wildtype groups---
TA_W_M.heatmap.matrix <- as.matrix(res.TA_W_M_filtered[ ,c(2)])
Heatmap(TA_W_M.heatmap.matrix)
Heatmap(TA_W_M.heatmap.matrix, cluster_columns = F) #Formats based on alphabetical order
#head(res.TA_W_M_filtered)
#class(res.TA_W_M_filtered)
#class(TA_W_M.heatmap.matrix)
#print(TA_W_M.heatmap.matrix)
#TA_W_M.heatmap.matrix <- t(TA_W_M.heatmap.matrix) # flip rows and columns around

#---Heatmap of the most significant genes for each individual from filtered-out dataset-----
var_genes <- apply(res.TA_W_M_filtered, 1, var)
head(var_genes)
select_var <- names(sort(var_genes))
head(select_var)
highly_variable_lcpm <- logTransCounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
Heatmap(highly_variable_lcpm)

#---Clear data and load packages-----------------------
rm(list = ls())
dev.off()

#======================================================
#============== Arun's Code ===========================

e12.13.res <- as.data.frame(as.matrix(subset(E12.5_E13.5_final, padj < 0.05)))
e12.13.res$absFC <- abs(e12.13.res$log2FoldChange)
e12.13.res <- subset(e12.13.res, absFC >= log(1.5, 2))
e12.13.res <- e12.13.res[order(e12.13.res$padj), ]

e13.14.res <- as.data.frame(as.matrix(subset(E13.5_E14.5_final, padj < 0.05)))
e13.14.res$absFC <- abs(e13.14.res$log2FoldChange)
e13.14.res <- subset(e13.14.res, absFC >= log(1.5, 2))
e13.14.res <- e13.14.res[order(e13.14.res$padj), ]

#List all genes DE between any two sequential stages
all.de.genes <- unique(c(rownames(e11.12.res), rownames(e12.13.res), rownames(e13.14.res)))

#Order samples by approximate order from PCA plot
sample.order.1 <- c(1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,11,20,21,22,23,24,25)
sample.order.2 <- c(5,3,4,6,2,1,7,9,8,10,14,12,15,17,18,16,13,11,19,20,21,25,24,22,23)

logTransCountsOrdered <- logTransCounts[, sample.order.1]
logTransCountsOrdered <- logTransCounts[, sample.order.2]

#Correlation analysis

#Function to remove genes with 0 variance between samples

var_filter <- function(df) {
  select_vec <- logical()
  for(i in 1:length(rownames(df))) {
    if(sd(df[i,]) == 0) {
      select_vec <- append(select_vec, F)
    } else {
      select_vec <- append(select_vec, T)
    }
  }
  return(df[select_vec,])
}

#Create normalized count table of only DE genes of interest
sig.counts <- logTransCountsOrdered[match(all.de.genes, rownames(logTransCountsOrdered)), ]

#Remove genes with 0 variance
sig.counts.filtered <- var_filter(sig.counts)

#Function to calculate pairwise correlation between each gene in a data frame
#Input is a data frame or matrix with genes as rows and samples as columns
#Output is a data frame with four columns
#Gene1, P-value adjusted with fdr, correlation (pearson correlation coefficient), Gene2
correlation_test <- function(df1) {
  cor_vec <-numeric() #Set vectors to be used to reduce memory use
  p_vec <- numeric()
  whole_vec <- numeric()
  name_vec <- character()
  final_df <- data.frame()
  count_r1 <- 1
  while(count_r1 <= length(rownames(df1)) - 1) {
    for(i in 1:(length(rownames(df1)) - count_r1)) {
      whole_vec <- as.numeric(cor.test(as.matrix(df1[count_r1,]), as.matrix(df1[i+count_r1,]), 't', 'pearson', exact=NULL)[1:5])
      cor_vec <- append(cor_vec, whole_vec[4])
      p_vec <- append(p_vec, whole_vec[3])
      name_vec <- append(name_vec, rownames(df1)[i + count_r1])
      new_pvec <- p.adjust(p_vec, 'fdr')
      og_name_vec <- rownames(df1)[count_r1]
      new_df <- data.frame(name_vec, new_pvec, cor_vec)
      new_df$og_name_vec <- og_name_vec
      final_df <- rbind(final_df, new_df)
      p_vec <- numeric()
      count_r1 <- count_r1 + 1
      name_vec <- character()
      cor_vec <- numeric()
    }
    colnames(final_df) <- c("Gene2", "pAdj", "Correlation", "Gene 1")
    final_df <- final_df[final_df$pAdj <= 0.1,]
    return(final_df)
  }}

#calculate correlation dataframe
corr_df <- correlation_test(sig.counts.filtered)

#Create a function that fits a power law to the correlation dataframe calculated earlier
#Power fit
powerFit <- function(df, p.vec) {
  #Pre-define vectors for more efficient memory
  r.sq.vec <- numeric()
  p.thresh.vec <- numeric()
  edge.vec <- integer()
  node.vec <- integer()
  ave.deg.vec <- numeric()
  for(i in 1:length(p.vec)) {
    #Iterate through the p-value vector and calculate network statistics
    new.df <- subset(df, pAdj <= p.vec[i])
    tmp.genes <- c(as.character(new.df$Gene2), as.character(new.df$`Gene 1`))
    n.nodes <- length(unique(tmp.genes))
    n.edges <- length(new.df$Gene2)
    ave.deg <- 2 * n.edges / n.nodes
    nodes <- table(as.integer(table(tmp.genes)))
    node.df <- data.frame(as.integer(names(nodes)), as.integer(nodes))
    colnames(node.df) <- c("degree", "nodes")
    tmp.mod <- lm(formula = log(degree, 10) ~ log(nodes, 10), data = node.df)
    r.squared <- as.numeric(summary(tmp.mod)[8])
    r.sq.vec <- append(r.sq.vec, r.squared)
    p.thresh.vec <- append(p.thresh.vec, p.vec[i])
    edge.vec <- append(edge.vec, n.edges)
    node.vec <- append(node.vec, n.nodes)
    ave.deg.vec <- append(ave.deg.vec, ave.deg)
    new.df <- data.frame()
  }
  df.final <- data.frame(p.thresh.vec, node.vec, edge.vec, ave.deg.vec, r.sq.vec)
  colnames(df.final) <- c("Pval", "Nodes", "Edges", "AveDegree", "R.squared")
  return(df.final)
}

#Create a function that generates a vector of p values to be used to determine an appropriate cutoff for generating a coexpression network
pseq <- function(start, long) {
  count = 1
  p.vec <- numeric()
  new.num <- start
  while(count <= long) {
    new.num <- new.num/2
    p.vec <- append(p.vec, new.num)
    new.num <- new.num/5
    p.vec <- append(p.vec, new.num)
    count = count + 1
  }
  return(p.vec)
}

pvec <-pseq(0.1, 30) 
#Create a dataframe with R2 values vs. p-value cutoff threshold
model.fit <- powerFit(cor_df, pvec)


#Create correlatin dataframe based off powerfit to be exported to cytoscape
finalCorrDf <- subset(cor_df, pAdj <= 1e-16)
