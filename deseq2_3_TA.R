#
##---Clear data and load packages-----------------------
rm(list = ls())
dev.off()
source('https://bioconductor.org/biocLite.R')
library("DESeq2")
library("ggplot2")
library("ComplexHeatmap")
library("pamr")
library("MCL")

##---Set working directory to iteration 1---
#setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration1/RNAseq_TAtextfiles/")
#directory <- "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration1/RNAseq_TAtextfiles/"

##---Set working directory to iteration 2---
setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/")
directory <- "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/"

##---Set up DESeq2 data, based on names of HTSeq counts in working directory---
sampleFiles <- dir(pattern = 'sorted')
#print(sampleFiles)

#---sample group set up---
ConditionMatch <- regexpr(pattern = '[A-Z]+', dir(pattern = '.txt'))
#print(ConditionMatch)
sampleConditions <- regmatches(dir(pattern = '*.txt'), ConditionMatch)
#print(sampleConditions)
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleConditions)
#print(sampleTable)


##-----ONLY USE IF USING VERA'S TEXT FILES LOCATED IN ITERATION 1-----
# #-reassignment of the files based on the Foxo3 gene expression patterning,
# #-which is higher in the flox/Z mice than the MCK mice. this gene expression 
# #-is a better indicator of genotypes than the labeling provided.
# #-DESeq2 object and rld arrary was generated first before this line to observe
# #-the Foxo3 gene expression.
# sampleTable$condition <- c(rep("TF", 1), rep("TM",1), rep("TF", 1), rep("TM",4))
# #print(sampleTable)
  
##---Calculate DESeq2 from HTSeq count tables-------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)
#print(ddsHTSeq)

#-Filter out genes with zero counts--------------------------
ddsHTSeqFiltered <- ddsHTSeq [ rowSums(counts(ddsHTSeq)) > 0, ]
ddsHTSeqFiltered <- DESeq(ddsHTSeqFiltered)
#print(ddsHTSeqFiltered)

##---Visualizations of Data Transformations--------------------------------
#-Generate log2 normalized count matrix and variance stabilizing matrix
rld <- rlog(ddsHTSeqFiltered, blind = FALSE)
vsd <- varianceStabilizingTransformation(ddsHTSeqFiltered, blind = FALSE)
#class(vsd)
#head(colData(vsd),3)
logTransCounts <- assay(rld)
#head(logTransCounts)

#---Singling out specific genes based on log counts---------------------
logTransCounts[grep("Pitx2", rownames(logTransCounts)), ] #Wildtype and mutant samples have similar counts
logTransCounts[grep("Foxo3", rownames(logTransCounts)), ] #Wildtype samples will have higher counts than mutants
#logTransCounts[grep("Arf1", rownames(logTransCounts)), ]


#---Comparing heatmaps of the three normalizing methods------------------
library("RColorBrewer")
install.packages("gplots")
library("gplots")

select <- order(rowMeans(counts(ddsHTSeqFiltered,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

#pdf("readcountsTranform.pdf")
heatmap.2(counts(ddsHTSeqFiltered,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6), main = "Read Counts Transformation")
dev.off()

#pdf("rlogTransform.pdf")
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6), main = "rLog Transformation")
dev.off()

#pdf("VariStablizeTransform.pdf")
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6), main = "Variance Stablizing Transformation")
dev.off()


#---Comparing samples to each other for correlation patterning-----------
#-using rlog transformation (rld)---
#TO-DO: print off the figures
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <-  with(colData(ddsHTSeqFiltered), paste(condition, type, sep = " : "))
heatmap.2(mat, trace = "none", col = rev(hmcol), margins = c(13,13))
dev.off()

#-using variance stabilizing transformation (VSD)---
#TO-DO: print off the figures
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <-  with(colData(ddsHTSeqFiltered), paste(condition, type, sep = " : "))
heatmap.2(mat, trace = "none", col = rev(hmcol), margins = c(13,13), main = "Correlation between Samples based on Variance Stabilizing Transformation")
dev.off()

##---MDS plot using Euclidean distances 
#-rLog Transformation (rld)----
#TO-DO: print off the figures
distsRL <- dist(t(assay(rld)))
DistMatrix <- as.matrix(distsRL)
mdsData <- data.frame(cmdscale(DistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3) + ggtitle("MDS using Euclidean distance and rLog Transformation")
dev.off()

#-Variance Stablizing Transformation (vsd)----
#TO-DO: print off the figures
distsRL <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(distsRL)
mdsData <- data.frame(cmdscale(DistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3) + ggtitle("MDS using Euclidean distances and Variance Stabilizing Transformation")
dev.off()

##---Poisson Distance Plot------------
#TO-DO: print off the figures
install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsHTSeqFiltered)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(ddsHTSeqFiltered)))
ggplot(mdsPois, aes(X1,X2,color=condition)) + geom_point(size=3)


##---Principle component analysis based on log2 normalized count matrix---
#TO-DO: print off the figures
OGPCAN <-prcomp(logTransCounts, center = T, scale = F, tol = 0)
#print(OGPCAN)
OGPCAN_matrix <- as.data.frame(OGPCAN$rotation)
#OGPCAN_matrix <- OGPCAN_matrix[, 1:7]
#print(OGPCAN_matrix)

#---PCA plot replicates set up---
OGPCAN_matrix$Condition <- sampleTable$condition
#rep("m1", 1), rep("m2", 1), rep("m3", 1)) 
#print(OGPCAN_matrix)

#---Plot PCA---
ggplot(OGPCAN_matrix, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size  = 16, face = "bold"),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14)) +
  #scale_y_continuous(limits = c(0.3, 0.4)) +
  #scale_x_continuous(limits = c(-5, 5)) +
  scale_color_discrete(name = "Sample") +
  ggtitle("Principle Component Analysis based on rlog transformation")
#print(ddsHTSeqFiltered)
dev.off()


#---Calculate differentially expressed genes from DESeq2---------------------
res.TA_W_M <- results(ddsHTSeqFiltered, contrast = c("condition", "TM", "TF"))
mcols(res.TA_W_M, use.names=TRUE)
#summary(res.TA_W_M)

#---MA plot from results---------------------------------------
#-Bland-Altman plot that visualizes the differences between measurements taken in two samples, 
#-by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values
plotMA.TA_W_M <- plotMA(res.TA_W_M, ylim=c(-2,2))
mcols(res.TA_W_M, use.names = TRUE)

##---Volcano Plot of the results----
#TO-DO: print off the figures
#-Measuring the effect of fold change and the statistical significance
with(res.TA_W_M, plot(log2FoldChange, -log10(pvalue), pch=10, main="Volcano plot", xlim=c(-5,7)))

#-Add colored points: red if padj<0.05, orange if log2FC>1, green if both
with(subset(res.TA_W_M, padj<0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res.TA_W_M, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res.TA_W_M, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

#-Label points with the textxy function from the calibrate plot
#TO-DO: find a way to separate out the labels
#install.packages("calibrate")
#library(calibrate)
#with(subset(res.TA_W_M, -log10(pvalue)>120 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=rownames(res.TA_W_M), cex=1))
with(subset(res.TA_W_M), identify(log2FoldChange, -log10(pvalue), labels=rownames(res.TA_W_M))) #Need to click on graphic to label the outliers
dev.off()


##---Filtering the Results (DE genes) into tables---
#-Calculate differentially expressed genes from DESeq2 object based on adjusted p-value
res.TA_W_M.05 <- results(ddsHTSeqFiltered, alpha=0.05)
table(res.TA_W_M.05$padj < 0.05)

#-Calculate differentially expressed genes from DESeq2 object based on log fold change equal to 0.5 (2^0.5)
res.TA_W_MLFC1 <- results(ddsHTSeqFiltered, lfcThreshold=0.5)
table(res.TA_W_MLFC1$padj < 0.1)

#-Set p-value less than 0.05
#res.TA_W_M_filtered <- as.data.frame(as.matrix(subset(res.TA_W_M, pvalue <= 0.05)))
#head(res.TA_W_M_filtered)
#dim(res.TA_W_M_filtered)

#-Set adjusted p-value to less than 0.05
res.TA_W_M_filtered <- as.data.frame(as.matrix(subset(res.TA_W_M, padj <= 0.05)))
#head(res.TA_W_M_filtered)
#dim(res.TA_W_M_filtered)

#-Subset data based on (1) adjusted p-value less than 0.05 (2) absolute value of the log2 fold change greater than 0.5
res.TA_W_M_filtered$absFC <- abs(res.TA_W_M_filtered$log2FoldChange)
#head(res.TA_W_M_filtered)
#nrow(res.TA_W_M_filtered)
#ncol(res.TA_W_M_filtered)
#dim(res.TA_W_M_filtered)

res.TA_W_M_filtered2 <- subset(res.TA_W_M_filtered, absFC > 1)
#head(res.TA_W_M_filtered2)
#dim(res.TA_W_M_filtered2)

#res.TA_W_M_filtered3 <- subset(res.TA_W_M_filtered2, padj < 0.05)
#head(res.TA_W_M_filtered3)
#dim(res.TA_W_M_filtered3)


#---Print out the filtered data as a text file---
write.table(res.TA_W_M_filtered2, "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis/res.TA_W_M_filtered_padjfoldchange_20180820.txt", sep ="\t")
write.table(res.TA_W_M_filtered, "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis/res.TA_W_M_filtered_padj_20180820.txt", sep ="\t")
write.table(res.TA_W_M, "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis/res.TA_W_M_20180820.txt", sep = "\t")


##---Heatmap of the most significant fold-change genes--------------
install.packages("pheatmap")
library("pheatmap")

# #Heatmap of the significant (padj< 0.05) and fold-change (>1) DE genes based on RLD
# Map <- assay(rld)[order(res.TA_W_M_filtered2$padj), ]
# Map <- Map - rowMeans(Map)
# df <- as.data.frame(colData(rld)[,c("condition")])
# pheatmap(Map)
# dev.off()
# 
# #Heatmap of the significant (padj<0.05) DE genes based on RLD
# Map <- assay(rld)[order(res.TA_W_M_filtered$padj), ]
# Map <- Map - rowMeans(Map)
# df <- as.data.frame(colData(rld)[,c("condition")])
# pheatmap(Map)
# dev.off()
# 
# #Heatmap of the significant (padj< 0.05) and fold-change (>1) DE genes based on VSD
# Mat <- assay(vsd)[order(res.TA_W_M_filtered2$padj), ]
# Mat <- Mat - rowMeans(Mat)
# df <- as.data.frame(colData(vsd)[,c("condition")])
# pheatmap(Mat)
# dev.off()

#Heatmap of the significant (padj<0.05) DE genes based on VSD
library(RColorBrewer)
Mat <- assay(vsd)[order(res.TAL_W_M_filtered$padj), ]
Mat <- Mat - rowMeans(Mat)
df <- as.data.frame(colData(vsd)[,c("condition")])
pheatmap(Mat, color= colorRampPalette(c("#990066", "#ffffff", "#009933"))(14),show_rownames = F, show_colnames = F)
dev.off()


# #---Heatmap plot of 200 gene expressions and significant Log2fold changes-----
# # TO DO: Work in progress...I am not seeing what I should be seeing
# var_genes <- apply(logTransCounts, 1, var)
# #head(var_genes)
# select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]
# #head(select_var)
# highly_variable_lcpm <- logTransCounts[select_var,]
# #dim(highly_variable_lcpm)
# #head(highly_variable_lcpm)
# Heatmap(highly_variable_lcpm)
# 
# #---Heatmap of the most significant genes for each individual from filtered-out dataset-----
# var_genes <- apply(res.TA_W_M_filtered, 1, var)
# #head(var_genes)
# select_var <- names(sort(var_genes))
# #head(select_var)
# highly_variable_lcpm <- logTransCounts[select_var,]
# #dim(highly_variable_lcpm)
# #head(highly_variable_lcpm)
# Heatmap(highly_variable_lcpm)

#---Clear data and load packages-----------------------
rm(list = ls())
dev.off()

