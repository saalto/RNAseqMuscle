#
##-------Clear data and figure & load packages-----------------------
rm(list = ls())
dev.off()
source('https://bioconductor.org/biocLite.R')
library("DESeq2")
library("ggplot2")
library("ComplexHeatmap")
library("pamr")
library("MCL")
library("AnnotationDbi")
library("Mus.musculus")

##---Set working directory to iteration 1-----------------------------
#setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration1/RNAseq_SOLtextfiles")
#directory <- "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration1/RNAseq_SOLtextfiles/"

##---Set working directory to iteration 2-----------------------------
setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus")
directory <- "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus/"

##---Set up DESeq2 data, based on file names of HTSeq counts in working directory---
sampleFiles <- dir(pattern = 'sorted')
#print(sampleFiles)

#---sample group set up---
ConditionMatch <- regexpr(pattern = '[A-Z]+', dir(pattern = '.txt'))
#print(ConditionMatch)
sampleConditions <- regmatches(dir(pattern = '*.txt'), ConditionMatch)
#print(sampleConditions)
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleConditions)
#print(sampleTable)

##----ONLY USE IF USING ARUN'S TEXT FILES LOCATED IN ITERATION 1----------
# #-removal of files referring to GA and TA
# sampleTable <- sampleTable[c(1:7, 9:13),] #Use only with Arun's text files
# #print(sampleTable)

# # -reassignment of the files based on the Foxo3 gene expression patterning,
# # -which is higher in the MCK mice than in the flox/Z mice. this gene expression
# # -is a better indicator of genotype than the labeling provided
# # -DESeq2 object must be created and rld array generated in order to observe
# # -Foxo3 gene expressions.
# # Use only with Arun's text files
# sampleTable$condition <- c(rep("MCK", 3), rep("FLZ", 4), rep("MCK", 2), rep("FLZ", 3))
# print(sampleTable)

##---Calculate DESeq2 from HTSeq count tables-------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)
#print(ddsHTSeq)

ddsHTSeqFiltered <- ddsHTSeq [ rowSums(counts(ddsHTSeq)) > 0, ] #Filtering out genes with zero counts
ddsHTSeqFiltered <- DESeq(ddsHTSeqFiltered)
#print(ddsHTSeqFiltered)

##---Visualizations of Data Transformations-------------------------------------
#-Generate log2 normalized count matrix and variance stabilizing matrix
rld <- rlog(ddsHTSeqFiltered, blind = FALSE)
vsd <- varianceStabilizingTransformation(ddsHTSeqFiltered, blind = FALSE)
logTransCounts <- assay(rld)
#head(logTransCounts)

#---Singling out specific genes: Use this section to confirm the expression levels 
# of known genes and determine if samples were mislabeled-
logTransCounts[grep("Pitx2", rownames(logTransCounts)), ] #expression levels should be similar between wildtype and mutant
logTransCounts[grep("Foxo3", rownames(logTransCounts)), ] #expression levels should be different between wildtype and mutant
#logTransCounts[grep("Arf1", rownames(logTransCounts)), ]

 
##---Comparing heatmaps of the three normalizing methods------------------
#TO-DO: Print out the figures
library("RColorBrewer")
#install.packages("gplots")
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


##---Comparing samples to each other for correlation patterning 
#-using rlog transformation (rld)---
#TO-DO: print off the figures
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <-  with(colData(ddsHTSeqFiltered), paste(condition, type, sep = " : "))
heatmap.2(mat, trace = "none", col = rev(hmcol), margins = c(13,13), main = "Correlation between Samples based on rLog Transformation")
dev.off()

#-using variance stabilizing transformation (VSD)---
#TO-DO: print off the figures
distsVSD <- dist(t(assay(vsd)))
mat <- as.matrix(distsVSD)
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
distsVSD <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(distsVSD)
mdsData <- data.frame(cmdscale(DistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3) + ggtitle("MDS using Euclidean distances and Variance Stabilizing Transformation")
dev.off()

##---Poisson Distance Plot------------
#TO-DO: print off the figures
#install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsHTSeqFiltered)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(ddsHTSeqFiltered)))
ggplot(mdsPois, aes(X1,X2,color=condition)) + geom_point(size=3) + ggtitle("Poisson Distance Plot of the Read Counts")
dev.off()


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
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Principle Component Analysis based on rlog transformation")
#print(ddsHTSeqFiltered)
dev.off()


##---Calculate differentially expressed genes from DESeq2---
res.SOL_W_M <- results(ddsHTSeqFiltered, contrast = c("condition", "M", "F"))
#summary(res.SOL_W_M)

##---MA plot of results---
#TO-DO: print off the figures
plotMA.SOL_W_M <- plotMA(res.SOL_W_M, ylim=c(-2,2))
mcols(res.SOL_W_M, use.names = TRUE)
dev.off()

##---Volcano Plot of the results----
with(res.SOL_W_M, plot(log2FoldChange, -log10(pvalue), pch=10, main="Volcano plot", xlim=c(-5,7)))

#-Add colored points: red if padj<0.05, orange if log2FC>1, green if both
with(subset(res.SOL_W_M, padj<0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.SOL_W_M, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res.SOL_W_M, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#-Label points with the textxy function from the calibrate plot
with(subset(res.SOL_W_M), identify(log2FoldChange, -log10(pvalue), labels=rownames(res.SOL_W_M))) #Need to click on graphic to label the outliers
dev.off()


##---Flipped Volcano Plot of results---
#TO-DO: flip the volcano plot, change the colors
#-Measuring the effect of fold change and the statistical significance
with(res.SOL_W_M, plot(log2FoldChange, log10(pvalue), pch=10, main="Volcano plot", xlim=c(-5,7)))

# #-Add colored points: red if padj<0.05, orange if log2FC>1, green if both
# with(subset(res.SOL_W_M, padj<0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(res.SOL_W_M, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
# with(subset(res.SOL_W_M, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
# 
# #-Label points with the textxy function from the calibrate plot
# with(subset(res.SOL_W_M), identify(log2FoldChange, -log10(pvalue), labels=rownames(res.SOL_W_M))) #Need to click on graphic to label the outliers
# dev.off()


##---Filtering the Results (DE genes) into tables---
#-Calculate differentially expressed genes from DESeq2 object based on adjusted p-value
res.SOL_W_M.05 <- results(ddsHTSeqFiltered, alpha=0.05)
table(res.TA_W_M.05$padj < 0.05)

#-Calculate differentially expressed genes from DESeq2 object based on log fold change equal to 0.5 (2^0.5)
res.SOL_W_MLFC1 <- results(ddsHTSeqFiltered, lfcThreshold=0.5)
table(res.TA_W_MLFC1$padj < 0.05)

#-Subset data based on (1) adjusted p-value less than 0.05 (2) absolute value of the log2 fold change greater than 0.5
res.SOL_W_M_filtered2 <- subset(res.SOL_W_M, padj < 0.05)
res.SOL_W_M_filtered2$absFC <- abs(res.SOL_W_M_filtered2$log2FoldChange)
#head(res.SOL_W_M_filtered3)
#dim(res.SOL_W_M_filtered3)

res.SOL_W_M_filtered3 <- subset(res.SOL_W_M_filtered2, absFC > 0.5)
#head(res.SOL_W_M_filtered2)
#dim(res.SOL_W_M_filtered2)

#-Print out the filtered data as a text file
write.table(res.SOL_W_M_filtered2, "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis_tables/res.SOL_W_M_filtered_padj_20180821.txt", sep ="\t")
write.table(res.SOL_W_M_filtered3, "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis_tables/res.SOL_W_M_filtered_padjfoldchange_20180821.txt", sep ="\t")
write.table(res.SOL_W_M, "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis_tables/res.SOL_W_M_Nofilter_20180820.txt", sep = "\t")


##---Heatmap of the most significant fold-change genes--------------
#install.packages("pheatmap")
library("pheatmap")

#Heatmap of the significant (padj<0.05) DE genes based on VSD
#library(RColorBrewer)
Mat <- assay(vsd)[order(res.SOL_W_M_filtered2$padj), ]
Mat <- Mat - rowMeans(Mat)
df <- as.data.frame(colData(vsd)[,c("condition")])
pheatmap(Mat, color= colorRampPalette(c("#990066", "#ffffff", "#009933"))(14), show_rownames = F, show_colnames = F)
dev.off()

# #Heatmap of the significant (padj<0.05) DE genes based on RLD
# Map <- assay(rld)[order(res.SOL_W_M_filtered2$padj), ]
# Map <- Map - rowMeans(Map)
# df <- as.data.frame(colData(rld)[,c("condition")])
# pheatmap(Map)
# dev.off()


##----for Heatmap plots of 200 gene expressions and significant Log2fold changes-----
# TO DO: Work in progress...I am not seeing what I should be seeing
# var_genes <- apply(logTransCounts, 1, var)
# #head(var_genes)
# select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]
# #head(select_var)
# highly_variable_lcpm <- logTransCounts[select_var,]
# #dim(highly_variable_lcpm)
# #head(highly_variable_lcpm)
# Heatmap(highly_variable_lcpm)
# dev.off()
# 
# # make a matrix and heatmap of the most significant and filtered out genes for each individual
# # NOTE: I believe I am comparing the genes twice. I do not trust this code completely.
# var_genes <- apply(res.SOL_W_M_filtered2, 1, var)
# #head(var_genes)
# select_var <- names(sort(var_genes))
# #head(select_var)
# highly_variable_lcpm <- logTransCounts[select_var,]
# #dim(highly_variable_lcpm)
# #head(highly_variable_lcpm)
# Heatmap(highly_variable_lcpm)
# dev.off()


##---Annotating and exporting Results using edgeR and DESeq objects---
#---NOTE: mm10 reference genome does not have keys-
# require(Mus.musculus)
# columns(Mus.musculus)
# keytypes(Mus.musculus)
# keys(Mus.musculus,"ENSEMBL")
# #TO-DO: the gene names need to match the keys
# keys <- row.names(res.SOL_W_M_filtered2)
# select(Mus.musculus, keys, column = "SYMBOL", keytype = "ENSEMBL")
# res.SOL_W_M_filtered2$symbol <- mapIds(Mus.musculus, 
#                                        keys = row.names(res.SOL_W_M_filtered2),
#                                        column ="SYMBOL",
#                                        keytype = "ENSEMBL",
#                                        multiVals = "first")
# ## 'select()' returned 1:many mapping between keys and columns
# #$y$genes$symbol <- res.SOL_W_M_filtered2$symbol #For EdgeR
# 
# res.SOL_W_M_filtered2$entrez <- mapIds(Mus.musculus,
#                                        keys = row.names(res.SOL_W_M_filtered2),
#                                        column = "ENTREZID",
#                                        keytype = "ENSEMBL",
#                                        multiVals = "first")
# ## 'select()' returned 1:many mapping between keys and columns
# #y$genes$entrez <- res.SOL_W_M_filtered2$entrez #For EdgeR
# 
# res.SOL_W_M_filtered2$symbol <- mapIds(Homo.sapiens,
#                                        keys=row.names(res.SOL_W_M_filtered2),
#                                        column="GENENAME",
#                                        keytype="ENSEMBL",
#                                        multiVals="first")
# ## 'select()' returned 1:many mapping between keys and columns
# #y$genes$symbol <- res.SOL_W_M_filtered2$symbol #For EdgeR
# 
# resOrdered <- res.SOL_W_M_filtered2[order(res.SOL_W_M_filtered2$padj),]
# head(resOrdered)

##---Annotation Using DESeq object---
#---NOTE: mm10 reference genome does not have keys---
# install.packages("goseq")
# library(goseq)
# 
# assayed.genes <- rownames(res.SOL_W_M)
# de.genes <-rownames(res.SOL_W_M)[which(res.SOL_W_M$padj < 0.05)]
# gene.vector = as.integer(assayed.genes%in%de.genes)
# names(gene.vector) = assayed.genes
# head(gene.vector)
# 
# supportedOrganisms()
# supportedOrganisms()[supportedOrganisms()$Genome=="mm10",]
# pwf = nullp(genes, "mm10", "ensGene")
# pwf = nullp(genes, "mm10", "knownGene")
# pwf = nullp(genes, "mm10", "geneSymbol")
# head(pwf)


##---Perform Differential Correlation Analysis---
install.packages("DGCA")
library(DGCA)
#data(Gene names as row names and sample names in the columns, gene expression value in the corresponding cell)
#data(DESeq object to identify sample names to tissue type)

# ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
#                      compare = c("F", "M"),
#                      adjust = "none", heatmapPlot = TRUE, nPerm = 0, nPairs = 100)
# head(ddcor_res)


#---Clear data and load packages-----------------------
rm(list = ls())
dev.off()
