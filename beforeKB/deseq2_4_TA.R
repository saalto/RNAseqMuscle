## This script uses files that were (1) renamed after identifying 
## mislabeling of tibialis anterior samples.and (2) removal of one 
## file that the heatmap showed confusing expression levels when 
## compared to the other mutants.

##---Clear data and load packages-----------------------
rm(list = ls())
dev.off()

source('https://bioconductor.org/biocLite.R')
library("DESeq2")
library("ggplot2")
library("ComplexHeatmap")
library("pamr")
library("MCL")

##---Set working directory to iteration 2---
setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/renamed_files/")
directory <- "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/renamed_files/"

##---Set up DESeq2 data, based on names of HTSeq counts in working directory---
sampleFiles <- dir(pattern = 'sorted')
#print(sampleFiles)
sampleIdentifiers <- c("C1", "C2", "M1", "M2", "M3","M4", "M5")

#---sample group set up---
ConditionMatch <- regexpr(pattern = '[A-Z]+', dir(pattern = '.txt'))
#print(ConditionMatch)
sampleConditions <- regmatches(dir(pattern = '*.txt'), ConditionMatch)
#print(sampleConditions)
sampleTable <- data.frame(sampleName = sampleIdentifiers, fileName = sampleFiles, 
                          condition = sampleConditions)
#print(sampleTable)

##---Calculate DESeq2 from HTSeq count tables-------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, 
                                       design = ~ condition)
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

#-Plot the read count after rlog transformation
dists <-dist(t(logTransCounts))
plot(hclust(dists))
#head(logTransCounts)

#---Singling out specific genes based on log counts---------------------
#logTransCounts[grep("Arf1", rownames(logTransCounts)), ]
logTransCounts[grep("Pitx2", rownames(logTransCounts)), ] 
#Wildtype and mutant samples have similar counts
logTransCounts[grep("Foxo3", rownames(logTransCounts)), ] 
#Mutant samples will have higher counts than wildtypes
logTransCounts[grep("Lars2", rownames(logTransCounts)), ] 
#Mutants sample have higher counts than controls


#---Comparing heatmaps of the three normalizing methods------------------
library("RColorBrewer")
#install.packages("gplots")
library("gplots")

select <- order(rowMeans(counts(ddsHTSeqFiltered,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


#pdf("readcountsTranform.pdf")
heatmap.2(counts(ddsHTSeqFiltered,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6), 
          main = "Read Counts Transformation")

dev.off()

#pdf("rlogTransform.pdf")
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6), 
          main = "rLog Transformation")

dev.off()

#pdf("VariStablizeTransform.pdf")
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6), 
          main = "Variance Stablizing Transformation")

dev.off()


##---Comparing samples to each other for correlation patterning
#-using rlog transformation (rld)---
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <-  with(colData(ddsHTSeqFiltered),
                                        paste(condition, type, sep = " : "))
heatmap.2(mat, trace = "none", col = rev(hmcol), margins = c(13,13))
          #main = "Correlation between Samples based on rLog Transformation")

dev.off()

#-using variance stabilizing transformation (VSD)---
distsVSD <- dist(t(assay(vsd)))
mat <- as.matrix(distsVSD)
rownames(mat) <- colnames(mat) <-  with(colData(ddsHTSeqFiltered),
                                        paste(condition, type, sep = " : "))
heatmap.2(mat, trace = "none", col = rev(hmcol), margins = c(13,13))
          #main = "Correlation between Samples based on Variance Stabilizing Transformation")

dev.off()


##---MDS plot using Euclidean distances
#-rLog Transformation (rld)----
distsRL <- dist(t(logTransCounts))
DistMatrix <- as.matrix(distsRL)
mdsData <- data.frame(cmdscale(DistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
mds$condition <- c("Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant")
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3)
  #ggtitle("MDS using Euclidean distance and rLog Transformation")

dev.off()

#-Variance Stablizing Transformation (vsd)----
distsVSD <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(distsVSD)
mdsData <- data.frame(cmdscale(DistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
mds$condition <- c("Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant")
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3)
  #ggtitle("MDS using Euclidean distances and Variance Stabilizing Transformation")

dev.off()


##---Poisson Distance Plot------------
#install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsHTSeqFiltered)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(ddsHTSeqFiltered)))
mdsPois$condition <- c("Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant")
ggplot(mdsPois, aes(X1,X2,color=condition)) + geom_point(size=3)
  #ggtitle("Poisson Distance Plot of the Read Counts")

dev.off()


##---Principle component analysis based on log2 normalized count matrix---
# Load factoextra for visualization
#install.packages("factoextra")
library(factoextra)

# Compute PCA
OGPCAN <- prcomp(logTransCounts, center = T, scale = F, tol = 0)
#print(OGPCAN)

# Visualize eigenvalues (scree plot).
#Show the percentage of variance explained by each principal component.
fviz_eig(OGPCAN)

# Eigenvalues
eig.val <- get_eigenvalue(OGPCAN)
#head(eig.val)

# Results for individuals
res.ind <- get_pca_ind(OGPCAN)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation

# Graph of individuals. Individuals with a similar profile are grouped together
fviz_pca_ind(OGPCAN,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # text overlapping
)

# Results for Variables
res.var <- get_pca_var(OGPCAN)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation

grp <- c("Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant")

#Graph of variables.
# Positive correlated variables point to the same side of the plot.
# Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(OGPCAN,
             col.var = grp, # Color by contributions to the PC
             palette = c("#00AFBB", "#FC4E07"),
             legend.title = "Tibialis Anterior",
             repel = TRUE    # text overlapping
)

# Biplot of individuals and variables
fviz_pca_biplot(OGPCAN, repel = FALSE, arrowsize =2,
                col.ind = "#696969",  # Individuals color,
                col.var = grp, 
                legend.title="Tibialis Anterior", 
                title = NULL,
                palette = c("#0000ff", "#00ffff")
                ) + scale_y_reverse()

# PCA plot replicates set up
OGPCAN_matrix <- as.data.frame(OGPCAN$rotation)
#print(OGPCAN_matrix)
OGPCAN_matrix$Condition <- c("Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant")
#print(OGPCAN_matrix)

# Plot PCA
ggplot(OGPCAN_matrix, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3) + geom_text(aes(label=row.names(OGPCAN_matrix)),vjust=0, hjust=0)+
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size  = 16, face = "bold"),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14)) +
  scale_color_discrete(name = "Sample") +
  xlab(paste0("PC1: ", sprintf("%.3f", eig.val$variance.percent[1]), "% variance")) +
  ylab(paste0("PC2: ", sprintf("%.3f", eig.val$variance.percent[2]), "% variance"))
#ggtitle("Principle Component Analysis based on rlog transformation")

dev.off()


##---Calculate differentially expressed genes from DESeq2---------------------
res.TA_W_M <- results(ddsHTSeqFiltered, contrast = c("condition", "TM", "TF"))
#summary(res.TA_W_M)

##---MA plot from results---------------------------------------
#-Bland-Altman plot that visualizes the differences between measurements 
# taken in two samples, by transforming the data onto M (log ratio) and 
# A (mean average) scales, then plotting these values
plotMA.TA_W_M <- plotMA(res.TA_W_M, ylim=c(-2,2))

##---Volcano Plots----------------------------------------------------------------
with(res.TA_W_M, plot(log2FoldChange, -log10(pvalue), 
                      pch=10,
                      xlab=expression(paste("log"[2]*Delta,"FC (M/C)")), 
                      ylab=expression(paste("-log"[10]*"(p-value)"))))

#-Add colored points: blue if padj<0.05, red if log2FC>1, green if both
with(subset(res.TA_W_M, padj<0.05), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="light grey"))
with(subset(res.TA_W_M, abs(log2FoldChange)>1), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="dark grey"))
with(subset(res.TA_W_M, padj<0.05 & log2FoldChange >1), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="yellow"))
with(subset(res.TA_W_M, padj<0.05 & log2FoldChange < (-1)), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

#-Label points with the textxy function from the calibrate plot
with(subset(res.TA_W_M),
     identify(log2FoldChange, -log10(pvalue), labels=rownames(res.TA_W_M), cex=0.6))

dev.off()


##---Filtering the Results (DE genes) into tables---
#-Subset data based on (1) adjusted p-value less than 0.05 AND 
# (2) absolute value of the log2 fold change greater than 0.5
res.TA_W_M_filtered2 <- subset(res.TA_W_M, padj < 0.05)
res.TA_W_M_filtered2$absFC <- abs(res.TA_W_M_filtered2$log2FoldChange)
#head(res.TA_W_M_filtered2)
#dim(res.TA_W_M_filtered2)

res.TA_W_M_filtered3 <- subset(res.TA_W_M_filtered2, absFC > 1)
res.TA_W_M_filtered4 <- subset(res.TA_W_M, pvalue < 0.05)
#head(res.TA_W_M_filtered3)
#dim(res.TA_W_M_filtered3)

#---Print out the filtered data as a text file---
write.table(res.TA_W_M_filtered2, 
            "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis/res.TA_W_M_filtered_padj_20180918.txt", 
            sep ="\t")
write.table(res.TA_W_M_filtered3, 
            "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis/res.TA_W_M_filtered_padjfoldchange_20180918.txt", 
            sep ="\t")
write.table(res.TA_W_M, 
            "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis/res.TA_W_M_20180918.txt", 
            sep = "\t")

write.table(res.TA_W_M_filtered4, 
            "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/RNAseq_analysis/res.TA_W_M_filtered_pvalue_20180928.txt", 
            sep ="\t")

##---Heatmap of the most significant fold-change genes--------------
#install.packages("pheatmap")
library("grid")
library("pheatmap")

# Edit body of pheatmap:::draw_colnames, customizing it to your liking
# For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, 
                 hjust = 0, rot = 0, gp = gpar(...))
  return(res)}

# 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# Heatmap of the significant (padj<0.05) DE genes based on VSD
Mat <- assay(vsd)[order(res.TA_W_M_filtered2$padj), ]
Mat <- Mat - rowMeans(Mat)
df <- as.data.frame(colData(vsd)[,c("condition")])
pheatmap(Mat, color= colorRampPalette(c("#0000ff", "#000000", "#ffff00"))(5), 
         breaks = c(-2, -1, -0.2, 0.2, 1, 2), cluster_col= F, 
         treeheight_row = 0, fontsize = 15,
         show_rownames = F, show_colnames = T)

dev.off()

## Heatmap of the significant (padj<0.05) DE genes based on VSD
pheatmap(Mat,
         color= colorRampPalette(c("#FF0000", "#000000", "#00ff00"))(100), 
         #breaks = c(-2, -1,-0.2, 0.2, 1, 2),
         cluster_col = F, show_rownames = F, show_colnames = T)

##---Clear data and delete figure----------------------
rm(list = ls())
  dev.off()