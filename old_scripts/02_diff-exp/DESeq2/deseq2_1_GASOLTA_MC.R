## This script uses files from gastrocnemius, soleus, and tibialis
## that were (1) renamed after identifying mislabeled samples, 
## (2) removed if expression levels were confusing compared to 
## the other mutants.

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
setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/combined/")
directory <- "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/combined/"

##---Set up DESeq2 data, based on names of HTSeq counts in working directory---
sampleFiles <- dir(pattern = 'sorted')
print(sampleFiles)
sampleIdentifiers <- c("CT1", "CT2","CG1", "CG2", "MT1", "MT2", "MT3", "MT4", "MT5", "MG1", "MG2")

#---sample group set up---
ConditionMatch <- regexpr(pattern = '[A-Z]+', dir(pattern = '.txt'))
print(ConditionMatch)
sampleConditions <- regmatches(dir(pattern = '*.txt'), ConditionMatch)
print(sampleConditions)
sampleTable <- data.frame(sampleName = sampleIdentifiers, fileName = sampleFiles, 
                          condition = sampleConditions)
print(sampleTable)

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
mds$condition <- c("Control", "Control", "Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant")
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3)
#ggtitle("MDS using Euclidean distance and rLog Transformation")
dev.off()

#-Variance Stablizing Transformation (vsd)----
distsVSD <- dist(t(assay(vsd)))
DistMatrix <- as.matrix(distsVSD)
mdsData <- data.frame(cmdscale(DistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
mds$condition <- c("Control", "Control", "Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant")
ggplot(mds, aes (X1, X2, color=condition)) + geom_point(size=3) + geom_text(aes(label=row.names(OGPCAN_matrix)),vjust=0, hjust=0)
#ggtitle("MDS using Euclidean distances and Variance Stabilizing Transformation")
dev.off()


##---Poisson Distance Plot------------
#install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsHTSeqFiltered)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(ddsHTSeqFiltered)))
mdsPois$condition <- c("Control", "Control", "Control", "Control", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant")
ggplot(mdsPois, aes(X1,X2,color=condition)) + geom_point(size=3) + geom_text(aes(label=row.names(OGPCAN_matrix)),vjust=0, hjust=0)

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

grp <- c("GA Control", "GA Control", "GA Mutant", "GA Mutant", 
         "SOL Control", "SOL Control", "SOL Control", "SOL Control",
         "SOL Control", "SOL Control", "SOL Mutant", "SOL Mutant",
         "SOL Mutant", "SOL Mutant","SOL Mutant", "TA Control",
         "TA Control", "TA Mutant","TA Mutant","TA Mutant","TA Mutant")

#Graph of variables.
# Positive correlated variables point to the same side of the plot.
# Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(OGPCAN,
             col.var = grp, # Color by contributions to the PC
             palette = c("#FC4E07", "#00AFBB", "#FC4E07", "#00AFBB", 
                         "#FC4E07", "#00AFBB", "#00AFBB"),
             repel = TRUE    # text overlapping
)

# Biplot of individuals and variables
fviz_pca_biplot(OGPCAN, repel = FALSE, arrowsize =2,
                col.ind = "#696969",  # Individuals color,
                col.var = grp, 
                legend.title="Tissue Sample", 
                title = NULL,
                palette = c("#FC4E07", "#00AFBB", "#FC4E07", "#00AFBB", 
                            "#FC4E07", "#00AFBB", "#00AFBB"))

# PCA plot replicates set up
OGPCAN_matrix <- as.data.frame(OGPCAN$rotation)
#print(OGPCAN_matrix)
OGPCAN_matrix$Condition <- c("TA Control", "TA Control", "GA Control", "GA Control",
                             "TA Mutant","TA Mutant","TA Mutant","TA Mutant", "TA Mutant","GA Mutant", "GA Mutant")
#print(OGPCAN_matrix)

# Plot PCA
ggplot(OGPCAN_matrix, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size  = 16, face = "bold"),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14)) +
  scale_color_discrete(name = "Sample") +
  xlab(paste0("PC1: ", sprintf("%.3f", eig.val$variance.percent[1]), "% variance")) +
  ylab(paste0("PC2: ", sprintf("%.3f", eig.val$variance.percent[2]), "% variance")) +
  stat_ellipse(type= "t")+ 
  stat_ellipse(type="norm", linetype = 2)
#ggtitle("Principle Component Analysis based on rlog transformation")

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

##---Clear data and delete figure-----------------------
rm(list = ls())
dev.off()