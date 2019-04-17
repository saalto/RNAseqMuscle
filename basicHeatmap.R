# --- Run edgeR, obtain the edgeR table, then proceed ---

library(RColorBrewer)
library(pheatmap)

newdata <- deedg$table[which(deedge$table$Pvalue <0.05),]

SignGenes <- cbind(newdata$logFC, newdata$logFC)

colorz = rev(brewer.pal(10,"RdBu")) 
colorz=colorRampPalette(colorz)(20)  
pheatmap(SignGenes, color= colorz, cluster_rows = FALSE) 