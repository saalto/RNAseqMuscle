# --- Heat map of Soleus data set renamed and without Sample C6 due to low z-scores---
%AllSoleusWO6_ReadCounts2 <- AllSoleusWO6_ReadCounts[row.names(AllSoleusWO6_ReadCounts) %in% MetabolicGenesFromWhole3,]
df1 <- AllSoleusWO6_DEMet_ReadCounts[apply(AllSoleusWO6_DEMet_ReadCounts,1,prod)>0,]
df2 <- log(df1)
AllSoleusWO6_zscore <- apply(df2,1,scale)
breakList = seq(min(AllSoleusWO6_zscore),max(AllSoleusWO6_zscore), by=0.2)

library(pheatmap)
library(RColorBrewer)
colorz = rev(brewer.pal(10,"RdBu"))
colorz = colorRampPalette(colorz)(length(breakList))
pheatmap(AllSoleusWO6_zscore, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)

# --- Heat map of TA data set renamed ---
AllTibialisReadCounts5 <- read.delim(file = "C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/tibialis/renamed_files/tibialis-data-combined1.txt", header=FALSE, sep = ",")
# View(AllTibialisReadCounts5)
#AllTibialisReadCounts6 <- cbind(C1=AllTibialisReadCounts5$V2, C2=AllTibialisReadCounts5$V3, M1=AllTibialisReadCounts5$V4, M2=AllTibialisReadCounts5$V5, M3=AllTibialisReadCounts5$V6, M4=AllTibialisReadCounts5$V7, M5=AllTibialisReadCounts5$V8)
#row.names(AllTibialisReadCounts6) <- AllTibialisReadCounts5$V1
#df4 <- AllTibialisReadCounts6[row.names(AllTibialisReadCounts6) %in% MetabolicGenesFromWhole3,]
df5 <- AllTibialis_renamed_DEMet_ReadCounts[apply(AllTibialis_renamed_DEMet_ReadCounts,1,prod)>0,]
df6 <- log(df5)
zscores6 <- apply(df6,1,scale)
#breakList2 = seq(min(zscores6),max(zscores6), by=0.2)
pheatmap(zscores6, breaks = breakList,color = colorz, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)

# --- Heat map of Ga data set ---
AllGastroReadCounts5 <- read.delim(file = "C:/Users/sarah/OneDrive/Documents/2019/01_2019_Winter/PROJECT/RNASeq_Analysis/2019_03_28/Ga/AllGastro_ReadCounts.txt", header=FALSE, sep = ",")
AllGastroReadCounts6 <- cbind(C1=AllGastroReadCounts5$V2, C2=AllGastroReadCounts5$V3, M1=AllGastroReadCounts5$V4, M2=AllGastroReadCounts5$V5, M3=AllGastroReadCounts5$V6)
row.names(AllGastroReadCounts6) <- AllGastroReadCounts5$V1
df7 <- AllGastroReadCounts6[row.names(AllGastroReadCounts6) %in% MetabolicGenesFromWhole3,]
df8 <- AllGastro_DEMet_ReadCounts[apply(AllGastro_DEMet_ReadCounts,1,prod)>0,]
df9 <- log(df8)
df10 <- df9[apply(df9,1,prod)>1,]
zscores7 <- apply(df10,1,scale)
breakList3 = seq(min(zscores7),max(zscores7), by=0.2)
pheatmap(zscores7, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)


# --- Heat map of Ga data set renamed---
#AllGastroReadCounts5 <- read.delim(file = "C:/Users/sarah/OneDrive/Documents/2019/01_2019_Winter/PROJECT/RNASeq_Analysis/2019_03_28/Ga/AllGastro_renamed_ReadCounts.txt", header=FALSE, sep = ",")
#AllGastroReadCounts6 <- cbind(C1=AllGastroReadCounts5$V2, C2=AllGastroReadCounts5$V3, C3=AllGastroReadCounts5$V4, M1=AllGastroReadCounts5$V5, M2=AllGastroReadCounts5$V6)
#row.names(AllGastroReadCounts6) <- AllGastroReadCounts5$V1
#df7 <- AllGastroReadCounts6[row.names(AllGastroReadCounts6) %in% MetabolicGenesFromWhole3,]
df8 <- AllGastro_renamedremoved_DEMet_ReadCounts[apply(AllGastro_renamedremoved_DEMet_ReadCounts,1,prod)>0,]
df9 <- log(df8)
df10 <- df9[apply(df9,1,prod)>1,]
zscores7 <- apply(df10,1,scale)
breakList3 = seq(min(zscores7),max(zscores7), by=0.2)
pheatmap(zscores7, breaks= breakList, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)



AllGastroWOM2_catenrich <- catenrich(all.genes = avector, sig.genes = row.names(AllGastroWOM2_pantherFilter), entrez_hash = mouse_hash, kegg_hash = kegg_hash, enrich_type = "GO:BP", p_cut = 0.05, genome = "mm10")
AllTibialisWOM4_catenrich <- catenrich(all.genes = avector, sig.genes = row.names(AllTibialisWOM4_pantherFilter), entrez_hash = mouse_hash, kegg_hash = kegg_hash, enrich_type = "GO:BP", p_cut = 0.05, genome = "mm10")
AllTibialisWOM4_catenrich_shorten <- AllTibialisWOM4_catenrich$overrep[AllTibialisWOM4_catenrich$overrep$category %in% go_term_list,]
AllGastroWOM2_catenrich_shorten <- AllGastroWOM2_catenrich$overrep[AllGastroWOM2_catenrich$overrep$category %in% go_term_list,]

df1 <- AllSoleusReadCounts[row.names(AllSoleusReadCounts) %in% AllSoleus_SSMetGenes,]
df2 <- df1[apply(df1,1,prod)>0,]
df3 <- log(df2)
zscore1 <- apply(df3,1,scale)
pheatmap(zscore1, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)
