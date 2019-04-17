MetabolicHeatmap <- cbind(MetabolicPathway$logFC, MetabolicPathway$logFC)
PurineHeatmap <- cbind(PurineMetabolism$logFC, PurineMetabolism$logFC)
InsulinHeatmap <- cbind(InsulinSignaling$logFC, InsulinSignaling$logFC)
PhosphatidylHeatmap <- cbind(Phosphatidylinositol$logFC, Phosphatidylinositol$logFC)
PyrimindineHeatmap <- cbind(Pyrimindine$logFC, Pyrimindine$logFC)
AminoAcidSugarHeatmap <- cbind(AminoAcidSugar$logFC, AminoAcidSugar$logFC)
GlycolysisHeatmap <- cbind(Glycolysis$logFC, Glycolysis$logFC)
GlycerophospholipidHeatmap <- cbind(Glycerophospholipid$logFC, Glycerophospholipid$logFC)
FructoseHeatmap <- cbind(Fructose$logFC, Fructose$logFC)
InositolHeatmap <- cbind(Inositol$logFC, Inositol$logFC)

colorz = rev(brewer.pal(10,"RdBu")) 
colorz=colorRampPalette(colorz)(20)
anothername <- seq(-5, 5, 0.5)


pheatmap(MetabolicHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(PurineHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(InsulinHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(PhosphatidylHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(PyrimindineHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(AminoAcidSugarHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(GlycolysisHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(GlycerophospholipidHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(FructoseHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
pheatmap(InositolHeatmap, breaks = anothername, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0)
