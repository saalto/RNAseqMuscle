NewMetabolicTable <- NewMetabolicTable[order(NewMetabolicTable$logFC, decreasing = TRUE),]
NewOxidPhospTable <- NewOxidPhospTable[order(NewOxidPhospTable$logFC, decreasing = TRUE),]
NewHuntintonTable <- NewHuntintonTable[order(NewHuntintonTable$logFC, decreasing = TRUE),]
NewAlzheimerTable <- NewAlzheimerTable[order(NewAlzheimerTable$logFC, decreasing = TRUE),]
NewParkinsonTable <- NewParkinsonTable[order(NewParkinsonTable$logFC, decreasing = TRUE),]
NewPurineTable <- NewPurineTable[order(NewPurineTable$logFC, decreasing = TRUE),]
NewPyrimidineTable <- NewPyrimidineTable[order(NewPyrimidineTable$logFC, decreasing = TRUE),]
NewNGlycanTable <- NewNGlycanTable[order(NewNGlycanTable$logFC, decreasing = TRUE),]
NewPhosphatidyTable <- NewPhosphatidyTable[order(NewPhosphatidyTable$logFC, decreasing = TRUE),]
NewTCATable <- NewTCATable[order(NewTCATable$logFC, decreasing = TRUE),]

NewAlzheimerBind <- cbind(NewAlzheimerTable$logFC, NewAlzheimerTable$logFC)
NewHuntingtonBind <- cbind(NewHuntintonTable$logFC, NewHuntintonTable$logFC)
NewMetabolicBind <- cbind(NewMetabolicTable$logFC, NewMetabolicTable$logFC)
NewNGlycanBind <- cbind(NewNGlycanTable$logFC, NewNGlycanTable$logFC)
NewOxidativePhospBind <- cbind(NewOxidPhospTable$logFC, NewOxidPhospTable$logFC)
NewParkinsonBind <- cbind(NewParkinsonTable$logFC, NewParkinsonTable$logFC)
NewPhosphatidyBind <- cbind(NewPhosphatidyTable$logFC, NewPhosphatidyTable$logFC)
NewPurineTableBind <- cbind(NewPurineTable$logFC, NewPurineTable$logFC)
NewPyrimidineBind <- cbind(NewPyrimidineTable$logFC, NewPyrimidineTable$logFC)
NewTCABind <- cbind(NewTCATable$logFC, NewTCATable$logFC)

pheatmap(NewMetabolicBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewOxidativePhospBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewHuntingtonBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewAlzheimerBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewParkinsonBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewPurineTableBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewPyrimidineBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewNGlycanBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewPhosphatidyBind, color=colorz, cluster_rows=FALSE)
pheatmap(NewTCABind, color=colorz, cluster_rows=FALSE)

colorz = rev(brewer.pal(20,"RdBu")) 
colorz=colorRampPalette(colorz)(30)
 