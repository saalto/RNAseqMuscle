source('C:/Users/sarah/Dropbox/ME/Pitx2 Muscle RNASeq/R/db-access/GO_genes.R')
library(org.Mm.eg.db)
library(GOplot)

# table of significant DE genes with info
# gene_frame <- readRDS('/Users/kevinbrown/projects/pathway-analysis/results/2019.02.12/sig-met-table.rds')
gene_frame <- AllSharedSoleus 

# sym-to-entrez and reverse hash
double_hash <- mouse_hash

# GO ids to assemble in the chord plot
#go_ids <- c("GO:0006807")

go_ids <- intersectSOlTb2$category

#Over-represented GO terms
overrep_file <- AllSoleusWO6_enrich

# filename to save chordplot
file_name <- 'C:/Users/sarah/Documents/chord-plot.pdf'

# ----- SHOULD NOT NEED TO CHANGE STUFF DOWN HERE ----

# set overall color scale (but make it symmetric) using max/min in table
max_fc <- max(gene_frame$logFC)
min_fc <- min(gene_frame$logFC)
if (max_fc > abs(min_fc)){
  fc_max <- max_fc
  fc_min <- -1.0*max_fc
} else {
  fc_max <- -1.0*min_fc
  fc_min <- min_fc
}

# assemble the data
#goplot_data <- GO_assemble(gene_frame,double_hash,overrep_file, go_ids,annot ="org.Mm.eg.db")
# make the plot
pdf(file_name, width = 20, height = 20)
pl <- GOChord(goplot_data$chord,space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('red','white','blue'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
pl <- pl + coord_fixed() + theme(legend.position = "bottom")
print(pl)
dev.off()
