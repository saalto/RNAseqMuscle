source('C:/Users/sarah/Dropbox/ME/Pitx2 Muscle RNASeq/R/db-access/KEGG_genes.R')
library(org.Mm.eg.db)
library(GOplot)

# table of significant DE genes with info
# gene_frame <- readRDS('/Users/kevinbrown/projects/pathway-analysis/results/2019.02.12/sig-met-table.rds')
gene_frame <- QLT_edgeRSoleusSignMet

# sym-to-entrez and reverse hash
double_hash <- readRDS('C:/Users/sarah/Dropbox/ME/Pitx2 Muscle RNASeq/R/mouse-hash.rds')

# KEGG id -> name hash
kegg_hash <- readRDS('C:/Users/sarah/Dropbox/ME/Pitx2 Muscle RNASeq/R/kegg-hash.rds')

# KEGG ids to assemble in the chord plot
# Nucleotide graphic
#kegg_ids <- c("00230", "00240")
# PIs graphic
#kegg_ids <- c("04910", "04070", "00564", "00562")
# Sugars graphic
kegg_ids <- c("00520","00010", "00051")


# filename to save chordplot
file_name <- 'C:/Users/sarah/OneDrive/Documents/2019/01_2019_Winter/PROJECT/RNASeq_Analysis/chord-plot10.pdf'

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
goplot_data <- KEGG_assemble(gene_frame,double_hash,kegg_hash,kegg_ids)
# make the plot
pdf(file_name, width = 11, height = 10)
pl <- GOChord(goplot_data$chord,space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('blue','white','red'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
pl <- pl + coord_fixed()
print(pl)
dev.off()
