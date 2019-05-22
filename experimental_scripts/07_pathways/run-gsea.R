# use clusterProfiler to do the GSEA
library(clusterProfiler)
library(org.Mm.eg.db)

# KEGG or GO?
enrich_type <- "GO"
# p value cutoff?
p_cut <- 0.05

# read in the table of sig. DE metabolic genes
#sig_met <- readRDS('/Users/kevinbrown/projects/pathway-analysis/results/2019.02.12/sig-met-table.rds')
sig_met <- AllSoleusWO6_DEMetGenes
# need the hash table from sym->entrez
mouse_hash <- readRDS('C:/Users/sarah/Dropbox/ME/Pitx2 Muscle RNASeq/R/mouse-hash.rds')

# ---- SHOULD NOT NEED TO CHANGE STUFF DOWN HERE ----

# drop rows of the table in which there is no entrez_id
pruned_sig <- sig_met[hash::has.key(row.names(sig_met),mouse_hash$stoe),]

# create the named vector
sig_genes <- pruned_sig$logFC
names(sig_genes) <- hash::values(mouse_hash$stoe,keys=row.names(pruned_sig),USE.NAMES=FALSE)

# sort by foldchange in decreasing order
sig_genes <- sort(sig_genes,decreasing=TRUE)

# run GSEA
if (enrich_type == "KEGG") {
  x <- gseKEGG(geneList=sig_genes,organism="mmu")
} else {
  x <- gseGO(geneList=sig_genes,OrgDb=org.Mm.eg.db)
}

# collect information into something usable
items_to_take <- which(attr(x,'result')$p.adjust < p_cut)
results <- data.frame("ID" = attr(x,'result')$ID[items_to_take], "Description" = attr(x,'result')$Description[items_to_take],
    "adjP" = attr(x,'result')$p.adjust[items_to_take],"setSize" = attr(x,'result')$setSize[items_to_take],stringsAsFactors=FALSE)

# also return lists of genes; list elements correspond to pathway IDs
genes <- list()
for (i in 1:length(results$ID)){
  # bust up the list of genes (which will be ENTREZ ids)
  id_entrez <- strsplit(attr(x,'result')$core_enrichment[i],"/")[[1]]
  genes[[results$ID[i]]] <- hash::values(mouse_hash$etos,keys=id_entrez,USE.NAMES=FALSE)
}
