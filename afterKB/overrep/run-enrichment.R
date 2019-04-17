# uses cluster profiler to do enrichment analysis (currently just KEGG)
library(clusterProfiler)

# KEGG or GO?
enrich_type <- "KEGG"
# organism?
organism <- "mmu"
# p value cutoff?
p_cut <- 0.05

# read in the table of sig. DE metabolic genes
#sig_met <- readRDS('/Users/kevinbrown/projects/pathway-analysis/results/2019.02.12/sig-met-table.rds')
sig_met <- AllSoleusWO6_DEMetGenes
# universe in this case should be all METABOLIC genes
#met_genes <- readRDS('/Users/kevinbrown/projects/pathway-analysis/R/metabolic-gene-list.rds')
met_genes <- MetabolicGenesFromWhole3
# need the hash table from sym->entrez
mouse_hash <- readRDS('C:/Users/sarah/Dropbox/ME/Pitx2 Muscle RNASeq/R/mouse-hash.rds')

# ---- SHOULD NOT NEED TO MAKE CHANGES DOWN HERE ----

# drop rows of the table in which there is no entrez_id
pruned_sig <- sig_met[hash::has.key(row.names(sig_met),mouse_hash$stoe),]

# make the entrez universe - this is dumb, but I can't figure out another way to do it
entrez_universe <- c()
for (i in 1:length(met_genes)){
  s <- met_genes[i]
  if (!is.null(mouse_hash$stoe[[s]])){
    entrez_universe <- c(entrez_universe,mouse_hash$stoe[[s]])
  }
}

if (enrich_type == "KEGG"){
  x <- enrichKEGG(hash::values(mouse_hash$stoe,keys=row.names(pruned_sig),USE.NAMES=FALSE),organism=organism,universe=entrez_universe)
}

# collect information into something usable
items_to_take <- which(attr(x,'result')$p.adjust < 0.05)
results <- data.frame("ID" = attr(x,'result')$ID[items_to_take], "Description" = attr(x,'result')$Description[items_to_take],
    "adjP" = attr(x,'result')$p.adjust[items_to_take],"Count" = attr(x,'result')$Count[items_to_take],stringsAsFactors=FALSE)

# also return lists of genes; list elements correspond to pathway IDs
genes <- list()
for (i in 1:length(results$ID)){
  # bust up the list of genes (which will be ENTREZ ids)
  id_entrez <- strsplit(attr(x,'result')$geneID[i],"/")[[1]]
  genes[[results$ID[i]]] <- hash::values(mouse_hash$etos,keys=id_entrez,USE.NAMES=FALSE)
}
