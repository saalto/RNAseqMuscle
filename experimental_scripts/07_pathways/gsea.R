#' Gene Set Enrichment Analysis (using either KEGG or GO terms)
#' @param gene_frame: should be a data frame with gene symbols for rows, with at least one column (logFC)
#' @param entrez_hash: list containing two hashes: symbol-to-entrez (stoe) and entrez-to-symbol (etos)
#' @param organism: three-letter KEGG string code for organism
#' @param orgDb: AnnotationDBI to use to fetch terms
#' @param enrich_type: "GO" or "KEGG"
#' @param p_cut: (adjusted) p-value cutoff for results
#' @export
gsea <- function(gene_frame,entrez_hash,organism="mmu",orgDb=org.Mm.eg.db,enrich_type="GO",p_cut=0.05) {
  require(clusterProfiler)
  require(org.Mm.eg.db)
  # drop rows of the table in which there is no entrez_id
  pruned_sig <- gene_frame[hash::has.key(row.names(gene_frame),entrez_hash$stoe),]
  # create the named vector
  sig_genes <- pruned_sig$logFC
  names(sig_genes) <- hash::values(entrez_hash$stoe,keys=row.names(pruned_sig),USE.NAMES=FALSE)
  # sort by foldchange in decreasing order
  sig_genes <- sort(sig_genes,decreasing=TRUE)
  # run GSEA
  if (enrich_type == "KEGG") {
    x <- gseKEGG(geneList=sig_genes,organism=organism)
  } else {
    x <- gseGO(geneList=sig_genes,OrgDb=orgDb)
  }

  # if we didn't find anything, don't continue
  if (dim(attr(x,'result'))[1] == 0) {
    return(NULL)
  }

  # collect information into something usable
  results <- list()

  # data
  items_to_take <- which(attr(x,'result')$p.adjust < p_cut)
  dtable <- data.frame("ID" = attr(x,'result')$ID[items_to_take], "Description" = attr(x,'result')$Description[items_to_take],
      "adjP" = attr(x,'result')$p.adjust[items_to_take],"setSize" = attr(x,'result')$setSize[items_to_take],stringsAsFactors=FALSE)
  # list of genes keyed on term/pathway ID
  genes <- list()
  for (i in 1:length(dtable$ID)){
    # bust up the list of genes (which will be ENTREZ ids)
    id_entrez <- strsplit(attr(x,'result')$core_enrichment[i],"/")[[1]]
    genes[[dtable$ID[i]]] <- hash::values(entrez_hash$etos,keys=id_entrez,USE.NAMES=FALSE)
  }
  # collect
  results$dtable <- dtable
  results$genes <- genes

  return(results)
}
