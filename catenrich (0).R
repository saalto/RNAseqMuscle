##' Term enrichment counting (using goseq) including corrections for gene length bias
##' @param all.genes: the gene universe (names only) to test against
##' @param sig.genes: the DE genes to test
##' @param entrez_hash: a double hash giving symbol-to-entrez ($stoe) and vice versa ($etos)
##' @param kegg_hash: hash for KEGG ids -> descriptors (not used if asking for GO enrichmetn)
##' @param enrich_type: "GO:BP", "GO:MF", "GO:CC" or "KEGG"
##' @param p_cut: (adjusted) p-value cutoff for results
##' @param genome: name of genome to test against
##' @export
catenrich <- function(all.genes,sig.genes,entrez_hash,kegg_hash,enrich_type="KEGG",p_cut=0.05,genome="mm10") {
  require(goseq)
  cats <- c(enrich_type)

  # remove any genes without entrez IDs
  all.withid <- all.genes[hash::has.key(all.genes,entrez_hash$stoe)]
  sig.withid <- sig.genes[hash::has.key(sig.genes,entrez_hash$stoe)]

  # remove the DE genes from all the genes (if they are in there)
  pruned.withid <- all.withid[!all.withid %in% sig.withid]

  # make sure there are no repeated entries in the pruned vector
  pruned.withid <- unique(pruned.withid)

  # stick them into a named binary vector
  genes <- c(rep(1,length(sig.withid)),rep(0,length(pruned.withid)))

  # convert to entrez and put names on the vector
  names(genes) <- hash::values(entrez_hash$stoe,keys=c(sig.withid,pruned.withid),USE.NAMES=FALSE)

  # fit weighting function, but don't generate the plot
  pwf <- nullp(genes,genome=genome,id="knownGene",plot.fit=FALSE)

  # run goseq
  x <- goseq(pwf,genome=genome,id="knownGene",test.cats=cats)

  # replace p values with adjusted pvalues
  x$over_represented_pvalue < p.adjust(x$over_represented_pvalue,method="BH")
  x$under_represented_pvalue <- p.adjust(x$under_represented_pvalue,method="BH")
  colnames(x)[2] <- "overrep_adjP"
  colnames(x)[3] <- "underrep_adjP"

  # now pull out the stuff we want
  results <- list()
  if (enrich_type != "KEGG"){
    results$overrep <- x[x$overrep_adjP < p_cut,c("category","overrep_adjP","term","ontology")]
    results$underrep <- x[x$underrep_adjP < p_cut,c("category","underrep_adjP","term","ontology")]
  } else {
    results$overrep <- x[x$overrep_adjP < p_cut,c("category","overrep_adjP")]
    results$overrep <- cbind(results$overrep,term=hash::values(kegg_hash,keys=results$overrep$category,USE.NAMES=FALSE),stringsAsFactors=FALSE)
    results$underrep <- x[x$underrep_adjP < p_cut,c("category","underrep_adjP")]
    results$underrep <- cbind(results$underrep,term=hash::values(kegg_hash,keys=results$underrep$category,USE.NAMES=FALSE),stringsAsFactors=FALSE)
  }

  # make the go (or kegg) to gene list
  go_to_gene <- list()
  overrep_terms <- results$overrep$category
  for (i in 1:length(overrep_terms)){
    go_to_gene[[overrep_terms[i]]] <- c()
  }
  all_go_terms <- getgo(names(genes),genome=genome,id="knownGene",fetch.cats=enrich_type)
  for (g in 1:length(all_go_terms)){
    if (!is.null(all_go_terms[g][1])) {
      term_bool = overrep_terms %in% all_go_terms[g]
      for (t in 1:length(term_bool)){
        if (t == TRUE) {
          go_to_gene[[overrep_terms[t]]] <- c(go_to_gene[[overrep_terms[t]]],entrez_hash$etos)
        }
      }
    }
  }
  results$gotogenes <- go_to_gene

  return(results)
}
