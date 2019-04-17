#' GO term counter for input list of genes
#' Documentation coming soon!
#' @param gene_list: list of gene symbols
#' @param sym_hash: R hash() giving symbol-to-entrez mapping
#' @param annot: AnnotationDBI object to use
#' @param ontology: ontology term to fetch (BP, MF, or CC)
#' @param specific: ignore general terms?
#' @export
GO_counter <- function(gene_list,sym_hash,annot=org.Mm.eg.db,ontology="BP",specific=TRUE){
  require(hash)
  require(GO.db)
  # because of missing ENTREZ IDs, some requested genes in gene_list may not be matchable.
  # we therefore create a modified gene_list without these
  mod_genes <- gene_list[unname(has.key(gene_list,sym_hash))]
  # use the hash to get entrez ids for the gene list
  entrez_ids <- hash::values(sym_hash,keys=mod_genes,USE.NAMES=FALSE)
  # now get a big data frame with all the GO IDs; get more general terms if specific == FALSE
  if(specific){
    all_go <- select(annot,keys=entrez_ids,columns="GO")
  }
  else {
    all_go <- select(annot,keys=entrez_ids,columns= "GOALL")
  }
  # normalize the column names
  colnames(all_go) <- c("ENTREZID","GO","EVIDENCE","ONTOLOGY")
  # subset the data frame to:
  # - only select rows with the appropriate ontology
  # - ignore any terms whose evidence comes from ND, NR, or IEA
  all_go <- all_go[all_go$EVIDENCE != "ND" & all_go$EVIDENCE != "NR" & all_go$EVIDENCE != "IEA" & all_go$ONTOLOGY == ontology,]
  # counts of all terms which appear
  all_go <- as.data.frame(table(all_go$GO),stringsAsFactors=FALSE)
  # fetch descriptors and bind to the data frame
  all_go <- cbind(all_go,unname(Term(all_go$Var1)))
  # remove row numbers
  go_terms <- all_go$Var1
  all_go <- all_go[,2:3]
  rownames(all_go) <- go_terms
  # better column names
  colnames(all_go) <- c("count","descriptor")
  # sort in decreasing count order
  all_go <- all_go[order(all_go$count,decreasing=TRUE),]
}
