#' KEGG pathway counter for input list of genes
#' KEGG terms are identified using the genes and organized based on the count and p-value.
#' @param gene_list: list of gene symbols
#' @param sym_hash: R hash() giving symbol-to-entrez mapping
#' @param kegg_hash: R hash() with ID-to-descriptor mapping for KEGG
#' @param annot: AnnotationDBI object to use
#' @export
KEGG_counter <- function(gene_list,sym_hash,kegg_hash,annot=org.Mm.eg.db){
  require(hash)
  require(GO.db)
  # because of missing ENTREZ IDs, some requested genes in gene_list may not be matchable.
  # we therefore create a modified gene_list without these
  mod_genes <- gene_list[unname(has.key(gene_list,sym_hash))]
  # use the hash to get entrez ids for the gene list
  entrez_ids <- hash::values(sym_hash,keys=mod_genes,USE.NAMES=FALSE)
  # now get a big data frame with all the KEGG IDs
  all_kegg <- select(annot,keys=entrez_ids,columns="PATH")
  # counts of all terms which appear
  all_kegg <- as.data.frame(table(all_kegg$PATH),stringsAsFactors=FALSE)
  # obtain the list of pathway descriptors from the hash and bind them to the data frame
  kegg_desc <- hash::values(kegg_hash,keys=all_kegg$Var1,USE.NAMES=FALSE)
  all_kegg <- cbind(all_kegg,kegg_desc)
  # remove row numbers
  kegg_ids <- all_kegg$Var1
  all_kegg <- all_kegg[,2:3]
  rownames(all_kegg) <- kegg_ids
  # better column names
  colnames(all_kegg) <- c("count","descriptor")
  # sort in decreasing count order
  all_kegg <- all_kegg[order(all_kegg$count,decreasing=TRUE),]
}
