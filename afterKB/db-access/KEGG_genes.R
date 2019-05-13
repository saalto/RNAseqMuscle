#' Returns a list of genes, taken from an input list, that have a given KEGG pathway number.  Also makes
#'  chord data matrix (from chord_dat) for use with GOplot
#' @param gene_list: list of gene symbols
#' @param double_hash: R hash() giving symbol-to-entrez ($stoe) mapping AND entrez-to-symbol ($etos) mapping
#' @param kegg_id: kegg_id to search for
#' @param annot: AnnotationDBI object to use
#' @param make_chord: makes data structure for
#' @export
KEGG_genes <- function(gene_list,double_hash,kegg_id,annot=org.Mm.eg.db){
  require(hash)
  require(GO.db)
  # because of missing ENTREZ IDs, some requested genes in gene_list may not be matchable.
  # we therefore create a modified gene_list without these
  mod_genes <- gene_list[unname(has.key(gene_list,double_hash$stoe))]
  # use the hash to get entrez ids for the gene list
  entrez_ids <- hash::values(double_hash$stoe,keys=mod_genes,USE.NAMES=FALSE)
  # now get a big data frame with all the KEGG IDs
  all_kegg <- select(annot,keys=entrez_ids,columns="PATH")
  # remove PATH NAs from the result
  all_kegg <- all_kegg[!is.na(all_kegg$PATH),]
  # filter the frame to remove anything without the kegg_id we want
  all_kegg <- all_kegg[all_kegg$PATH == kegg_id,]
  # now use the entrez_ids to get back the gene names
  gene_syms <- hash::values(double_hash$etos,keys=all_kegg$ENTREZID,USE.NAMES=FALSE)
}

#' script for assembling KEGG data for use with GOplot. returns a data frame of the circ
#'  type (see GOplot docs) and a matrix for use with GOplot (see chord_dat docs)
#' @param gene_frame: data frame with fold change/p value info, as well as gene names
#' @param double_hash: symbol-to-entrez and entrez-to-symbol hash map
#' @param kegg_hash: KEGG id to pathway name/descriptor
#' @param kegg_ids: vector of IDs to check against
#' @param annot: organism database
#' @export
KEGG_assemble <- function(gene_frame,double_hash,kegg_hash,kegg_ids,annot=org.Mm.eg.db){
  require(GOplot)
  # these will hold the results as we loop through IDs
  ID = c()
  term = c()
  count = c()
  genes = c()
  logFC = c()
  # go through KEGG ids one by one
  for (i in 1:length(kegg_ids)){
    k_id = kegg_ids[i]
    # pull the gene list for the kegg ID
    g <- KEGG_genes(row.names(gene_frame),double_hash,k_id,annot)
    genes <- c(genes,g)
    # now info
    count <- c(count,rep(length(genes),length(g)))
    ID <- c(ID,rep(k_id,length(g)))
    term <- c(term,rep(hash::values(kegg_hash,keys=k_id,USE.NAMES=FALSE),length(g)))
    # pull logFC
    logFC <- c(logFC,gene_frame[row.names(gene_frame) %in% g,]$logFC)
  }
  # category is a dummy if we aren't dealing with GO ids
  category <- rep('NONE',length(genes))
  # this mimics the data frame produced by circ_dat
  circ <- data.frame(category,ID,term,count,genes,logFC,stringsAsFactors=FALSE)
  # produce the one that chord data make; need the "genes" data frame
  u_genes <- unique(circ$genes)
  u_logFC <- gene_frame[row.names(gene_frame) %in% u_genes,]$logFC
  genes_for_chord = data.frame("ID"=u_genes,"logFC"=u_logFC,stringsAsFactors=FALSE)
  # this is the one that chord_dat makes
  chord <- chord_dat(circ,genes_for_chord,unique(circ$term))
  # put them both in a list
  goplot_data <- list("circ" = circ, "chord" = chord)
}
