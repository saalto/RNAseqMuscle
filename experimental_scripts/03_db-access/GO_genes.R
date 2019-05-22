#' Returns a list of genes, taken from an input list, that have a given GO term ID
#' @param gene_list: list of gene symbols
#' @param double_hash: R hash() giving symbol-to-entrez ($stoe) mapping AND entrez-to-symbol ($etos) mapping
#' @param go_id: go_id to search for
#' @param annot: AnnotationDBI object to use
#' @export

GO_genes <- function(gene_frame,double_hash,go_ids,annot=org.Mm.eg.db){
  require(hash)
  require(org.Mm.eg.db)
  # because of missing ENTREZ IDs, some requested genes in gene_list may not be matchable.
  # we therefore create a modified gene_list without these
  mod_genes <- gene_frame[unname(has.key(row.names(gene_frame),double_hash$stoe)),]
  # use the hash to get entrez ids for the gene list
  entrez_ids <- hash::values(double_hash$stoe,keys=row.names(mod_genes),USE.NAMES=FALSE)
  # now get a big data frame with all the GO IDs
  all_go <- AnnotationDbi::select(annot,keys=entrez_ids,columns="GOALL", keytype="ENTREZID")
  # remove NAs from the result
  all_go <- all_go[!is.na(all_go$GO),]
  # filter the frame to remove anything other than BP we want
  all_go <- all_go[all_go$ONTOLOGY == "BP",]
  #filter the frame to remove anything other than the GO term we want
  all_go <- all_go[all_go$GO %in% go_ids,]
  # now use the entrez_ids to get back the gene names
  gene_syms <- hash::values(double_hash$etos,keys=all_go$ENTREZID,USE.NAMES=FALSE)
  all_go$GENES <- gene_syms
  return(all_go)
}


#' script for assembling GO data for use with GOplot. returns a data frame of the circ
#'  type (see GOplot docs) and a matrix for use with GOplot (see chord_dat docs)
#' @param gene_frame: data frame with fold change/p value info, as well as gene names
#' @param double_hash: symbol-to-entrez and entrez-to-symbol hash map
#' @param overrep_file: Over-represented GO term list generated from catenrich.R function
#' @param go_ids: vector of IDs to check against
#' @param annot: organism database
#' @export
GO_assemble <- function(gene_frame,double_hash,overrep_file,go_ids, annot=org.Mm.eg.db){
  require(GOplot)
  # these will hold the results as we loop through IDs
  ID = c()
  term = c()
  count = c()
  genes = c()
  logFC = c()
  all_go <- GO_genes(gene_frame,double_hash,go_ids,annot=org.Mm.eg.db)
  
  # go through GO ids one-by-one and pull gene information
  for (i in 1:length(go_ids)){
    go_id = go_ids[i]
    # pull the gene list for specific GO ID
    g <- all_go$GENES[all_go$GO %in% go_id]
    genes <- c(genes,g)
    # now info
    count <- c(count,rep(length(genes),length(g)))
    ID <- c(ID,rep(go_id,length(g)))
  }
  
  # category is BP when dealing with GO ids
  category <- rep('BP',length(genes))
  
  # pull logFC and term
  for(i in 1:length(genes)){
    g = genes[i]
    ident = ID[i]
    logFC <- c(logFC, gene_frame[row.names(gene_frame) %in% g,]$logFC)
    overrep_term <- overrep_file$overrep[overrep_file$overrep$category %in% ident,]
    term <- c(term, overrep_term$term)
  }
  
  # this mimics the data frame produced by circ_dat
  circ <- data.frame(category,ID,count,term,genes,logFC,stringsAsFactors=FALSE)
  # produce the one that chord data make; need the "genes" data frame
  u_genes <- unique(circ$genes)
  u_logFC <- gene_frame[row.names(gene_frame) %in% u_genes,]$logFC
  genes_for_chord = data.frame("ID"=u_genes,"logFC"=u_logFC,stringsAsFactors=FALSE)
  #genes_for_chord = genes_for_chord[genes_for_chord$ID %in% genelist$GENE,]
  # this is the one that chord_dat makes
  chord <- chord_dat(circ,genes_for_chord,unique(circ$term))
  # put them both in a list
  goplot_data <- list("circ" = circ, "chord" = chord)
}
