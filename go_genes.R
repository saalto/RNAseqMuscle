go_genes <- function(gene_list,double_hash,go_id,annot=org.Mm.eg.db,ontology="BP",specific=FALSE){
  require(hash)
  require(GO.db)
  # because of missing ENTREZ IDs, some requested genes in gene_list may not be matchable.
  # we therefore create a modified gene_list without these
  mod_genes <- gene_list[unname(has.key(gene_list,double_hash$stoe))]
  # use the hash to get entrez ids for the gene list
  entrez_ids <- hash::values(double_hash$stoe,keys=mod_genes,USE.NAMES=FALSE)  
  # now get a big data frame with all the GO IDs; get more general terms if specific == FALSE
  if(specific){
    all_go <- select(annot,keys=entrez_ids,columns="GO")
  }
  else {
    all_go <- select(annot,keys=entrez_ids,columns="GOALL")
  }
  # normalize the column names
  colnames(all_go) <- c("ENTREZID","GO","EVIDENCE","ONTOLOGY")
  # subset the data frame to:
  # - only select rows with the appropriate ontology
  # - ignore any terms whose evidence comes from ND, NR, or IEA
  all_go <- all_go[all_go$EVIDENCE != "ND" & all_go$EVIDENCE != "NR" & all_go$EVIDENCE != "IEA" & all_go$ONTOLOGY == ontology,]
  # filter the frame to remove anything without the go_id we want
  all_go <- all_go[all_go$GO == go_id,]
  # remove GO NAs from the result
  all_go <- all_go[!is.na(all_go$GO),]
  # now use the entrez_ids to get back the gene names
  gene_syms <- hash::values(double_hash$etos,keys=all_go$ENTREZID,USE.NAMES=FALSE)
  # remove duplicate entries
  gene_syms <- unique(gene_syms)
  # use the hash to get entrez ids for the gene list
  entrez_ids <- hash::values(double_hash$stoe,keys=gene_syms,USE.NAMES=FALSE)  

}
