#' edgeR DE gene analysis on input count file
#' @param data_file: csv files with genes names and replicate counts
#' @param groups: contrast vector for samples in the data read from data_file
#' @param adj_p: threshold for significance
#' @param drop_low: drop genes with < 10 total counts?
#' @param qlf_test: if true, assess significance with quasi-likelihood F test. Otherwise, use likelihood ratio.
#' @export
edgeR_DE <- function(data_file,groups,adj_p=0.05,drop_low=TRUE,qlf_test=TRUE){
  require("edgeR")

  # this will hold the results
  deres <- list()

  cts <- read.csv(data_file,header=FALSE,stringsAsFactors=FALSE)
  # separate data into gene names and just counts
  genes <- cts[,c(1)]
  cts <- cts[,c(2:ncol(cts))]

  # make sure groups is a factor
  groups <- factor(groups)

  # set the rownames as the gene symbols
  rownames(cts) <- genes

  # pre-filtering of stuff with very few counts
  if(drop_low){
    keep <- rowSums(cts) >= 10
    cts <- cts[keep,]
  }
  # save the gene universe
  universe <- rownames(cts)

  # now do the calculations
  y <- DGEList(counts=cts,group=groups)
  y <- calcNormFactors(y)
  design <- model.matrix(~groups)
  y <- estimateDisp(y,design)

  # here are the tests for significance: use either the quasi-likelihood F test
  #   or a likelihood ratio test
  if(qlf_test){
    fit <- glmQLFit(y,design)
    lt <- glmQLFTest(fit)
  } else {
    fit <- glmFit(y,design)
    lt <- glmLRT(fit)
  }

  # now extract the genes we want
  tagtable <- topTags(lt,n=dim(cts)[1])

  # threshold for significance
  deindex <- which(tagtable$table$FDR < adj_p, arr.ind=TRUE)

  # return the results in list
  deres$table <- tagtable$table
  deres$sig_genes <- row.names(tagtable$table)[deindex]
  deres$universe <- universe

  return(deres)
}
