#' DESeq2 DE gene analysis on input count file
#' Documentation coming soon!
#' @param data_file
#' @param col_names
#' @param condition
#' @param adj_p
#' @param drop_low
#' @export
DESeq2_DE <- function(data_file,col_names,condition,adj_p=0.05,drop_low=TRUE){
  require("DESeq2")
  # will hold the results
  deres <- list()

  # read in count data
  cts <- read.csv(data_file,header=FALSE)

  # separate into gene names and just counts
  genes <- cts[,c(1)]
  cts <- cts[,c(2:12)]

  # change row/col names (function will fail here if Excel was used to open and save the data file)
  rownames(cts) <- genes
  colnames(cts) <- col_names

  coldata <- data.frame(colnames(cts), "condition" = condition)

  dds <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design = ~ condition)

  # pre-filtering of stuff with very few counts
  if(drop_low) {
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
  }

  # run the analysis, with the requested p-value cutoff
  dds <- DESeq(dds)
  res <- results(dds,alpha=adj_p)

  # find significant genes
  deindex <- which(res$padj < adj_p, arr.ind=TRUE)
  deres$sig_genes <- row.names(res)[deindex]
  deres$universe <- row.names(res)
  deres$table <- res
  return(deres)
}
