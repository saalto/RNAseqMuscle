#' Subsets an edgeR data table using a list of genes (row names), and returns the:
#'  +gene symbol with the maximum absolute logFC ($maxgene) and its sign ($maxsign)
#'  +mean and std. dev. of logFC in the group ($FCmean,$FCsd)
#'  +mean and std. dev. of absolute logFC in the group ($absFCmean,$absFCsd)
#' @param gene_list: list of gene symbols to use to subset the table
#' @param data_table: edge R table
#' @export
avg_by_group <- function(gene_list,data_table){
  # holds the results
  res <- list()
  # subset the data table
  sub_table <- data_table[row.names(data_table) %in% gene_list,]
  # find the maximum
  maxloc <- which.max(abs(sub_table$logFC))
  # this is the associated symbol and the sign of its change (1=positive, 0=zero, or -1=negative, respectively)
  res$maxgene <- row.names(sub_table)[maxloc]
  res$maxsign <- sign(sub_table$logFC[maxloc])
  # statistics for the log FC's
  res$FCmean <- mean(sub_table$logFC)
  res$FCsd <- sd(sub_table$logFC)
  # same for absolute log FC's
  res$absFCmean <- mean(abs(sub_table$logFC))
  res$absFCsd <- sd(abs(sub_table$logFC))
  return(res)
}
