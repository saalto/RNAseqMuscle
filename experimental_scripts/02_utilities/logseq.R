#' Function for generating a sequence of logarithmically spaced points
#' @param from starting point for sequence
#' @param to end point for sequence
#' @param length.out number of points in the interval
#' @return a vector of logarithmically spaced points
#' @export
logseq <- function(from,to,length.out) {
  s <- exp(seq(log(from),log(to),length.out = length.out))
  return(s)
}
