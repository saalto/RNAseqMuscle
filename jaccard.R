#' Jaccard index for two vectors of strings
#' @param x a vector of strings (duplicates allowed)
#' @param y another vector of strings (duplicates allowed)
#' @return the Jaccard index 
#' @export
jaccard <- function(x,y) {
  I <- length(intersect(x,y))
  S <- I/(length(unique(x))+length(unique(y))-I)
  return(S)
}
