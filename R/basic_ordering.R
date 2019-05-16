#' Nearest-Neighbor Basic Ordering
#' 
#' Create gridded observation ordering based on ordering in one variable, and randomizing the rest.
#' @usage 
#' basic_ordering(X, col)
#' @param X a numeric matrix, containing numeric, or integer values; should not include primary ordering variable (i.e. time). These points should be replicated at each value of the primary ordering variable. Each column corresponds to a dimension of the space.
#' @param col an integer, specifying which column of X to order by (all other column orders are randomized).
#' @return Returns an object of class \code{Xord}, which is a list containing the following elements:
#'   \item{X}{A matrix containing the original X matrix}
#'   \item{order}{A vector giving the ordering of observations}
#' @export

basic_ordering <- function(X, col = 1) {
  Xdf <- data.frame(X, label = 1:dim(X)[1])
  Xdf <- Xdf[sample(1:dim(X)[1], dim(X)[1]), ]
  ord <- order(Xdf[, col])
  Xdf <- Xdf[ord, ]
  out <- list(X = X, order = as.numeric(Xdf$label), bins = NULL, score = NULL)
  class(out) <- "Xord"
  return(out)
}