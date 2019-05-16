#' Nearest-Neighbor Distance Matrix Construction
#' 
#' Creates a sparse numeric matrix containing pairwise distances between neighbors.
#' @param Xord an object of class \code{Xord}, returned by functions such as multistart_ordering or bin_ordering.
#' @param m an integer specifying the maximum number of neighbors to use for each point in \code{X}.
#' @return Returns a numeric distance matrix with sparsity determined by the number of nearest-neighbors.
#' @importFrom plgp distance

mtx_D_const <- function(Xord, m) {
  ord <- Xord$order
  X <- Xord$X
  X.ord <- X[ord, , drop = FALSE]
  N <- dim(X.ord)[1]

  Kmat <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    Kmat[max(1, (i - m)):i, i] <- Kmat[i, max(1, (i - m)):i] <- 
      sqrt(plgp::distance(X.ord[i, , drop = FALSE], X.ord[max(1, (i - m)):i, , drop = FALSE]))
  }
  return(Kmat)
}