#' Covariance Matrix Construction
#' 
#' Creates a sparse numeric matrix 
#' @param prim a numeric vector containing the primary ordering variable (i.e. time).
#' @param Xord an object of class \code{Xord}, returned by functions such as multistart_ordering or bin_ordering.
#' @param nn an integer specifying the maximum number of neighbors to use for each point in \code{X}.
#' @param kernf kernel function to use for kernel calculations when constructing.
#' @param ... Further arguments to pass to the distance function, such as lengthscales and a nugget.
#' @return Returns a numeric covariance matrix with sparsity determined by the number of nearest-neighbors, and non-zero elements calculated using the selected \code{distance} function.

mtx_elem_const <- function(prim, Xord, nn, kernf, d, g) {
  ord <- order(Xord$order)
  X <- Xord$X
  X.ord <- X[ord, , drop = FALSE]
  N <- dim(X.ord)[1]

  Kmat <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    Kmat[max(1, (i - nn)):i, i] <- Kmat[i, max(1, (i - nn)):i] <- 
      plgp::covar.sep(X.ord[i, , drop = FALSE], X.ord[max(1, (i - nn)):i, , drop = FALSE], d = c(1.2, 0.8), g = 0.1)
  }
  system.time({test <- chol(Kmat, sparse = TRUE)})
  system.time({test <- chol(Kmat, pivot = TRUE, tol = .Machine$double.eps)})
  system.time({Kmati <- chol2inv(chol(Kmat, pivot = TRUE, tol = .Machine$double.eps))})
  system.time(sum(diag(Kmat %*% Kmati)))
  system.time({Kmati <- solve(Kmat)})
    
}