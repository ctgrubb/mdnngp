#' Parallel Implementations of Ordering using a Multistart Scheme
#' 
#' Performs a user-defined number of parallel implementations of an ordering scheme, and then selects the best design under the pooled-variance criterion
#' @usage
#' multistart_ordering(X, type = c("random", "binned"), multistart = 10, ...)
#' @param X a matrix, containing numeric, or integer values; should not include primary ordering variable (i.e. time). These points should be replicated at each value of the primary ordering variable. Each column corresponds to a dimension of the space.
#' @param type a string corresponding to which type of ordering to perform. Options inclode "random", which randomly orders each observation, and "binned", which attempts to uniformly order across dimensions using bins.
#' @param multistart an integer specifying the number of independent, parallelized iterations of ordering.
#' @param ... optional arguments regarding parallelization.
#' @return Returns an object of class \code{Xord}, which is a list containing the following elements:
#'   \item{X}{A matrix containing the original X matrix}
#'   \item{order}{A vector giving the ordering of observations}
#'   \item{bins}{A matrix containing the bin assignments for each observation}
#'   \item{score}{Pooled variance of orderings across bins; used to judge which ordering is best when conducted multiple times}
#' @export
#' @import parallel

multistart_ordering <- function(X, type = "binned", multistart = 10, bins = NULL, iters = NULL, ...) {
  if(!(type == "binned" || type == "random")) stop("Type not recognized")
  if(!hasArg("par")) {
    par <- "PSOCK"
    if(Sys.info()[1] == "Linux") par <- "FORK"
  }
  if(!hasArg("ncores")) ncores <- min(max(1, parallel::detectCores() / 2), multistart)
  if(type == "binned") {
    if(is.null(bins)) stop("Binned method requires a bins parameter")
    if(is.null(iters)) stop("Binned method requires a number of iteration")
  }
  cl <- parallel::makeCluster(ncores, type = par, outfile = "")
  on.exit(stopCluster(cl), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  clusterSetRNGStream(cl, .Random.seed)
  clusterExport(cl, list("bins", "iters"), envir = environment())
  out <- clusterEvalQ(cl, bin_ordering(X, bins, iters))
  scores <- rep(NA, multistart)
  for(i in 1:multistart) scores[i] <- out[[i]]$score
  wm <- which.min(scores)
  out <- out[[wm]]
  class(out) <- "Xord"
  return(out)
}