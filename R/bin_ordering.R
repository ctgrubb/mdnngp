#' Nearest-Neighbor Binned Ordering
#'
#' Create gridded observation ordering based on multi-dimensional binning in non-primary ordering variable.
#' @usage
#' bin_ordering(X, bins, iters = 1000)
#' @param X a numeric matrix, containing numeric, or integer values; should not include primary ordering variable (i.e. time). These points should be replicated at each value of the primary ordering variable. Each column corresponds to a dimension of the space.
#' @param bins an integer specifying the number of bins per dimension of X, or a numeric vector containing the number of bins to use for each dimension of X. If vector, length(bins) should be equal to dim(X)[2]. Number of bins in any dimension must be less than or equal to the number of unique values in that dimension.
#' @param iters Number of iterations of stochastic swapping to perform.
#' @return Returns an object of class \code{Xord}, which is a list containing the following elements:
#'   \item{X}{A matrix containing the original X matrix}
#'   \item{order}{A vector giving the ordering of observations}
#'   \item{bins}{A matrix containing the bin assignments for each observation}
#'   \item{score}{Pooled variance of orderings across bins; used to judge which ordering is best when conducted multiple times}
#' @export
#' @importFrom plgp distance

bin_ordering <- function(X, bins, iters = 1000) {
  if(!(length(bins) == 1 || length(bins) == dim(X)[2])) {
    stop("Length of bins does not match size of X.")
  }
  bins <- rep(bins, length.out = dim(X)[2])
  n_uniqs <- apply(X, 2, function(x) {length(unique(x))})
  if(any(n_uniqs < bins)) {
    stop("Number of unique values in at least 1 dimensions is less than the requested number of bins.")
  }
  Xdf <- data.frame(X, label = sample(1:dim(X)[1], dim(X)[1]))
  for(i in 1:dim(X)[2]) {
      Xdf[[paste0("bin", i)]] <- as.numeric(cut(X[, i], breaks = bins[i], labels = 1:bins[i]))
  }
  ts <- matrix(NA, nrow = iters, ncol = dim(X)[2])
  ts[1, ] <- apply(Xdf[, (dim(X)[2] + 2):(2 * dim(X)[2] + 1)], 2, 
                   function(x) {kruskal.test(Xdf$label, x)$statistic})
  for(i in 2:iters) {
    swap <- sample(1:dim(Xdf)[1], 1)
    temp <- Xdf$label[swap]
    swap_to <- sample((1:dim(Xdf)[1])[-swap], 1)
    Xdf$label[swap] <- Xdf$label[swap_to]
    Xdf$label[swap_to] <- temp
    ts[i, ] <- apply(Xdf[, (dim(X)[2] + 2):(2 * dim(X)[2] + 1)], 2, 
                     function(x) {kruskal.test(Xdf$label, x)$statistic})
    if(!all(ts[i, ] < ts[i - 1, ])) {
      Xdf$label[swap_to] <- Xdf$label[swap]
      Xdf$label[swap] <- temp
      ts[i, ] <- ts[i - 1, ]
    }
  }
  meanlist <- list()
  meanmatrix <- matrix(NA, nrow = dim(X)[2], ncol = 2)
  for(i in 1:dim(X)[2]) {
    meanlist[[paste0(i)]] <- as.numeric(aggregate(Xdf$label, by = list(bins = Xdf[[paste0("bin", i)]]), FUN = mean)[, 2])
    meanmatrix[i, ] <- c(var(meanlist[[paste0(i)]]), length(meanlist[[paste0(i)]]) - 1)
  }
  pv <- sum(meanmatrix[, 1] * meanmatrix[, 2]) / sum(meanmatrix[, 2])
  out <- list(X = as.matrix(Xdf[, 1:dim(X)[2]]), 
              order = as.numeric(Xdf[, dim(X)[2] + 1]), 
              bins = as.matrix(Xdf[, (dim(X)[2] + 2):(2 * dim(X)[2] + 1)]),
              score = pv)
  class(out) <- "Xord"
  return(out)
}