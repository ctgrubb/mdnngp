#' Maximin Design
#' 
#' Create a stochastically optimized maximin design.
#' @usage 
#' maximin(n = 100, m = 2, T = 10000, X.init = NULL)
#' @param n an integer, specifying the number of samples to generate.
#' @param m an integer, specifying the number of dimensions.
#' @param T an integer, specifying the number of stochastic optimization iterations to perform.
#' @param X.init a matrix, representing an initial design. If not null, values of n and m are obtained from this initial design.
#' @return Returns a list containing the following elements:
#'   \item{X}{A matrix, representing the optimized design}
#'   \item{prog}{A vector, containing the value of the criterion (minimum distance between points) over the \code{T} iterations}
#'   \item{time}{The amount of time taken to complete the optimization (in seconds)}
#' @export

maximin <- function(n = 100, m = 2, T = 10000, X.init = NULL) {
  tic <- proc.time()[3]
  if(!is.null(X.init)) {
    X <- X.init
    n <- nrow(X)
    m <- ncol(X)
  } else {
    X <- matrix(runif(n * m), ncol = m)
  }
  d <- distance(X)
  d.low <- d[lower.tri(d)]
  wmd <- which.min(d.low)
  md <- d.low[wmd]
  prog <- rep(0, T)
  for(t in 1:T) {
    row <- sample(ltoi(wmd, n), 1)
    xold <- X[row, , drop = FALSE]
    X[row, ] <- runif(m)
    d.vec <- distance(X, X[row, , drop = FALSE])
    mdprime <- min(d.vec[-row])
    if(mdprime > md) {
      d[row, ] <- d[, row] <- d.vec
      d.low <- d[lower.tri(d)]
      wmd <- which.min(d.low)
      md <- d.low[wmd]
    } else {
      X[row, ] <- xold
    }
    prog[t] <- md
  }
  toc <- proc.time()[3]
  return(list(X = X, prog = prog, time = toc - tic))
}