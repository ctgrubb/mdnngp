#' Gaussian Process Spatial Linear Regression
#' 
#' Performs markov-chain monte carlo on the posterior of a linear regression model containing an unobserved spatial random effect that comes from a Gaussian Process.
#' @usage 
#' gpslm(formula, data = parent.frame(), coords, m = NULL, col = NULL, kernel = "exponential", priors, starting, tuning, n.samples = 1000)
#' @param formula a formula, representing the linear model to fit.
#' @param data an optional data frame containing the model variables.
#' @param coords an \emph{n x 2} numeric matrix containing the coordinates at each observation.
#' @param m an integer, specifying the number of neighbors to use for the nearest-neighbor kernel. If NULL, a dense kernel is used.
#' @param col an integer, specifying the column to order by if using a nearest-neighbor kernel. If \code{m = NULL}, this argument is ignored.
#' @param kernel a string, indicating which kernel should be used for the spatial Gaussian process. Currently, possible the only possible option is "exponential".
#' @param priors a list, containing the priors to use for the model. Should contain elements "phi", a vector of length 2; "sigma.sq", a vector of length 2; "tau.sq", a vector of length 2; "beta", a vector of length 2.
#' @param starting a list, containing starting values to use for the MCMC. Should contain elements "phi", "sigma.sq", and "tau.sq"; can also contain element "beta" if desired, otherwise an initial estimate is calculated.
#' @param tuning a list, containing the tuning parameters to use for the MCMC. Should contain elements "phi" and "sigma.sq"; full conditionals are used for all other parameters.
#' @param n.samples an integer, specifying the number of posterior samples to generate.
#' @details This function fits the  model \eqn{y = X\beta + w + \epsilon}, where \eqn{w ~ GP(0, \sigma^2 K)}, \eqn{\epsilon ~ N(0, \tau^2)}. For the exponential kernel, \eqn{K = exp(-\phi D)}, where \eqn{D} is the Euclidian distance matrix between coordinates. A flat prior is applied to \eqn{\phi}, inverse-gamma priors are applied to both \eqn{\sigma^2} and \eqn{\tau^2}, and a normal prior is applied to \eqn{\beta}.
#' @return Returns a list containing the following elements:
#'   \item{w.samples}{a matrix containing the posterior samples for the unobserved spatial process}
#'   \item{beta.samples}{a matrix containing the posterior samples for beta}
#'   \item{phi.samples}{a vector containing the posterior samples for phi}
#'   \item{tau.sq.samples}{a vector containing the posterior samples for tau.sq}
#'   \item{sigma.sq.samples}{a vector containing the posterior samples for sigma.sq}
#'   \item{accept.prob}{a vector of acceptance probabilities for parameters requiring MH}
#'   \item{time}{the total time required to perform the MCMC}
#' @export
#' @import invgamma truncnorm mvtnorm
#' @importFrom plgp distance

gpslm <- function(formula, data = parent.frame(), coords, m = NULL, col = NULL, kernel = "exponential", 
                  priors, starting, tuning, n.samples = 1000) {
  
  tic <- proc.time()[3]
  
  phi.tune <- tuning$phi
  sigma.sq.tune <- tuning$sigma.sq
  
  phi.prior <- priors$phi
  sigma.sq.prior <- priors$sigma.sq
  tau.sq.prior <- priors$tau.sq
  beta.prior <- priors$beta
  
  X <- model.matrix(formula, data)
  p <- ncol(X)
  Y <- eval(formula[[2]], envir = data)
  N <- length(Y)
  
  
  if(is.null(m)) {
    Xt <- t(X)
    D <- as.matrix(dist(coords))
  } else {
    ord <- order(coords[, col])
    undo_ord <- match(1:N, ord)
    coords <- coords[ord, ]
    X <- X[ord, ]
    Xt <- t(X)
    Y <- Y[ord]
    D <- as.matrix(dist(coords))
    n.indx <- list()
    n.indx[[1]] <- NA
    n.indx[[2]] <- 1
    for(i in 3:length(ord)) {
      n.indx[[i]] <- as.numeric(names(sort(D[i, 1:(i - 1)])))[1:min(m, i - 1)]
    }
    A <- matrix(0, nrow = N, ncol = N)
    I <- B <- diag(N)
  }
  
  w.samples <- matrix(NA, nrow = n.samples, ncol = N)
  colnames(w.samples) <- paste0("w", 1:N)
  beta.samples <- matrix(NA, nrow = n.samples, ncol = p)
  colnames(beta.samples) <- paste("beta", 0:(p - 1))
  phi.samples <- rep(NA, n.samples)
  tau.sq.samples <- rep(NA, n.samples)
  sigma.sq.samples <- rep(NA, n.samples)
  accept.prob <- rep(0, 2)
  names(accept.prob) <- c("phi", "sigma.sq")
  
  phi.samples[1] <- starting$phi
  sigma.sq.samples[1] <- starting$sigma.sq
  tau.sq.samples[1] <- starting$tau.sq
  
  denseK <- sigma.sq.samples[1] * exp(-phi.samples[1] * D)
  if(!is.null(m)) {
    B[1, 1] <- denseK[1, 1]
    for(i in 1:(N - 1)) {
      n <- n.indx[[i + 1]]
      A[i + 1, n] <- solve(denseK[n, n], denseK[n, i + 1])
      B[i + 1, i + 1] <- denseK[i + 1, i + 1] - denseK[i + 1, n] %*% A[i + 1, n]
    }
    K <- solve(I - A) %*% B %*% t(solve(I - A))
    Ki <- t(I -A) %*% solve(B) %*% (I - A)
  } else {
    K <- denseK
    Ki <- solve(K)
  }
  
  Vi <- diag(1 / tau.sq.samples[1], N)
  beta.samples[1, ] <- solve(Xt %*% X, Xt %*% Y)
  
  Vs <- solve(Vi + Ki)
  w.samples[1, ] <- Vs %*% Vi %*% (Y - X %*% beta.samples[1, ])
  
  Vb <- diag(beta.prior[2], p)
  Vbi <- solve(Vb)
  mub <- rep(beta.prior[1], p)
  
  for(i in 2:n.samples) {
    V <- diag(tau.sq.samples[i - 1], N)
    Vi <- solve(V)
    denseK <- sigma.sq.samples[i - 1] * exp(-phi.samples[i - 1] * D)
    if(!is.null(m)) {
      B[1, 1] <- denseK[1, 1]
      for(j in 1:(N - 1)) {
        n <- n.indx[[j + 1]]
        A[j + 1, n] <- solve(denseK[n, n], denseK[n, j + 1])
        B[j + 1, j + 1] <- denseK[j + 1, j + 1] - denseK[j + 1, n] %*% A[j + 1, n]
      }
      K <- solve(I - A) %*% B %*% t(solve(I - A))
      Ki <- t(I -A) %*% solve(B) %*% (I - A)
    } else {
      K <- denseK
      Ki <- solve(K)
    }
    Vb_star <- solve(Vbi + Xt %*% Vi %*% X)
    mub_star <- Vbi %*% mub + Xt %*% Vi %*% (Y - w.samples[i - 1, ])
    beta.samples[i, ] <- mdnngp:::rmvn(1, Vb_star %*% mub_star, Vb_star + sqrt(.Machine$double.eps))
    Vs <- solve(Vi + Ki)
    w.samples[i, ] <- mdnngp:::rmvn(n = 1, mu = Vs %*% Vi %*% (Y - X %*% beta.samples[i, ]), V = Vs + sqrt(.Machine$double.eps))
    error <- Y - X %*% beta.samples[i, ] - w.samples[i, ]
    tau.sq.samples[i] <- rinvgamma(1, shape = N / 2 + tau.sq.prior[1], rate = tau.sq.prior[2] + 0.5 * t(error) %*% error)
    
    # MH proposals
    phi.proposal <- rtruncnorm(1, phi.prior[1], phi.prior[2], phi.samples[i - 1], phi.tune)
    sigma.sq.proposal <- rtruncnorm(1, 0, Inf, sigma.sq.samples[i - 1], sigma.sq.tune)
    
    # MH acceptance - phi
    target_new <- as.numeric(dmvnorm(as.numeric(w.samples[i, ]), sigma = sigma.sq.samples[i - 1] * 
                                       exp(-phi.proposal * D), log = TRUE) + 
                               dunif(phi.proposal, phi.prior[1], phi.prior[2], log = TRUE))
    
    proposal_new <- log(dtruncnorm(phi.proposal, phi.prior[1], phi.prior[2], phi.samples[i - 1], phi.tune))
    
    target_old <- as.numeric(dmvnorm(as.numeric(w.samples[i, ]), sigma = sigma.sq.samples[i - 1] * 
                                       exp(-phi.samples[i - 1] * D), log = TRUE) + 
                               dunif(phi.samples[i - 1], phi.prior[1], phi.prior[2], log = TRUE))
    
    proposal_old <- log(dtruncnorm(phi.samples[i - 1], phi.prior[1], phi.prior[2], phi.proposal, phi.tune))
    
    if(log(runif(1)) < target_new + proposal_old - target_old - proposal_new) {
      accept.prob[1] <- accept.prob[1] + 1
      phi.samples[i] <- phi.proposal
    } else {
      phi.samples[i] <- phi.samples[i - 1]
    }
    
    # MH acceptance - sigma.sq
    target_new <- as.numeric(dmvnorm(as.numeric(w.samples[i, ]), sigma = sigma.sq.proposal * 
                                       exp(-phi.samples[i] * D), log = TRUE) + 
                               dinvgamma(sigma.sq.proposal, shape = sigma.sq.prior[1], rate = sigma.sq.prior[2], log = TRUE))
    
    proposal_new <- log(dtruncnorm(sigma.sq.proposal, 0, Inf, sigma.sq.samples[i - 1], sigma.sq.tune))
    
    target_old <- as.numeric(dmvnorm(as.numeric(w.samples[i, ]), sigma = sigma.sq.samples[i - 1] * 
                                       exp(-phi.samples[i] * D), log = TRUE) + 
                               dinvgamma(sigma.sq.samples[i - 1], shape = sigma.sq.prior[1], rate = sigma.sq.prior[2], log = TRUE))
    
    proposal_old <- log(dtruncnorm(sigma.sq.samples[i - 1], 0, Inf, sigma.sq.proposal, sigma.sq.tune))
    
    if(log(runif(1)) < target_new + proposal_old - target_old - proposal_new) {
      accept.prob[2] <- accept.prob[2] + 1
      sigma.sq.samples[i] <- sigma.sq.proposal
    } else {
      sigma.sq.samples[i] <- sigma.sq.samples[i - 1]  
    }
    
  }
  
  if(!is.null(m)) {
    w.samples <- w.samples[, undo_ord]
  }
  
  toc <- proc.time()[3]
  out <- list(w.samples = w.samples, beta.samples = beta.samples, phi.samples = phi.samples, 
              tau.sq.samples = tau.sq.samples, sigma.sq.samples = sigma.sq.samples, 
              accept.prob = accept.prob, time = toc - tic)
  return(out)
  
}