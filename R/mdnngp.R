#' Gaussian Process Spatial Linear Regression
#' 
#' Performs markov-chain monte carlo on the posterior of a linear regression model containing an unobserved spatial random effect that comes from a Gaussian Process.
#' @usage 
#' gpslm(formula, data = parent.frame(), coords, m = NULL, kernel = "exponential", priors, starting, tuning, n.samples = 1000)
#' @param formula a formula, representing the linear model to fit.
#' @param data an optional data frame containing the model variables.
#' @param coords an \emph{n x 2} numeric matrix containing the coordinates at each observation.
#' @param m an integer, specifying the number of neighbors to use for the nearest-neighbor kernel. If NULL, a dense kernel is used.
#' @param kernel a string, indicating which kernel should be used for the spatial Gaussian process. Currently, possible the only possible option is "exponential".
#' @param priors a list, containing the priors to use for the model. Should contain elements "phi.unif", a vector of length 2; "sigma.sq.ig", a vector of length 2; "tau.sq.ig", a vector of length 2; "beta.norm", a vector of length 2.
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
#' @import invgamma truncnorm
#' @importFrom plgp distance

gpslm <- function(formula, data = parent.frame(), coords, m = NULL, kernel = "exponential", priors, starting, tuning, n.samples = 1000) {
  
  tic <- proc.time()
  
  phi.tune <- tuning$phi
  sigma.sq.tune <- tuning$sigma.sq
  
  phi.prior <- priors$phi
  sigma.sq.prior <- priors$sigma.sq
  tau.sq.prior <- priors$tau.sq
  beta.prior <- priors$beta
  
  X <- model.matrix(formula, data)
  p <- ncol(X)
  Xt <- t(X)
  Y <- eval(formula[[2]], envir = data)
  N <- length(Y)
  D <- sqrt(plgp:::distance(coords))
  
  w.samples <- matrix(NA, nrow = n.samples, ncol = N)
  colnames(w.samples) <- paste0("w", 1:N)
  beta.samples <- matrix(NA, nrow = n.samples, ncol = p)
  colnames(beta.samples) <- paste("beta", 0:(p - 1))
  phi.samples <- rep(NA, n.samples)
  tau.sq.samples <- rep(NA, n.samples)
  sigma.sq.samples <- rep(NA, n.samples)
  accept.prob <- rep(NA, 2)
  names(accept.prob) <- c("phi", "sigma.sq")
  
  phi.samples[1] <- starting$phi
  sigma.sq.samples[1] <- starting$sigma.sq
  tau.sq.samples[1] <- starting$tau.sq
  
  K <- sigma.sq.samples[1] * exp(-phi.samples[1] * D)
  Ki <- solve(K)
  Vi <- diag(1 / tau.sq.samples, N)
  beta.samples[1, ] <- solve(Xt %*% X, Xt %*% Y)
  
  Vs <- solve(Vi + Ki)
  w.samples[1, ] <- Vs %*% Vi %*% (Y - X %*% beta.samples[1, ])
  
  Vb <- diag(beta.prior[2], p)
  Vbi <- solve(Vb)
  mub <- rep(beta.prior[1], p)
  
  for(i in 2:n.samples) {
    V <- diag(tau.sq.samples[i - 1], N)
    Vi <- solve(V)
    K <- sigma.sq.samples[i - 1] * exp(-phi.samples[i - 1] * D)
    Ki <- solve(K)
    Vb_star <- solve(Vbi + Xt %*% Vi %*% X)
    mub_star <- Vbi %*% mub + Xt %*% Vi %*% (Y - w.samples[i - 1, ])
    beta.samples[i, ] <- mdnngp:::rmvn(1, Vb_star %*% mub_star, Vb_star) # c(1, 5)
    Vs <- solve(Zt %*% Di %*% Z + Csi)
    W.storage[[i]] <- mdnngp:::rmvn(1, Vs %*% Zt %*% Di %*% (y - X %*% Univ.storage[i, 4:5]), Vs)
    error <- y - X %*% Univ.storage[i, 4:5] - Z %*% W.storage[[i]]
    Univ.storage[i, 3] <- rinvgamma(1, shape = length(y) / 2 + tau.sq.prior[1], rate = tau.sq.prior[2] + 0.5 * t(error) %*% error) # 0.1
    Univ.storage[i, 1] <- rtruncnorm(1, phi.prior[1], phi.prior[2], Univ.storage[i - 1, 1], phi.tune)
    target_new <- as.numeric(dmvnorm(as.numeric(W.storage[[i]]), sigma = Univ.storage[i - 1, 2] * 
                                       exp(-Univ.storage[i, 1] * D.Coords), log = TRUE) + 
                               dunif(Univ.storage[i, 1], phi.prior[1], phi.prior[2], log = TRUE))
    proposal_new <- log(dtruncnorm(Univ.storage[i, 1], phi.prior[1], phi.prior[2], Univ.storage[i - 1, 1], phi.tune))
    target_old <- as.numeric(dmvnorm(as.numeric(W.storage[[i]]), sigma = Univ.storage[i - 1, 2] * 
                                       exp(-Univ.storage[i - 1, 1] * D.Coords), log = TRUE) + 
                               dunif(Univ.storage[i - 1, 1], phi.prior[1], phi.prior[2], log = TRUE))
    proposal_old <- log(dtruncnorm(Univ.storage[i - 1, 1], phi.prior[1], phi.prior[2], Univ.storage[i, 1], phi.tune))
    if(log(runif(1)) > target_new + proposal_old - target_old - proposal_new) {
      Univ.storage[i, 1] <- Univ.storage[i - 1, 1]
    }
    Univ.storage[i, 2] <- rtruncnorm(1, 0, Inf, Univ.storage[i - 1, 2], sigma.sq.tune)
    target_new <- as.numeric(dmvnorm(as.numeric(W.storage[[i]]), sigma = Univ.storage[i, 2] * 
                                       exp(-Univ.storage[i, 1] * D.Coords), log = TRUE) + 
                               dinvgamma(Univ.storage[i, 2], shape = sigma.sq.prior[1], rate = sigma.sq.prior[2], log = TRUE))
    proposal_new <- log(dtruncnorm(Univ.storage[i, 2], 0, Inf, Univ.storage[i - 1, 2], sigma.sq.tune))
    target_old <- as.numeric(dmvnorm(as.numeric(W.storage[[i]]), sigma = Univ.storage[i - 1, 2] * 
                                       exp(-Univ.storage[i, 1] * D.Coords), log = TRUE) + 
                               dinvgamma(Univ.storage[i - 1, 2], shape = sigma.sq.prior[1], rate = sigma.sq.prior[2], log = TRUE))
    proposal_old <- log(dtruncnorm(Univ.storage[i - 1, 2], 0, Inf, Univ.storage[i, 2], sigma.sq.tune))
    if(log(runif(1)) > target_new + proposal_old - target_old - proposal_new) {
      Univ.storage[i, 2] <- Univ.storage[i - 1, 2]
    }
  }
  
  
}