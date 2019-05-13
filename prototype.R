library(invgamma)
library(mvtnorm)
library(truncnorm)

phi.start <- 12
sigma.sq.start <- 1
tau.sq.start <- 0.1

phi.tune <- 0.2
sigma.sq.tune <- 0.2
tau.sq.tune <- 0.2

phi.prior <- c(3, 30)
sigma.sq.prior <- c(2, 2)
tau.sq.prior <- c(2, .1)
beta.prior <- c(0, 1000)

X <- X.comb
Xt <- t(X)
y <- y
Coords <- total.set
D.Coords <- sqrt(plgp:::distance(Coords))

Z <- diag(1, length(y))
Zt <- t(Z)
D.start <- diag(tau.sq.start, length(y))

C.start <- sigma.sq.start * plgp::covar.sep(Coords, d = rep(phi.start, 2), g = sqrt(.Machine$double.eps))
beta.start <- solve(Xt %*% X, Xt %*% y)
Vs <- solve(Zt %*% D.start %*% Z + solve(C.start))
w.start <- Vs %*% Zt %*% solve(D.start) %*% (y - X %*% beta.start)

Vb <- diag(beta.prior[2], 2)
Vbi <- solve(Vb)
mub <- rep(beta.prior[1], 2)

iters <- 1000

Univ.storage <- matrix(NA, nrow = iters, ncol = 5)
colnames(Univ.storage) <- c("phi", "sigma.sq", "tau.sq", "beta0", "beta1")
W.storage <- list()

Univ.storage[1, ] <- c(phi.start, sigma.sq.start, tau.sq.start, beta.start)
W.storage[[1]] <- w.start

for(i in 2:iters) {
  D <- diag(Univ.storage[i - 1, 3], length(y))
  Di <- solve(D)
  Cs <- Univ.storage[i - 1, 2] * exp(-Univ.storage[i - 1, 1] * D.Coords)
  Csi <- solve(Cs)
  Vb_star <- solve(Vbi + Xt %*% Di %*% X)
  mub_star <- Vbi %*% mub + Xt %*% Di %*% (y - Z %*% W.storage[[i - 1]])
  Univ.storage[i, 4:5] <- mdnngp:::rmvn(1, Vb_star %*% mub_star, Vb_star) # c(1, 5)
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

rowMeans(do.call("cbind", W.storage))
image(x, x, matrix(rowMeans(do.call("cbind", W.storage)), ncol = length(x)))
image(x, x, matrix(w, ncol = length(x)))
plot(1:iters, Univ.storage[, 1], type = "l")
plot(1:iters, Univ.storage[, 2], type = "l")
plot(1:iters, Univ.storage[, 3], type = "l")
plot(1:iters, Univ.storage[, 4], type = "l")
plot(1:iters, Univ.storage[, 5], type = "l")
