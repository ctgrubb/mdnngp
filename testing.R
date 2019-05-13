rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}


gridl <- 50
x <- seq(0, 1, length.out = gridl)
coords <- as.matrix(expand.grid(x, x))
n <- dim(coords)[1]
X <- cbind(1, rnorm(n))
B <- as.matrix(c(1, 5))
sigma.sq <- 1
tau.sq <- 0.1
phi <- 12
D <- as.matrix(dist(coords))
R <- exp(-phi * D)
w <- rmvn(1, rep(0, n), sigma.sq * R)
y <- rnorm(n, X %*% B + w, sqrt(tau.sq))

cols <- heat.colors(128)
image(x, x, matrix(w, ncol = gridl), useRaster = TRUE, col = cols)

library(spNNGP)

n.samples <- 500

starting <- list("phi"=5, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

m.s <- spNNGP(y~X-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2)

round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)
plot(apply(m.s$p.w.samples, 1, median), w)

m.r <- spNNGP(y~X-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2)

round(summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)

