library(spNNGP)

rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}


gridl <- 100
x <- seq(0, 1, length.out = 100)
coords <- as.matrix(expand.grid(x, x))
n <- dim(coords)[1]
X <- matrix(1, nrow = n, ncol = 1)
B <- matrix(4, nrow = 1, ncol = 1)
sigma.sq <- 2
tau.sq <- 1
phi <- 6
D <- as.matrix(dist(coords))
R <- exp(-phi * D)
w <- rmvn(1, rep(0, n), sigma.sq * R)
y <- rnorm(n, X %*% B + w, sqrt(tau.sq))


starting <- list("phi" = 5, "sigma.sq" = 5, "tau.sq" = 1)
tuning <- list("phi" = 0.2, "sigma.sq" = 0.2, "tau.sq" = 0.2)
priors <- list("phi.Unif" = c(3/1, 3/0.01), "sigma.sq.IG" = c(2, 5), "tau.sq.IG" = c(2, 1))
cov.model <- "exponential"
m <- spNNGP(y ~ 1, coords = coords, starting = starting, method = "response", n.neighbors = 10,
              tuning = tuning, priors = priors, cov.model = "exponential", return.neighbors = TRUE,
              n.samples = 5000, n.omp.threads = 24, n.report = 500)

round(summary(m$p.beta.samples)$quantiles, 2)
round(summary(m$p.theta.samples)$quantiles, 2)

test <- cbind(m$coords.ord, 1:n)

plot(test[, 1], test[, 2], type = "n")
text(test[, 1], test[, 2], label = as.character(test[, 3]))

plot(coords)

ypred <- spPredict(m, X, coords, n.omp.threads = 24, n.report = 500)
ypred.median <- apply(ypred$p.y.0, 1, median, na.rm = TRUE)

cols <- heat.colors(128)
image(x, x, matrix(y, ncol = gridl), useRaster = TRUE, col = cols)
image(x, x, matrix(X %*% B + w, ncol = gridl), useRaster = TRUE, col = cols)
image(x, x, matrix(ypred.median, ncol = gridl), useRaster = TRUE, col = cols)

