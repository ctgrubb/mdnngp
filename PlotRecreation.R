library(raster)
library(spBayes)
library(ggplot2)
library(scales)
library(mdnngp)

N <- 2500 # Size of synthetic data in each dimension, i.e. total size gridl^d
train.prop <- 0.8 # Proportion for training

total.set <- matrix(runif(N * 2), ncol = 2)

train.size <- round(train.prop * dim(total.set)[1])
test.size <- round((1 - train.prop) * dim(total.set)[1])

train <- sample(c(rep(TRUE, train.size), rep(FALSE, test.size)))
train.coords <- total.set[train, ]
test.coords <- total.set[!train, ]

X.train <- cbind(1, rnorm(train.size))
X.test <- cbind(1, rnorm(test.size))
X.comb <- rbind(X.train, X.test)

B <- as.matrix(c(1, 5))

sigma.sq <- 1
tau.sq <- 0.1
phi <- 12

D <- sqrt(laGP::distance(total.set))
R <- exp(-phi * D)
w <- mdnngp:::rmvn(1, rep(0, train.size + test.size), sigma.sq * R)
y <- rnorm(train.size + test.size, X.comb %*% B + w, sqrt(tau.sq))

w.train <- w[train]
w.test <- w[!train]

y.train <- y[train]
y.test <- y[!train]

interp.plot(res = 200, X = train.coords, Z = w.train, Zlim = c(-3, 3), save = TRUE, 
            filename = "~/git-repos/classlib/rsm/hw/finalproject/p1.png", width = 720,
            height = 720)

full_fit <- spLM(y.train ~ X.train - 1, coords = train.coords, 
                 starting = list("phi" = 10, "sigma.sq" = 10, "tau.sq" = 10), 
                 tuning = list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1), 
                 priors = list("beta.Norm" = list(rep(0, 2), diag(1000, 2)), 
                               "phi.Unif" = c(3, 30), 
                               "sigma.sq.IG" = c(2, 2), 
                               "tau.sq.IG" = c(2, 0.1)),
                 cov.model = "exponential", n.samples = 1000, verbose = FALSE, n.report = 50)
full_rs <- spRecover(full_fit, start = 100, verbose = FALSE)

interp.plot(res = 200, X = train.coords, Z = rowMeans(full_rs$p.w.recover.samples), 
            Zlim = c(-3, 3), filename = "~/git-repos/classlib/rsm/hw/finalproject/p2.png", 
            save = TRUE, width = 720, height = 720)

knot_fit <- spLM(y.train ~ X.train - 1, coords = train.coords, knots = c(8, 8, -0.0625),
                 starting = list("phi" = 10, "sigma.sq" = 10, "tau.sq" = 10), 
                 tuning = list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1), 
                 priors = list("beta.Norm" = list(rep(0, 2), diag(1000, 2)), 
                               "phi.Unif" = c(3, 30), 
                               "sigma.sq.IG" = c(2, 2), 
                               "tau.sq.IG" = c(2, 0.1)),
                 cov.model = "exponential", n.samples = 1000, verbose = FALSE, n.report = 50)
knot_rs <- spRecover(knot_fit, start = 100, verbose = FALSE)

interp.plot(res = 200, X = train.coords, Z = rowMeans(knot_rs$p.w.recover.samples), 
            Zlim = c(-3, 3), filename = "~/git-repos/classlib/rsm/hw/finalproject/p3.png", 
            save = TRUE, width = 720, height = 720)

gpslm_fit <- gpslm(formula = y.train ~ X.train - 1, coords = train.coords, m = 10, col = 2, n.samples = 500,
                   priors = list("phi" = c(3, 30), 
                                 "sigma.sq" = c(2, 2),
                                 "tau.sq" = c(2, 0.1), 
                                 "beta" = c(0, 1000)),
                   starting = list("phi" = 12, "sigma.sq" = 1, "tau.sq" = 0.1),
                   tuning = list("phi" = 0.1, "sigma.sq" = 0.1))

interp.plot(res = 200, X = train.coords, Z = colMeans(gpslm_fit$w.samples), 
            Zlim = c(-3, 3), filename = "~/git-repos/classlib/rsm/hw/finalproject/p4.png", 
            save = TRUE, width = 720, height = 720)

gpslm_fit <- gpslm(formula = y.train ~ X.train - 1, coords = train.coords, m = 20, col = 2, n.samples = 500,
                   priors = list("phi" = c(3, 30), 
                                 "sigma.sq" = c(2, 2),
                                 "tau.sq" = c(2, 0.1), 
                                 "beta" = c(0, 1000)),
                   starting = list("phi" = 12, "sigma.sq" = 1, "tau.sq" = 0.1),
                   tuning = list("phi" = 0.1, "sigma.sq" = 0.1))

interp.plot(res = 200, X = train.coords, Z = colMeans(gpslm_fit$w.samples), 
            Zlim = c(-3, 3), filename = "~/git-repos/classlib/rsm/hw/finalproject/p5.png", 
            save = TRUE, width = 720, height = 720)





