rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  try <- try({D <- chol(V)})
  if(class(try) == "try-error") {
    D <- Matrix::nearPD(V, conv.norm.type = "F")$mat
    D <- chol(V)
  }
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

itol <- function (i, j, n) {
  k <- (2 * n - j) * (j - 1) / 2 + (i - j)
  return(k)
}

ltoi <- function (k, n) {
  j <- floor(((2 * n + 1) - sqrt((2 * n - 1) ^ 2 - 8 * (k - 1))) / 2)
  i <- j + k - (2 * n - j) * (j - 1) / 2
  return(c(i, j))
}
