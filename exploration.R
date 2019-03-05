x <- seq(0, 1, length.out = 1000)
xx <- expand.grid(x1 = x, x2 = x)
y <- sin(xx$x1 * 2 * pi) - xx$x2^2

image(x, x, matrix(y, ncol = length(x)))

      