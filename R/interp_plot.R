#' Interpolating Image Plots
#' 
#' Interpolate a potentially non-gridded surface and output as an interpolated raster plot.
#' @usage 
#' interp.plot(X, Z, bounds = c(0, 1), resolution = 100, save = FALSE, ...)
#' @param X a numeric matrix with 2 columns specifying the coordinates at the observed Z values.
#' @param Z a numeric vector of responses corresponding to the coordinates in X.
#' @param Zlim a numeric vector representing the lower and upper limit for Z to use in plotting.
#' @param bounds a vector of length 2 indicating the lower and upper bound for interpolation and plotting.
#' @param resolution an integer, specifying the resolution to interpolate at (number of unique points in each dimension).
#' @param save a logical, indicating whether to save the plot as a file, or output to the current graphical device.
#' @param ... extra parameters, to be passed to the png function
#' @export
#' @import ggplot2
#' @importFrom laGP aGP
#' @importFrom parallel detectCores

interp.plot <- function(bounds = c(0, 1), res = 100, X, Z, Zlim, save = FALSE, ...) {
x <- seq(bounds[1], bounds[2], length.out = res)
XX <- as.matrix(expand.grid(x, x))
fit <- laGP::aGP(X, Z, XX, method = "nn", end = 50, omp.threads = detectCores() / 2, verb = 0)
XXf <- data.frame(XX, fit$mean)
colnames(XXf) <- c("X", "Y", "Z")
if(save) png(...)
out <- ggplot(data = XXf, mapping = aes(x = X, y = Y, fill = Z)) + 
  geom_raster(interpolate = TRUE) + 
  scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"), limits = Zlim) + 
  theme_bw()
if(save) {
  print(out)
  dev.off()
}
return(out)
}