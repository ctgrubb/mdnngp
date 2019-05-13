library(data.table)

if(Sys.info()["nodename"] == "broce") {                                          
  sourcefolder <- "~/git-repos/mdnngp/R/"
  datafolder <- "/storage/chris/mdnngp/"
} else {                                                                         
  stop("Unknown host detected. Please declare folder locations.")                
}

nldas_data <- fread(paste0(datafolder, "nldas_data.csv"), header = TRUE, sep = ",")
modis_data <- fread(paste0(datafolder, "modis_data.csv"), header = TRUE, sep = ",")
modis_time <- fread(paste0(datafolder, "modis_time.csv"), header = TRUE, sep = ",")

modis_data <- cbind(modis_data, modis_time[, -c(1, 2, 3)])

nldas_df <- melt(nldas_data, id.vars = c(1, 2, 3), variable.name = "DateTime", 
                 variable.factor = FALSE, value.name = "Temperature")[1:(384 * 366), ]
modis_df <- melt(modis_data, id.vars = c(1, 2, 3), measure.vars = list(60:115, 4:59), 
                 variable.name = "Date", variable.factor = FALSE, value.name = c("Time", "Temperature"))

nldas_df$DateTime <- as.POSIXct(nldas_df$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
nldas_df$Hour <- as.numeric(format(nldas_df$DateTime, format = "%H"))
nldas_df$Day <- as.numeric(floor(difftime(nldas_df$DateTime, nldas_df$DateTime[1], units = "days")))

modis_df$DateTime <- as.POSIXct(paste0(rep(colnames(modis_time)[-c(1, 2, 3)], each = nrow(modis_data)), modis_df$Time), 
                                format = "%Y-%m-%d %H:%M", tz = "EST")
modis_df$Hour <- as.numeric(format(modis_df$DateTime, format = "%H", tz = "GMT")) + 
  as.numeric(format(modis_df$DateTime, format = "%M", tz = "GMT")) / 60
modis_df$Day <- as.numeric(floor(difftime(modis_df$DateTime, nldas_df$DateTime[1], units = "days")))

nldas_matrix <- as.matrix(nldas_df[, c(1, 2, 6, 7)])
modis_matrix <- as.matrix(modis_df[, c(1, 2, 8, 9)])

encode <- function(x) {
  ranges <- matrix(NA, nrow = ncol(x), ncol = 2)
  for(i in 1:ncol(x)) {
    ranges[i, ] <- range(x[, i], na.rm = TRUE)
    x[, i] <- (x[, i] - ranges[i, 1]) / (ranges[i, 2] - ranges[i, 1])
  }
  out <- list(X = x, r = ranges)
}

enc <- encode(rbind(nldas_matrix, modis_matrix))
X_n <- enc$X[1:nrow(nldas_matrix), ]
X_m <- enc$X[-(1:nrow(nldas_matrix)), ]
X_m <- X_m[complete.cases(X_m), ]
Y_n <- as.numeric(nldas_df$Temperature)
Y_m <- as.numeric(modis_df$Temperature[complete.cases(modis_df)])

library(laGP)

x <- seq(0, 1, length.out = 5)
ref <- as.matrix(expand.grid(x, x, .7, x))
ref2 <- as.matrix(expand.grid(x, x, x, x))
ref3 <- as.matrix(expand.grid(x, x, x, 0.08275862))

test1 <- aGPsep(X = X_m, Z = Y_m, XX = ref, omp.threads = 20, verb = 0)
test2 <- aGPsep(X = X_n, Z = Y_n, XX = ref2, omp.threads = 20, verb = 0)
test3 <- aGPsep(X = rbind(X_m, X_n), Z = c(Y_m, Y_n), XX = ref2, omp.threads = 20, verb = 0)

test4 <- aGPsep(X = X_n[X_n[, 4] < 0.09, ], Z = Y_n[X_n[, 4] < 0.09], XX = ref3, omp.threads = 20, verb = 0)

library(orthopolynom)

Legendre <- function(x, n, normalized = TRUE, intercept = FALSE, rescale = TRUE) {
  ## Create a design matrix for Legendre polynomials
  ## x - numeric
  ## n - see orthopolynom
  ## normalized - logical, see orthopolynom
  ## intercept - logical, add intercept
  tmp <- legendre.polynomials(n = n, normalized = normalized)
  if(!intercept) tmp <- tmp[2:length(tmp)]
  polynomial.values(polynomials = tmp, x = x, matrix = TRUE)
}

polynomial.values <- function(polynomials, x, matrix = FALSE) {
  ## Changed copy of polynomial.vales from orthopolynom in order
  ## to add matrix argument
  require(polynom)
  n <- length(polynomials)
  if(!matrix) {
    values <- vector(mode = "list", length = n)
  } else {
    values <- matrix(ncol = n, nrow = length(x))
  }
  j <- 1
  while(j <= n) {
    if(!matrix) {
      values[[j]] <- predict(polynomials[[j]], x)
    } else {
      values[, j] <- predict(polynomials[[j]], x)
    }
    j <- j + 1
  }
  values
}


lMref <- data.frame(cbind(Legendre(ref2[, 1], 3), Legendre(ref2[, 2], 3), Legendre(ref2[, 3], 3), Legendre(ref2[, 4], 3)))
colnames(lMref) <- paste0("L", 1:12)

lM <- data.frame(Y_n, cbind(Legendre(X_n[, 1], 3), Legendre(X_n[, 2], 3), Legendre(X_n[, 3], 3), Legendre(X_n[, 4], 3)))
colnames(lM) <- c("Y", paste0("L", 1:12))
test <- lm(Y ~ ., data = lM)
range(predict(test, lMref))

