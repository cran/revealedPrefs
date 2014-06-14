################################################################################
################################################################################
## Crawford-Pendakur Upper and Lower Bound functions
library(Rcpp)

## Function for Crawford-Pendakur-type Upper Bound algorithms
## for quantities x and prices p
cpUpper <- function(x, p, times= 1, method= c("fastfloyd", "deep", "floyd")) {
  method <- match.arg(method)
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  
  the.samples <- replicate(times, sample(0:(nrow(x)-1)))
  
  if (method == "floyd") {
    the.call <- .Call("CpUp", p%*%t(x), the.samples, PACKAGE = "revealedPrefs")
  } else if (method == "fastfloyd") {
    the.call <- .Call("FastUp", p%*%t(x), the.samples, 
                      PACKAGE = "revealedPrefs")
  } else 
    the.call <- .Call("DeepCpUp", x, p, the.samples, 
                      PACKAGE = "revealedPrefs")
  
  the.call$clustering <- as.numeric(the.call$clustering)
  the.call$cluster.pop <- 
    as.numeric(the.call$cluster.pop[the.call$cluster.pop != 0])
  class(the.call) <- "upperBound"
  the.call$n.types <- length(the.call$cluster.pop)
  the.call
}

## Function for Crawford-Pendakur-type Lower Bound algorithms
## for quantities x and prices p
cpLower <- function(x, p, times= 1) {
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  
  the.samples <- replicate(times, sample(0:(nrow(x)-1)))
  the.call <- .Call("CpLow", x, p, the.samples, PACKAGE= "revealedPrefs")
  the.call$violators <- the.call$violators + 1
  the.call$n.types <- length(the.call$violators)
  class(the.call) <- "lowerBound"
  the.call
}

################################################################################
################################################################################
## S3 methods

print.lowerBound <- function(x, ...) {
  cat("  Lower bound on the number of types :", x$n.types, "\n")
  cat("  (best of", length(x$hist.n.types), "run(s) of the algorithm)\n")
}

print.upperBound <- function(x, ...) {
  cat("  Upper bound on the number of types :", x$n.types, "\n")
  cat("  Cluster populations                :", x$cluster.pop, "\n")
  cat("  (best of", length(x$hist.n.types), "run(s) of the algorithm)\n")
}

summary.lowerBound <- function(object, ...) print(object)
summary.upperBound <- function(object, ...) print(object)
