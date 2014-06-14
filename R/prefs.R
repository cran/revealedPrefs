################################################################################
################################################################################
## Direct and indirect preferences
library(Rcpp)

## Function to compute all direct preferences
## prefs[i, j]= 0 iff i not prefered to j (p_i q_i < p_i q_j)
## prefs[i, j]= 1 iff equality with prices i (p_i q_i == p_i q_j)
## prefs[i, j]= 2 iff i strictly prefered to j (p_i q_i > p_i q_j)
directPrefs <- function(x, p) {
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  px <- p %*% t(x)
  prefs <- matrix(0, nrow= nrow(x), ncol= nrow(x))
  prefs[diag(px) == px] <- 1
  prefs[diag(px) > px] <- 2
  
  prefs
}

## Function to compute all indirect preferences
## prefs[i, j]= 0 iff i not indirectly prefered to j
## prefs[i, j]= 1 iff i indirectly prefered to j (only equalities)
## prefs[i, j]= 2 iff i indirectly strictly prefered to j 
##                    (with at least one strict preference)
indirectPrefs <- function(x, p) {
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  
  .Call("IndirectPrefs", p %*% t(x), PACKAGE= "revealedPrefs")
}
