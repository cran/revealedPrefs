################################################################################
################################################################################
## Rationality axioms functions
library(Rcpp)

################################################################################
## WARP

## Function to check WARP with exact algorithm (check all pairs)
## violation if (x_i p_i >= x_k p_i) AND (x_k p_k >= x_i p_k) AND (x_i != x_k)
checkWarp <- function(x, p) {
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  the.call <- .Call("CheckWarp", x, p, PACKAGE= "revealedPrefs")
  the.call$type <- "WARP"
  class(the.call) <- "axiomTest"
  return(the.call)
}

################################################################################
## SARP

## Function to check SARP:
##  - depth-first search with tabu list
##  - floyd-warshall 
## SARP violated if slack preference cycle of unequal quantities
checkSarp <- function(x, p, method= c("deep", "floyd")) {
  method <- match.arg(method)
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  
  if (method == "floyd") {
    res <- .Call("CheckSarp", x, p, PACKAGE= "revealedPrefs")
  } else {
    res <- .Call("DeepSarp", x, p, PACKAGE= "revealedPrefs")
    if (res$violation) {
      res$path <- res$path + 1
      cycle.start <- 
        which(res$path == res$path[length(res$path)])[1]
      res$violators <- 
        res$path[cycle.start:(length(res$path) - 1)]
      
      if (length(res$violators) == 2) {
        res$direct.violation <- TRUE
      } else res$direct.violation <- FALSE
    }
  }
  res$type <- "SARP"
  res$method <- method
  class(res) <- "axiomTest"
  return(res)
}

################################################################################
## GARP

## Function to check GARP with exact algorithm
## (Warshall-Floyd or depth-first tabu)
## GARP violated if strict cycle present
## for quantities x and prices p
checkGarp <- function(x, p, method= c("deep", "floyd")){
  method <- match.arg(method)
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  x <- as.matrix(x)
  p <- as.matrix(p)

  if (method == "floyd") {
    the.call <- .Call("CheckGarp", p%*%t(x), PACKAGE = "revealedPrefs")
  } else {
    the.call <- .Call("DeepGarp", x, p, sample(1:nrow(x)) - 1, 
                      PACKAGE= "revealedPrefs")
    if (the.call$violation) {
      the.call$path <- the.call$path + 1
      cycle.start <- 
        which(the.call$path == the.call$path[length(the.call$path)])[1]
      the.call$violators <- 
        the.call$path[cycle.start:(length(the.call$path) - 1)]
      the.call$strict <- 
        the.call$path.strict[cycle.start:length(the.call$path.strict)]
      
      if (length(the.call$violators) == 2) {
        the.call$direct.violation <- TRUE
      } else the.call$direct.violation <- FALSE
    }
  }
  the.call$type <- "GARP"
  the.call$method <- method
  class(the.call) <- "axiomTest"
  the.call
}


################################################################################
################################################################################
## S3 methods

print.axiomTest <- function(x, ...) {
  cat("  Axiomatic rationality test:", x$type, 
      ifelse(x$violation, "violation found.\n", "no violation.\n"))
}

summary.axiomTest <- function(object, ...) {
  cat("\n ", object$type, "rationality test:", 
      ifelse(object$violation, "violation found.\n", "no violation.\n"))
  
  cat("  Method:")
  if (object$type == "WARP") {
    cat(" Pairwise comparisons.\n")
  } else {
    if (object$method == "floyd") cat(" Floyd-Warshall algorithm.\n")
    if (object$method == "deep") cat(" Depth-first search with tabu list.\n")
  }

  if (object$violation) {
    cat("\n")
    if (object$type == "WARP") {
      cat("  Violating observations:", 
          paste(object$violators, rep(">=", 2), collapse= " "), 
          object$violators[1], "\n")
      cat("                        : (direct preferences)\n") 
      cat("                        :", 
          object$violators[1], "!=", object$violators[2], "\n")

      cat("\n  Other axioms:\n")
      cat("  * SARP      : violated (symmetry of direct preferences,", 
          "unequal quantities).\n")
      cat("  * GARP      : unknown (strict preferences not computed).\n")
    }
    
    if (object$type == "SARP") {
      cat("  Violating observations:", 
          paste(object$violators, ">=", collapse= " "), 
          object$violators[1], "\n")
      if (object$direct.violation | object$method == "deep") {
        cat("                        : (direct preferences)\n")
      } else cat("                        : (indirect preferences)\n")
      cat("                        : And not all quantities in cycle equal.\n")
      
      cat("\n  Other axioms:\n")
      cat("  * WARP      :")
      if (object$direct.violation) {
        cat(" violated (symmetry of direct preferences, unequal quantities).\n")
      } else {
        # In Floyd algorithm all direct violations are tested before indirect
        if (object$method == "floyd") 
          cat(" not violated (all pairwise checks passed).\n")
        # In depth-first search direct violations are not all checked first
        if (object$method == "deep") 
          cat(" unknown (stopped before all pairwise checks conducted).\n")
      }
      cat("  * GARP      : unknown (strict preferences not computed).\n")
    }

    if (object$type == "GARP") {
      signs <- rep(">=", )
      signs[object$strict] <- ">"
      cat("  Violating observations:", 
          paste(object$violators, signs, collapse= " "), 
          object$violators[1], "\n")
      if (object$method == "deep" | object$direct.violation) {
        cat("                        : (direct preferences)\n")
      } else cat("                        : (indirect preferences)\n")
      
      cat("\n  Other axioms:\n")
      cat("  * WARP      :")
      if (object$direct.violation) {
        cat(" violated (symmetry of direct preferences, unequal quantities).\n")
      } else cat(" unknown (equality of quantities not tested).\n")
      cat("  * SARP      : violated (symmetry of indirect preferences,",
          "unequal quantities).\n")
    }    
  } else { # If no violation detected
    cat("\n  Other axioms:\n")
    if (object$type == "WARP") {
      cat("  * SARP      : unknown (indirect preferences not computed).\n")
      cat("  * GARP      : unknown (strict preferences not computed).\n")
    }
    if (object$type == "SARP") {
      cat("  * WARP      : not violated.\n")
      cat("  * GARP      : unknown (strict preferences not computed).\n")
    }
    if (object$type == "GARP") {
      cat("  * WARP      : unknown (equality of quantities not tested).\n")
      cat("  * SARP      : unknown (equality of quantities not tested).\n")
    }
  }
  cat("\n")
}
