\name{revealedPrefs}
\alias{revealedPrefs}

\docType{package}
\title{
Revealed preferences and microeconomic rationality
}
\description{
Computation of (direct and indirect) revealed preferences, fast non-parametric tests of rationality axioms (WARP, SARP, GARP) and detection of GARP-consistent subpopulations.
}
\details{
\tabular{ll}{
Package: \tab revealedPrefs\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-06-14\cr
License: \tab GPL (>=2)\cr
}

This package is meant for the analysis of (quantity, price) data, eg. of bundles of goods and corresponding prices. It features fast algorithms that make the analysis of large datasets feasible.

Functions directPrefs and indirectPrefs compute revealed preferences.

Functions checkWarp, checkSarp, checkGarp perform fast non-parametric tests of rationality using the corresponding rationality axioms.

Functions cpLower and cpUpper generate Crawford-Pendakur type bounds on the number of subpopulations and provide the corresponding clusterings.
}
\author{Julien Boelaert \email{jubo.stats@gmail.com}
}
\references{
Varian, H. R. (1982) The Nonparametric Approach to Demand Analysis, 
\emph{Econometrica}, 50(4):945-973.

Crawford, I. and Pendakur, K. (2013). How many types are there?
\emph{The Economic Journal}, 123(567):77-95.
}
\keyword{ package }
\seealso{
See \link{directPrefs} for computation of preferences, \link{checkGarp} for rationality tests, and \link{cpUpper} for clustering of data into non-violating subsets.
}
\examples{
# Compute preferences and check rationality on a GARP-violating dataset:
data(noGarp)
indirectPrefs(noGarp$x, noGarp$p)
checkGarp(noGarp$x, noGarp$p)

# Cluster dataset into GARP-consistent subpopulations:
cpUpper(noGarp$x, noGarp$p)
}
