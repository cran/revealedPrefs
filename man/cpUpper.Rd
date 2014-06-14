\name{cpUpper}
\alias{cpUpper}
\alias{print.upperBound}
\alias{summary.upperBound}

\title{
Upper bound on the number of GARP-consistent subpopulations and clustering.
}
\description{
The \code{cpUpper} function computes a Crawford-Pendakur type upper bound on the number of GARP-consistent subpopulations, and performs clustering of the data.
}
\usage{
cpUpper(x, p, times= 1, method= c("fastfloyd", "deep", "floyd"))
\method{print}{upperBound}(x, ...)
\method{summary}{upperBound}(object, ...)
}
\arguments{
  \item{x}{data frame or matrix containing the observed quantities, where each row corresponds to an observation and the columns are types of goods, or an object of class \code{upperBound} to be used with \code{print},}
  \item{p}{data frame or matrix (of same dimensions as x) containing the corresponding prices,}
  \item{times}{number of times the algorithm is run (the final result is the best of \code{times} results,  ie. lowest number of clusters found),}
  \item{method}{character string: \code{"fastfloyd"} or \code{"deep"} or \code{"floyd"} (see Details section). Default \code{"fastfloyd"}. }
  \item{object}{an object of class \code{upperBound} as returned by \code{cpUpper},}
  \item{...}{additional arguments passed to the \code{print} and 
  \code{summary} methods (unused).}
}
\value{
\code{cpUpper} returns an object of class \code{upperBound} which contains the following elements:
  \item{\code{n.types}}{upper bound on the number of types,}
  \item{\code{clustering}}{numeric vector with length equal to number of observations, containing the cluster number of each observation,}
  \item{\code{cluster.pop}}{numeric vector containg the numbers of observations allocated to each cluster,}
  \item{\code{hist.n.types}}{numeric vector containing the history of numbers of clusters found during multiple runs of the algorithm.}
}
\details{
For each run of the algorithm, a random permutation of the observations is drawn, and one by one each observation is associated with the biggest cluster that can include it while preserving consistency with the GARP rationality axiom. If no cluster is compatible with a given observation a new cluster is created to accomodate it.

Three methods are available: 
\code{"fastfloyd"}{ (default) uses an iterative variant of the Floyd-Warshall algorithm, in which the check of consistency of the current observation with a given cluster is done in a single step of the Floyd-Warshall algorithm. Much faster than \code{"floyd"}.}

\code{"deep"}{ uses a single run of recursive depth-first search with tabu list for each check of an observation against a given cluster. Faster than \code{"fastfloyd"} on large datasets (eg. > 5000 observations).}

\code{"floyd"}{ uses the algorithm described in Crawford and Pendakur (2013): at each step the complete Floyd-Warshall algorithm is run to check whether each cluster can accomodate the current observation. Much slower than the two other algorithms.}
}
\note{
Warning: \code{cpUpper} can be very slow for large datasets (eg. more than 1000 observations).
}
\references{
Crawford, I. and Pendakur, K. (2013). How many types are there?
\emph{The Economic Journal}, 123(567):77-95.
}
\author{
Julien Boelaert \email{jubo.stats@gmail.com}
}
\seealso{See \link{cpLower} for the lower bound on number of types.}
\examples{
# Cluster GARP-violating data:
data(noGarp)
cp.up <- cpUpper(noGarp$x, noGarp$p)
cp.up$clustering
}

\keyword{methods}
