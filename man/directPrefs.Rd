\name{directPrefs}
\alias{directPrefs}
\alias{indirectPrefs}

\title{
Compute direct and indirect revealed preferences.
}
\description{
The \code{directPrefs} function computes direct revealed preferences, the \code{indirectPrefs} function computes indirect revealed preferences.
}
\usage{
directPrefs(x, p)
indirectPrefs(x, p)
}
\arguments{
  \item{x}{data frame or matrix containing the observed quantities, where each row corresponds to an observation and the columns are types of goods.}
  \item{p}{data frame or matrix (of same dimensions as x) containing the corresponding prices.}
}
\value{
A matrix of direct or indirect revealed preferences:
  \item{\code{prefs[i, j] = 0}}{if bundle \code{i} is not revealed prefered to bundle \code{j}} 
  \item{\code{prefs[i, j] = 1}}{if bundle \code{i} is revealed prefered to bundle \code{j}} 
  \item{\code{prefs[i, j] = 2}}{if bundle \code{i} is strictly revealed prefered to bundle \code{j}} 
}
\details{
Direct preferences are directly computed from matrix multiplication of prices and quantities. Indirect preferences are computed using a variant of the Floyd-Warshall algorithm.
}
\references{
Varian, H. R. (1984) \emph{Microeconomic Analysis}. New York/London:
Norton, 2nd edition, pp 141-143.
}
\author{
Julien Boelaert \email{jubo.stats@gmail.com}
}

\examples{
# Compute direct and indirect preferences of SARP-violating data:
data(noSarp)
directPrefs(noSarp$x, noSarp$p)
indirectPrefs(noSarp$x, noSarp$p)
}
