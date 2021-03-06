% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeToPoly.R
\name{treeToPoly}
\alias{treeToPoly}
\title{Convert trees to coefficient matrices}
\usage{
treeToPoly(
  trees,
  type = c("default", "yEvaluated", "tipLabel"),
  y,
  varLabels = FALSE,
  numThreads = -1
)
}
\arguments{
\item{trees}{a single phylo object or a list of phylo objects}

\item{type}{one of:
\describe{
  \item{\dQuote{real}}{tree distinguishing polynomials in two variables x (columns) and y (rows)}
  \item{\dQuote{yEvaluated}}{tree distinguishing polynomials with y evaluated at a specified argument}
  \item{\dQuote{tipLabel}}{complex coefficient polynomial that utilize binary trait tip labels on the phylo objects}
}}

\item{y}{the y value to evaluate the polynomial at when type is \dQuote{yEvaluated}, ignored otherwise}

\item{varLabels}{boolean for whether to add row and column names corresponding to the variables in the polynomial}

\item{numThreads}{number of threads to be used, the default (-1) will use the number of cores in the machine and numThreads = 0 will only use the main thread}
}
\value{
the resulting coefficient matrix or matrices of the form:
\describe{
  \item{\dQuote{real}}{a real matrix where the ith row, jth column represents the x^(j-1)*y^(i-1) coefficient}
  \item{\dQuote{yEvaluated}}{a vector where the kth column represents the x^(k-1) coefficient}
  \item{\dQuote{tipLabel}}{given trees with two unique tip labels \dQuote{a}, \dQuote{b} a complex matrix where the ith row, jth column represents the a^(i-1)*b^(j-1) coefficient}
}
}
\description{
Converts rooted full binary trees to tree distinguishing polynomials described with coefficient matrices.
}
\examples{
library(treenomial)
library(ape)

# generate a tree
tree <- rtree(n = 30, rooted = TRUE)

# a real coefficient matrix
treeToPoly(tree, varLabels = TRUE, numThreads = 0)

# complex coefficient vector for the tree
treeToPoly(tree, type = "yEvaluated", y = 1+1i, varLabels = TRUE, numThreads = 0)

# for a list of trees
treeToPoly(rmtree(4, 20), varLabels = TRUE, numThreads = 0)

}
