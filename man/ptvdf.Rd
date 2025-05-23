% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ptvdf}
\alias{ptvdf}
\title{Calculate p-values using t-distribution with Variable Degrees of Freedom}
\usage{
ptvdf(x, df, lower)
}
\arguments{
\item{x}{arma::mat (vector) - A matrix or vector of random variable values for which the p-values will be calculated.}

\item{df}{arma::mat (vector) - A matrix or vector of degrees of freedom for the t-distribution, matching the size of \code{x}.}

\item{lower}{bool - If \code{TRUE}, calculates the lower-tail probability (P(T <= x)); if \code{FALSE}, calculates the upper-tail probability.}
}
\value{
arma::mat (vector) - A matrix containing the computed cumulative distribution function (p-values) for each element in \code{x}. The result is returned as a 1xN matrix, where N is the number of elements in \code{x}.
}
\description{
This function computes the cumulative distribution function (p-values) for a given random variable \code{x} and corresponding degrees of freedom \code{df} using the t-distribution. The function can compute the lower or upper tail probabilities depending on the value of the \code{lower} argument.
}
\keyword{internal}
