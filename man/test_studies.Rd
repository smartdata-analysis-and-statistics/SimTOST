% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{test_studies}
\alias{test_studies}
\title{test_studies}
\usage{
test_studies(nsim, n, comp, param, param.d, arm_seed, ncores)
}
\arguments{
\item{nsim}{number of simulated studies}

\item{n}{sample size}

\item{comp}{index comparator}

\item{param}{list of parameters (mean,sd,tar)}

\item{param.d}{design parameters}

\item{arm_seed}{seed for each endpoint to get consistent in simulations across all comparators}

\item{ncores}{number of cores used for the calculation}
}
\value{
a logical matrix of size  (nsim) X (number of endpoints + 1) function only replicates test_bioq nsim times.
}
\description{
Internal function to estimate the bioequivalence test for nsim simulated studies given a sample size n
}
\keyword{internal}
