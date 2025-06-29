% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.simss.R
\name{summary.simss}
\alias{summary.simss}
\title{Summary for Simulation Results}
\usage{
\method{summary}{simss}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"simss"} returned by a sampleSize function.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A named numeric vector with the sample size per arm and the total (Total) sample size.
}
\description{
Generates a summary of the simulation results, including per-arm and total sample sizes.
}
\examples{
# Assume `res` is a result from `sampleSize()`
# summary(res)
}
