% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.r
\name{print.simss}
\alias{print.simss}
\title{Print Summary of Sample Size Estimation}
\usage{
\method{print}{simss}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"simss"}, typically generated by a sample size estimation function.}

\item{...}{Optional arguments to be passed from or to other methods.}
}
\value{
No return value, called for side effects. The function prints the summary results
of the sample size estimation to the console in a structured format.
}
\description{
Prints the summary results of the sample size estimation for bioequivalence trials,
including achieved power, total sample size, and power confidence intervals.
The function also details the study design, primary endpoint comparisons,
and applied multiplicity corrections.
}
\details{
This function displays key metrics from a sample size estimation analysis.
It provides an overview of the study design, treatment comparisons,
tested endpoints, significance level adjustments, and estimated sample size.
For studies with multiple primary endpoints, it describes the multiplicity correction applied.
}
\author{
Thomas Debray \email{tdebray@fromdatatowisdom.com}
}
