% Generated manually; source documentation is in R/simpower_methods.R
\name{simpower-methods}
\alias{print.simss}
\alias{summary.simss}
\alias{confint.simss}
\alias{plot.simss}
\alias{summary.simpower}
\alias{print.simpower}
\alias{confint.simpower}
\alias{plot.simpower}
\title{Methods for unified SimTOST result classes}
\description{Printing, summary, confidence-interval, and plotting methods for the unified \code{simss} and \code{simpower} result classes returned by \code{sampleSize()} and \code{simPower()}.}
\usage{
\method{print}{simss}(x, ...)
\method{summary}{simss}(object, ...)
\method{confint}{simss}(object, ...)
\method{plot}{simss}(x, target_power = 0.80, ...)
\method{print}{simpower}(x, ...)
\method{summary}{simpower}(object, ...)
\method{confint}{simpower}(object, parm, level = 0.95, ...)
\method{plot}{simpower}(x, target_power = 0.80, ...)
}
\arguments{\item{object}{A result returned by \code{simPower()} or \code{sampleSize()}.}\item{x}{A result returned by \code{simPower()} or \code{sampleSize()}.}\item{parm}{Unused compatibility argument.}\item{level}{Requested confidence level.}\item{target_power}{Target power shown as a horizontal reference line. The default is 0.80.}\item{...}{Additional arguments.}}
\value{The summary method returns a data frame; the confidence-interval method returns interval limits; plotting methods return a \code{ggplot} object; and printing methods return the input invisibly.}
