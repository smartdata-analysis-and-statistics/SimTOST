% Generated manually; source documentation is in R/run_simulations.R
\name{run_simulations}
\alias{run_simulations}
\title{Run a design-specific simulation}
\description{Dispatches to the compiled simulation routine for a selected trial design and test type.}
\usage{run_simulations(design = c("parallel", "2x2"), ctype = c("DOM", "ROM"), ...)}
\arguments{
\item{design}{Trial design, either \code{"parallel"} or \code{"2x2"}.}
\item{ctype}{Test type, either \code{"DOM"} or \code{"ROM"}.}
\item{...}{Arguments passed to the selected simulation routine.}
}
\value{The result returned by the selected simulation routine.}
