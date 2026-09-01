#' SimTOST: A Package for Sample Size Simulations
#'
#' The SimTOST package provides tools for simulating sample sizes, calculating power,
#' and assessing type-I error for various statistical scenarios.
#'
#' @docType package
#' @name SimTOST
#' @title Sample Size Estimation via Simulation
#' @author
#' Thomas Debray \email{tdebray@fromdatatowisdom.com} (author and maintainer)
#'
#' Other contributors:
#' \itemize{
#'   \item Tim Friede \email{tim.friede@med.uni-goettingen.de} \code{[contributor]}
#'   \item Johanna Munoz \email{johanna.munoz@fromdatatowisdom.com} \code{[contributor]}
#'   \item Dewi Amaliah \email{dewi.amaliah@fromdatatowisdom.com} \code{[contributor]}
#'   \item Wei Wei \email{wei.wei@biogen.com} \code{[contributor]}
#'   \item Marian Mitroiu \email{marian.mitroiu@biogen.com} \code{[contributor]}
#'   \item Scott McDonald \email{scott.mcdonald@fromdatatowisdom.com} \code{[contributor]}
#'   \item Biogen Inc \code{[copyright holder, funder]}
#' }
#' @references
#' Mielke, J., Jones, B., Jilma, B. & König, F. Sample Size for Multiple Hypothesis Testing in Biosimilar Development. \emph{Statistics in Biopharmaceutical Research} 10, 39–49 (2018).
#' @section Planning:
#' Each function name links to its full help page.
#' \itemize{
#'   \item \code{\link{sampleSize}}: simulation-based sample-size planning
#'   for continuous and count outcomes, including multiple endpoints and
#'   comparator families.
#'   \item \code{\link{simPower}}: simulated power for a fixed sample size,
#'   with support for continuous and count-outcome analyses.
#'   \item \code{\link{sampleSize_Mielke}}: Mielke et al.'s sample-size
#'   calculation for multiple, correlated hypotheses and k-out-of-m rules.
#' }
#' @section Continuous outcomes:
#' \itemize{
#'   \item \code{\link{simParallelEndpoints}}: generates correlated normal
#'   or log-normal endpoint data for a parallel-group design.
#'   \item \code{\link{get_par}}: prepares and validates endpoint,
#'   covariance, allocation, and hierarchy parameters for planning.
#'   \item \code{\link{run_simulations}}: dispatches a continuous
#'   simulation to the selected design and test combination.
#' }
#' @section Simulation-result methods:
#' \itemize{
#'   \item \code{\link{print.simss}} (print): prints a concise design, power,
#'   confidence-interval, and sample-size report.
#'   \item \code{\link{summary.simss}} (summary): returns and prints a
#'   structured summary of the simulation result.
#'   \item \code{\link{confint.simss}} (confint): extracts the stored Monte
#'   Carlo confidence interval for achieved power.
#'   \item \code{\link{plot.simss}} (plot): plots simulated power against
#'   sample size with confidence intervals and the target-power line.
#'   \item \code{\link{update.simss}} and \code{\link{update.simpower}}
#'   (update): reruns a result while replacing only explicitly supplied
#'   planning or simulation parameters.
#' }
"_PACKAGE"
