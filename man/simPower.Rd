% Primary fixed-sample-size power interface
\name{simPower}
\alias{simPower}
\title{Estimate power at a fixed sample size}
\description{Calculates simulated power for a prespecified sample size.}
\usage{
simPower(n, distribution = c("norm", "lnorm", "pois", "nbinom"),
         mu_list = NULL, varcov_list = NA, sigma_list = NA, cor_mat = NA,
         sigmaB = 0, rate_list = NULL, exposure = 1, dispersion = 0.1,
         Eper = c(0, 0), Eco = c(0, 0), rho = 0, TAR = NULL,
         arm_names = NA, ynames_list = NA, type_y = NA,
         list_comparator = NA, list_y_comparator = NA, alpha = 0.05,
         lequi.tol = NA, uequi.tol = NA, list_lequi.tol = NA,
         list_uequi.tol = NA, dtype = "parallel", ctype = "ROM",
         vareq = TRUE, k = NA, adjust = "no", dropout = NA,
         nsim = 5000, seed = 1234, ncores = 1,
         keep_sim_data = FALSE, .warn_redundant_bon = TRUE)
}
\arguments{
\item{n}{Integer sample size or vector of sample sizes. For parallel designs this is the base sample size used to derive arm sizes; for 2x2 designs it is per sequence. A vector returns a power curve across base sample sizes.}
\item{distribution}{Outcome distribution using R's names: \code{"norm"},
\code{"lnorm"}, \code{"pois"}, or \code{"nbinom"}. Longer labels such as
\code{"normal"}, \code{"lognormal"}, and \code{"poisson"} are accepted for
compatibility.}
\item{mu_list, varcov_list, sigma_list}{Continuous-outcome inputs.}
\item{cor_mat}{Endpoint correlation matrix, including for joint count simulations.}
\item{sigmaB, Eper, Eco}{Crossover-design parameters.}
\item{rate_list}{Named arm-rate list for count outcomes.}
\item{exposure, dispersion}{Count exposure and negative-binomial dispersion.}
\item{list_comparator, list_y_comparator}{A list of treatment-reference comparisons and
endpoint selections. Each comparator must contain exactly two arm names in the form
\code{c(test, reference)}. The first arm is always the test arm and the second arm
is always the reference arm. This ordering determines the direction of the estimand:
for \code{ctype = "DOM"}, it is test minus reference, and for \code{ctype = "ROM"},
it is test divided by reference. The same convention is used for count outcomes and
for all comparisons shown by the summaries and plots.}
\item{alpha, lequi.tol, uequi.tol, list_lequi.tol, list_uequi.tol}{Significance level and equivalence margins.}
\item{dtype, ctype, vareq, k, adjust, dropout}{Design and testing settings.
For a \code{k}-of-\code{m} rule, \code{adjust = "t"} applies Mielke's strong
k-out-of-m level \code{alpha / (m - k + 1)} separately within each
comparator's selected endpoint family. Legacy \code{"pc"} and
\code{"partial-conjunction"} labels are accepted as aliases.}
\item{nsim, seed, ncores}{Simulation count, seed, and number of cores. For count outcomes, trials are split into independent seeded chunks and each chunk is evaluated by the C++ count kernel.}
\item{keep_sim_data}{Logical. If \code{TRUE}, retain model-scale observations
for each simulated trial.}
\item{.warn_redundant_bon}{Logical. If \code{TRUE}, issue warnings for
redundant or uncalibrated adjustment configurations.}
}
\value{An object of class \code{simpower} containing estimated power and its Monte Carlo confidence interval.}
\details{Use \code{summary()}, \code{confint()}, and \code{plot()} to inspect the result.
When the order of the arms matters, specify \code{list_comparator} explicitly rather
than relying on automatically generated comparisons. For example,
\code{c("T", "R")} estimates \code{T - R} for DOM and \code{T / R} for ROM;
reversing the vector estimates the reverse comparison. The effective \code{m}
is the number of endpoints actually selected for each comparator.}
\examples{
simPower(n = 100, distribution = "Poisson",
         rate_list = list(TEST = .21, REF = .20),
         list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
         list_lequi.tol = list(TEST_vs_REF = .80),
         list_uequi.tol = list(TEST_vs_REF = 1.25),
         exposure = 10, nsim = 100, seed = 1)
}
