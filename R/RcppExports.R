# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Compute p-values for a t-distribution with Fixed Degrees of Freedom
#'
#' @description Computes p-values for a given set of random variables under a t-distribution with fixed degrees of freedom.
#'
#' @param x A numeric matrix (or vector) representing the random variables.
#' @param df A double specifying the degrees of freedom.
#' @param lower A logical value indicating whether to compute the lower-tail probability (\code{P(T <= x)}). If \code{FALSE}, the function returns the upper-tail probability (\code{P(T > x)}).
#'
#' @return A numeric matrix containing the computed cumulative distribution function (CDF) values (p-values).
#'
#' @keywords internal
#' @export
ptv <- function(x, df, lower) {
    .Call(`_SimTOST_ptv`, x, df, lower)
}

#' @title Calculate p-values using t-distribution with Variable Degrees of Freedom
#'
#' @description This function computes the cumulative distribution function (p-values) for a given random variable \code{x} and corresponding degrees of freedom \code{df} using the t-distribution. The function can compute the lower or upper tail probabilities depending on the value of the \code{lower} argument.
#'
#' @param x arma::mat (vector) - A matrix or vector of random variable values for which the p-values will be calculated.
#' @param df arma::mat (vector) - A matrix or vector of degrees of freedom for the t-distribution, matching the size of \code{x}.
#' @param lower bool - If \code{TRUE}, calculates the lower-tail probability (P(T <= x)); if \code{FALSE}, calculates the upper-tail probability.
#'
#' @return arma::mat (vector) - A matrix containing the computed cumulative distribution function (p-values) for each element in \code{x}. The result is returned as a 1xN matrix, where N is the number of elements in \code{x}.
#'
#' @keywords internal
#' @export
ptvdf <- function(x, df, lower) {
    .Call(`_SimTOST_ptvdf`, x, df, lower)
}

#' @title Check Equivalence for Multiple Endpoints
#'
#' @description This function evaluates whether equivalence criteria are met based on a predefined set of endpoints. It first checks whether all primary endpoints satisfy equivalence (if sequential testing is enabled). Then, it determines whether the required number of endpoints (\code{k}) meet the equivalence threshold. The function returns a binary decision indicating whether overall equivalence is established.
#'
#' @param typey An integer vector specifying the hierarchy of each endpoint, where \code{1} denotes a primary endpoint and \code{2} denotes a secondary endpoint.
#' @param adseq A boolean flag indicating whether sequential testing is enabled. If set to \code{TRUE}, all primary endpoints must pass equivalence before secondary endpoints are evaluated. If set to \code{FALSE}, primary and secondary endpoints are assessed independently.
#' @param tbioq A matrix containing the equivalence test results for each endpoint, where \code{1} indicates that equivalence is met and \code{0} indicates that equivalence is not met.
#' @param k An integer specifying the minimum number of endpoints required for overall equivalence.
#'
#' @details When sequential testing is enabled (\code{adseq = TRUE}), all primary endpoints must meet equivalence before secondary endpoints are considered. If sequential testing is disabled (\code{adseq = FALSE}), all endpoints are evaluated simultaneously without hierarchical constraints. The function then determines whether at least \code{k} endpoints meet the equivalence criteria. If the conditions are satisfied, the final equivalence decision (\code{totaly}) is \code{1}; otherwise, it is \code{0}.
#'
#' @return Returns a (1 × 1 matrix) containing a binary equivalence decision. A value of \code{1} indicates that equivalence is established, while \code{0} indicates that equivalence is not established.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @keywords internal
#' @export
check_equivalence <- function(typey, adseq, tbioq, k) {
    .Call(`_SimTOST_check_equivalence`, typey, adseq, tbioq, k)
}

#' @title Simulate a 2x2 Crossover Design and Compute Difference of Means (DOM)
#'
#' @description Simulates a two-sequence, two-period (2x2) crossover design and evaluate equivalence for the difference of means (DOM).
#'
#' @param n integer number of subjects per sequence
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaW matrix  within subject covar-variance matrix across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param sigmaB double between subject variance (assumed same for all endpoints)
#' @param dropout vector of size 2 with dropout proportion per sequence (0,1)
#' @param Eper vector of size 2 with period effect on period (0,1)
#' @param Eco vector of size 2 with carry over effect of arm c(Reference, Treatment).
#' @param typey vector with positions of primary endpoints
#' @param adseq boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seed seed for the simulation
#'
#' @return A numeric matrix containing the simulated hypothesis test results.
#' The first column represents the overall equivalence decision, where 1 indicates
#' success and 0 indicates failure. The subsequent columns contain the hypothesis
#' test results for each endpoint, followed by mean estimates for the reference and
#' treatment groups, and standard deviations for the reference and treatment groups.
#'
#' @export
test_2x2_dom <- function(n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed) {
    .Call(`_SimTOST_test_2x2_dom`, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed)
}

#' @title Simulate a 2x2 Crossover Design and Compute Ratio of Means (ROM)
#'
#' @description Simulates a two-sequence, two-period (2x2) crossover design and evaluate equivalence for the ratio of means (ROM).
#'
#' @param n integer number of subjects per sequence
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaW matrix  within subject covar-variance matrix across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param sigmaB double between subject variance (assumed same for all endpoints)
#' @param dropout vector of size 2 with dropout proportion per sequence (0,1)
#' @param Eper vector of size 2 with period effect on period (0,1)
#' @param Eco vector of size 2 with carry over effect of arm c(Reference, Treatment).
#' @param typey vector with positions of primary endpoints
#' @param adseq boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seed seed for the simulation
#'
#' @return A numeric matrix containing the simulated hypothesis test results. The first column represents the overall equivalence decision, where 1 indicates success and 0 indicates failure. The subsequent columns contain the hypothesis test results for each endpoint, followed by mean estimates for the reference and treatment groups, and standard deviations for the reference and treatment groups.
#'
#' @export
test_2x2_rom <- function(n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed) {
    .Call(`_SimTOST_test_2x2_rom`, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed)
}

#' @title Simulate a Parallel Design and Test Difference of Means (DOM)
#'
#' @description
#' Simulates a parallel-group design and performs equivalence testing using the difference of means (DOM) approach.
#' This function evaluates whether the treatment and reference groups are equivalent based on predefined
#' equivalence margins and hypothesis testing criteria.
#'
#' @param n integer number of subjects per arm
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaT matrix covar-variance matrix on treatment arm across endpoints
#' @param SigmaR matrix covar-variance matrix on reference arm across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param dropout vector of size 2 with dropout proportion per arm (T,R)
#' @param typey vector with positions of primary endpoints
#' @param adseq  boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seedT integer seed for the simulation on treatment arm
#' @param arm_seedR integer seed for the simulation on reference arm
#' @param TART double treatment allocation rate for the treatment arm
#' @param TARR double treatment allocation rate for the reference arm
#' @param vareq boolean assumed equivalence variance between arms for the t-test
#'
#' @return A numeric matrix containing the simulated hypothesis test results.
#' The first column represents the overall equivalence decision, where 1 indicates
#' success and 0 indicates failure. The subsequent columns contain the hypothesis
#' test results for each endpoint, followed by mean estimates for the reference and
#' treatment groups, and standard deviations for the reference and treatment groups.
#'
#' @details
#' The function simulates a parallel-group study design and evaluates equivalence
#' using the difference of means (DOM) approach. It accounts for dropout rates and
#' treatment allocation proportions while generating simulated data based on the
#' specified covariance structure. The test statistics are computed, and a final
#' equivalence decision is made based on the predefined number of required significant
#' endpoints (\code{k}). If sequential testing (\code{adseq}) is enabled, primary endpoints
#' must establish equivalence before secondary endpoints are evaluated.
#' When \code{vareq = TRUE}, the test assumes equal variances between groups and
#' applies Schuirmann's two one-sided tests (TOST).
#'
#' @export
test_par_dom <- function(n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq) {
    .Call(`_SimTOST_test_par_dom`, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq)
}

#' @title Simulate a Parallel Design and Test Ratio of Means (ROM)
#'
#' @description
#' Simulates a parallel-group design and performs equivalence testing using the ratio of means (ROM) approach.
#' This function evaluates whether the treatment and reference groups are equivalent based on predefined
#' equivalence margins and hypothesis testing criteria.
#'
#' @param n integer number of subjects per arm
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaT matrix covar-variance matrix on treatment arm across endpoints
#' @param SigmaR matrix covar-variance matrix on reference arm across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param dropout vector of size 2 with dropout proportion per arm (T,R)
#' @param typey vector with positions of primary endpoints
#' @param adseq boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seedT integer seed for the simulation on treatment arm
#' @param arm_seedR integer seed for the simulation on reference arm
#' @param TART double treatment allocation rate for the treatment arm
#' @param TARR double treatment allocation rate for the reference arm
#' @param vareq Boolean. If \code{TRUE}, assumes equal variance between arms and applies Schuirmann's two one-sided tests (TOST) for equivalence using a pooled variance.
#'
#' @return A numeric matrix containing the simulated hypothesis test results.
#' The first column represents the overall equivalence decision, where 1 indicates
#' success and 0 indicates failure. The subsequent columns contain the hypothesis
#' test results for each endpoint, followed by mean estimates for the reference and
#' treatment groups, and standard deviations for the reference and treatment groups.
#'
#' @details
#' The function simulates a parallel-group study design and evaluates equivalence
#' using the ratio of means (ROM) approach. It accounts for dropout rates and
#' treatment allocation proportions while generating simulated data based on the
#' specified covariance structure. The test statistics are computed, and a final
#' equivalence decision is made based on the predefined number of required significant
#' endpoints (\code{k}). If sequential testing (\code{adseq}) is enabled, primary endpoints
#' must establish equivalence before secondary endpoints are evaluated.
#' When \code{vareq = TRUE}, the test assumes equal variances between groups and
#' applies Schuirmann's two one-sided tests (TOST).
#'
#' @export
test_par_rom <- function(n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq) {
    .Call(`_SimTOST_test_par_rom`, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq)
}

#' @title Run Simulations for a Parallel Design with Difference of Means (DOM) test
#'
#' @description
#' This function simulates a parallel-group trial across multiple iterations.
#' It evaluates equivalence across multiple endpoints using the
#' Difference of Means (DOM) test.
#'
#' @param nsim Integer. The number of simulations to run.
#' @param n Integer. The sample size per arm (before dropout).
#' @param muT arma::vec. Mean vector for the treatment arm.
#' @param muR arma::vec. Mean vector for the reference arm.
#' @param SigmaT arma::mat. Covariance matrix for the treatment arm.
#' @param SigmaR arma::mat. Covariance matrix for the reference arm.
#' @param lequi_tol arma::rowvec. Lower equivalence thresholds for each endpoint.
#' @param uequi_tol arma::rowvec. Upper equivalence thresholds for each endpoint.
#' @param alpha arma::rowvec. Significance level for each endpoint.
#' @param dropout arma::vec. Dropout rates for each arm (T, R).
#' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
#' @param adseq Boolean. If \code{TRUE}, applies sequential (hierarchical) testing.
#' @param k Integer. Minimum number of endpoints required for equivalence.
#' @param arm_seed_T arma::ivec. Random seed vector for the treatment group (one per simulation).
#' @param arm_seed_R arma::ivec. Random seed vector for the reference group (one per simulation).
#' @param TART Double. Treatment allocation ratio (proportion of subjects in treatment arm).
#' @param TARR Double. Reference allocation ratio (proportion of subjects in reference arm).
#' @param vareq Boolean. If \code{TRUE}, assumes equal variances across treatment and reference groups.
#'
#' @details
#' Equivalence testing uses either the Difference of Means (DOM) test,
#' applying predefined equivalence thresholds and significance levels. When hierarchical testing (\code{adseq})
#' is enabled, all primary endpoints must demonstrate equivalence before secondary endpoints are evaluated.
#' Dropout rates are incorporated into the sample size calculation to ensure proper adjustment.
#' Randomization is controlled through separate random seeds for the treatment and reference groups,
#' enhancing reproducibility.
#'
#' @return
#' The function returns an arma::mat storing simulation results row-wise for consistency
#' with R's output format. The first row (\code{totaly}) contains the overall equivalence decision
#' (1 for success, 0 for failure). The subsequent rows include equivalence deicisons for each endpoint,
#' mean estimates for both treatment and reference groups, and corresponding standard deviations.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
run_simulations_par_dom <- function(nsim, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seed_T, arm_seed_R, TART, TARR, vareq) {
    .Call(`_SimTOST_run_simulations_par_dom`, nsim, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seed_T, arm_seed_R, TART, TARR, vareq)
}

#' @title Run Simulations for a Parallel Design with Ratio of Means (ROM) test
#'
#' @description
#' This function simulates a parallel-group trial across multiple iterations.
#' It evaluates equivalence across multiple endpoints using the
#' Ratio of Means (ROM) test.
#'
#' @param nsim Integer. The number of simulations to run.
#' @param n Integer. The sample size per arm (before dropout).
#' @param muT arma::vec. Mean vector for the treatment arm.
#' @param muR arma::vec. Mean vector for the reference arm.
#' @param SigmaT arma::mat. Covariance matrix for the treatment arm.
#' @param SigmaR arma::mat. Covariance matrix for the reference arm.
#' @param lequi_tol arma::rowvec. Lower equivalence thresholds for each endpoint.
#' @param uequi_tol arma::rowvec. Upper equivalence thresholds for each endpoint.
#' @param alpha arma::rowvec. Significance level for each endpoint.
#' @param dropout arma::vec. Dropout rates for each arm (T, R).
#' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
#' @param adseq Boolean. If \code{TRUE}, applies sequential (hierarchical) testing.
#' @param k Integer. Minimum number of endpoints required for equivalence.
#' @param arm_seed_T arma::ivec. Random seed vector for the treatment group (one per simulation).
#' @param arm_seed_R arma::ivec. Random seed vector for the reference group (one per simulation).
#' @param TART Double. Treatment allocation ratio (proportion of subjects in treatment arm).
#' @param TARR Double. Reference allocation ratio (proportion of subjects in reference arm).
#' @param vareq Boolean. If \code{TRUE}, assumes equal variances across treatment and reference groups.
#'
#' @details
#' Equivalence testing uses either the Ratio of Means (ROM) test,
#' applying predefined equivalence thresholds and significance levels. When hierarchical testing (\code{adseq})
#' is enabled, all primary endpoints must demonstrate equivalence before secondary endpoints are evaluated.
#' Dropout rates are incorporated into the sample size calculation to ensure proper adjustment.
#' Randomization is controlled through separate random seeds for the treatment and reference groups,
#' enhancing reproducibility.
#'
#' @return
#' The function returns an arma::mat storing simulation results row-wise for consistency
#' with R's output format. The first row (\code{totaly}) contains the overall equivalence decision
#' (1 for success, 0 for failure). The subsequent rows include equivalence decisions for each endpoint,
#' mean estimates for both treatment and reference groups, and corresponding standard deviations.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
run_simulations_par_rom <- function(nsim, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seed_T, arm_seed_R, TART, TARR, vareq) {
    .Call(`_SimTOST_run_simulations_par_rom`, nsim, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seed_T, arm_seed_R, TART, TARR, vareq)
}

#' @title Run Simulations for a 2x2 Crossover Design with Difference of Means (DOM) test
#'
#' @description
#' This function simulates a 2x2 crossover trial across multiple iterations.
#' It evaluates equivalence across multiple endpoints using the
#' Difference of Means (DOM) test.
#'
#' @param nsim Integer. The number of simulations to run.
#' @param n Integer. The sample size per period.
#' @param muT Numeric vector. Mean outcomes for the active treatment.
#' @param muR Numeric vector. Mean outcomes for the reference treatment.
#' @param SigmaW Numeric matrix. Within-subject covariance matrix for endpoints.
#' @param lequi_tol Numeric vector. Lower equivalence thresholds for each endpoint.
#' @param uequi_tol Numeric vector. Upper equivalence thresholds for each endpoint.
#' @param alpha Numeric vector. Significance levels for hypothesis testing across endpoints.
#' @param sigmaB Numeric. Between-subject variance for the crossover model.
#' @param dropout Numeric vector of size 2. Dropout rates for each sequence.
#' @param Eper Numeric vector. Expected period effects for each sequence.
#' @param Eco Numeric vector. Expected carryover effects for each sequence.
#' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
#' @param adseq Logical. If \code{TRUE}, applies sequential (hierarchical) testing.
#' @param k Integer. Minimum number of endpoints required for equivalence.
#' @param arm_seed Integer vector. Random seed for each simulation.
#'
#' @details
#' This function evaluates equivalence using the Difference of Means (DOM) test.
#' Equivalence is determined based on predefined lower (\code{lequi_tol}) and upper (\code{uequi_tol}) equivalence thresholds,
#' and hypothesis testing is conducted at the specified significance level (\code{alpha}).
#' If \code{adseq} is \code{TRUE}, primary endpoints must establish equivalence before secondary endpoints are evaluated.
#' The sample size per period is adjusted based on dropout rates, ensuring valid study conclusions.
#' The simulation incorporates within-subject correlation using \code{SigmaW} and accounts for between-subject variance with \code{sigmaB}.
#' Expected period effects (\code{Eper}) and carryover effects (\code{Eco}) are included in the model.
#' A fixed random seed (\code{arm_seed}) is used to ensure reproducibility across simulations.
#'
#' @return
#' A numeric matrix where each column stores simulation results:
#' The first row (\code{totaly}) represents the overall equivalence decision (1 = success, 0 = failure).
#' Subsequent rows contain equivalence decisions per endpoint,
#' mean estimates for the treatment group, mean estimates for the reference group,
#' standard deviations for treatment, and standard deviations for reference.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
run_simulations_2x2_dom <- function(nsim, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed) {
    .Call(`_SimTOST_run_simulations_2x2_dom`, nsim, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed)
}

#' @title Run Simulations for a 2x2 Crossover Design with Ratio of Means (ROM) test
#'
#' @description
#' This function simulates a 2x2 crossover trial across multiple iterations.
#' It evaluates equivalence across multiple endpoints using the
#' Ratio of Means (ROM) test.
#'
#' @param nsim Integer. The number of simulations to run.
#' @param n Integer. The sample size per period.
#' @param muT Numeric vector. Mean outcomes for the active treatment.
#' @param muR Numeric vector. Mean outcomes for the reference treatment.
#' @param SigmaW Numeric matrix. Within-subject covariance matrix for endpoints.
#' @param lequi_tol Numeric vector. Lower equivalence thresholds for each endpoint.
#' @param uequi_tol Numeric vector. Upper equivalence thresholds for each endpoint.
#' @param alpha Numeric vector. Significance levels for hypothesis testing across endpoints.
#' @param sigmaB Numeric. Between-subject variance for the crossover model.
#' @param dropout Numeric vector of size 2. Dropout rates for each sequence.
#' @param Eper Numeric vector. Expected period effects for each sequence.
#' @param Eco Numeric vector. Expected carryover effects for each sequence.
#' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
#' @param adseq Logical. If \code{TRUE}, applies sequential (hierarchical) testing.
#' @param k Integer. Minimum number of endpoints required for equivalence.
#' @param arm_seed Integer vector. Random seed for each simulation.
#'
#' @details
#' This function evaluates equivalence using the Ratio of Means (ROM) test.
#' Equivalence is determined based on predefined lower \code{lequi_tol} and upper \code{uequi_tol} equivalence thresholds,
#' and hypothesis testing is conducted at the specified significance level \code{alpha}.
#' If \code{adseq} is \code{TRUE}, primary endpoints must establish equivalence before secondary endpoints are evaluated.
#' The sample size per period is adjusted based on dropout rates, ensuring valid study conclusions.
#' The simulation incorporates within-subject correlation using \code{SigmaW} and accounts for between-subject variance with \code{sigmaB}.
#' Expected period effects \code{Eper} and carryover effects \code{Eco} are included in the model.
#' A fixed random seed \code{arm_seed} is used to ensure reproducibility across simulations.//'
#'
#' @return
#' A numeric matrix where each column stores simulation results:
#' The first row (\code{totaly}) represents the overall equivalence decision (1 = success, 0 = failure).
#' Subsequent rows contain equivalence decisions per endpoint,
#' mean estimates for the treatment group, mean estimates for the reference group,
#' standard deviations for treatment, and standard deviations for reference.
#'
#'  @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
run_simulations_2x2_rom <- function(nsim, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed) {
    .Call(`_SimTOST_run_simulations_2x2_rom`, nsim, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed)
}

