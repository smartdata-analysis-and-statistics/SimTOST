#' @title Simulated Test Statistic for Noninferiority/Equivalence Trials
#'
#' @description
#' Simulates test statistics for multiple hypothesis testing in biosimilar development,
#' following the approach described by Mielke et al. (2018). It calculates the necessary
#' sample size for meeting equivalence criteria across multiple endpoints while
#' considering correlation structures and applying multiplicity adjustments.
#'
#' @details
#' This function is designed for multiple-endpoint clinical trials, where success
#' is defined as meeting equivalence criteria for at least a subset of tests.
#' Simulated test statistics are based on multivariate normal distribution assumptions,
#' and the function supports k-out-of-m success criteria for regulatory approval.
#'
#' Type I error control is achieved through multiplicity adjustments as proposed by
#' Lehmann and Romano (2005) to ensure rigorous error rate management. This approach
#' is particularly relevant for biosimilar studies, where sample size estimation must
#' account for multiple comparisons across endpoints, doses, or populations.
#'
#' @param N Integer specifying the number of subjects per sequence.
#' @param m Integer specifying the number of endpoints.
#' @param k Integer specifying the number of endpoints that must meet equivalence
#' to consider the test successful.
#' @param R Matrix specifying the correlation structure between endpoints.
#' This should be an \code{m x m} matrix, e.g., generated using \code{variance.const.corr()}.
#' @param sigma Numeric specifying the standard deviation of endpoints.
#' Can be a vector of length \code{m} (one per endpoint) or a single value.
#' In a 2x2 crossover design, this represents within-subject variance.
#' In a parallel-group design, it represents the treatment group standard deviation.
#' @param true.diff Numeric specifying the assumed true difference between test and reference.
#' Can be a vector of length \code{m} or a single value.
#' @param equi.tol Numeric specifying the equivalence margins.
#' The interval is defined as \code{(-equi.tol, +equi.tol)}.
#' @param design Character specifying the study design.
#' Options are \code{"22co"} for a 2x2 crossover design or \code{"parallel"} for a parallel-group design.
#' @param alpha Numeric specifying the significance level.
#' @param adjust Character specifying the method for multiplicity adjustment.
#' Options include \code{"no"} for no adjustment, \code{"bon"} for Bonferroni correction,
#' and \code{"k"} for k-adjustment.
#'
#' @references
#' Kong, L., Kohberger, R. C., & Koch, G. G. (2004). Type I Error and Power in
#' Noninferiority/Equivalence Trials with Correlated Multiple Endpoints: An Example
#' from Vaccine Development Trials. \emph{Journal of Biopharmaceutical Statistics, 14}(4), 893–907.
#'
#' Lehmann, E. L., & Romano, J. P. (2005). Generalizations of the Familywise Error Rate.
#' \emph{The Annals of Statistics, 33}(2), 1138–1154.
#'
#' Mielke, J., Jones, B., Jilma, B., & König, F. (2018). Sample Size for Multiple
#' Hypothesis Testing in Biosimilar Development. \emph{Statistics in Biopharmaceutical Research, 10}(1), 39–49.
#'
#' @return
#' A numeric vector representing a realization of the simulated test statistic for the given setting.
#'
#'
#' @keywords internal
sign_Mielke <- function(N, m, k, R, sigma, true.diff, equi.tol = log(1.25),
                        design, alpha = 0.05, adjust = "no") {
  if (length(true.diff) == 1) {
    true.diff <- true.diff * rep(1, m)
  }
  if (length(sigma) == 1) {
    sigma <- sigma * rep(1, m)
  }
  if (length(equi.tol) == 1) {
    equi.tol <- equi.tol * rep(1, m)
  }
  # generate data and calculate Cholesky decomposition
  test.raw <- stats::rnorm(m)
  T.matrix <- chol(R)

  # adjustment for type 1-error rate
  if (adjust == "k") {
    alpha <- k*alpha/m
  } else if (adjust == "bon") {
    alpha <- alpha/m
  } else if (adjust == "no") {
    alpha <- alpha
  } # calculate the test statistics (as stated in the paper)
  if (design == "parallel") {
    Z1 <- t(T.matrix) %*% test.raw + sqrt(N/2)/sigma * (equi.tol - true.diff)
    Z2 <- -t(T.matrix) %*% test.raw + sqrt(N/2)/sigma * (equi.tol + true.diff)
  }
  if (design == "22co") {
    Z1 <- t(T.matrix) %*% test.raw + sqrt(2*N)/sigma * (equi.tol - true.diff)
    Z2 <- -t(T.matrix) %*% test.raw + sqrt(2*N)/sigma * (equi.tol + true.diff)
  }
  min.Z <- apply(cbind(Z1,Z2),1,min)
  #decide for each endpoints
  # if Z1 or Z2 is smaller
  critic <- stats::qnorm(1 - alpha)
  # decide how many endpoints are successful and if
  # that number is higher than the requested number k reject H0
  test.dec <- (sum(min.Z > critic) >= k)
  return(dec = test.dec)
}

#' @title Power Calculation for Hypothesis Testing in Equivalence Trials
#'
#' @description
#' Estimates the power of hypothesis testing in equivalence trials using
#' the method described by Mielke et al. This approach accounts for multiple
#' endpoints, correlation structures, and multiplicity adjustments.
#'
#' @param N Integer specifying the number of subjects per sequence.
#' @param m Integer specifying the number of endpoints.
#' @param k Integer specifying the number of endpoints that must meet equivalence
#' to consider the test successful.
#' @param R Matrix specifying the correlation structure between endpoints.
#' This should be an \code{m x m} matrix, e.g., generated using \code{variance.const.corr()}.
#' @param sigma Numeric specifying the standard deviation of endpoints.
#' Can be a vector of length \code{m} (one per endpoint) or a single value.
#' In a 2x2 crossover design, this represents within-subject variance.
#' In a parallel-group design, it represents the treatment group standard deviation.
#' @param true.diff Numeric specifying the assumed true difference between test and reference.
#' Can be a vector of length \code{m} or a single value.
#' @param equi.tol Numeric specifying the equivalence margins, with the interval defined as
#' \code{(-equi.tol, +equi.tol)}. Default is \code{log(1.25)}.
#' @param design Character specifying the study design.
#' Options are \code{"22co"} for a 2x2 crossover design or \code{"parallel"} for a parallel-group design.
#' @param alpha Numeric specifying the significance level. Default is \code{0.05}.
#' @param adjust Character specifying the method for multiplicity adjustment.
#' Options include \code{"no"} for no adjustment, \code{"bon"} for Bonferroni correction,
#' and \code{"k"} for k-adjustment.
#' @param nsim Integer specifying the number of simulations to perform. Default is \code{10,000}.
#'
#' @return
#' A numeric value representing the estimated power based on the simulations.
#'
#' @keywords internal
power_Mielke <- function(N, m, k, R, sigma, true.diff, equi.tol = log(1.25),
                         design, alpha=0.05, adjust="no", nsim = 10000) {
  mean(replicate(nsim, sign_Mielke(N = N, m = m, k = k, R = R, sigma = sigma,
                                   true.diff = true.diff, equi.tol = equi.tol,
                                   design = design, alpha = alpha,
                                   adjust = adjust)))
}

#' @title Power Calculation for Difference of Means (DOM) Hypothesis Test
#'
#' @description Computes the statistical power for testing the difference of means (DOM) between two groups using Monte Carlo simulations. The power is estimated based on specified sample sizes, means, standard deviations, and significance level.
#'
#' @param seed Integer. Seed for reproducibility.
#' @param mu_test Numeric. Mean of the test group.
#' @param mu_control Numeric. Mean of the control group.
#' @param sigma_test Numeric. Standard deviation of the test group.
#' @param sigma_control Numeric. Standard deviation of the control group.
#' @param N_test Integer. Sample size of the test group.
#' @param N_control Integer. Sample size of the control group.
#' @param lb Numeric. Lower bound for the equivalence margin.
#' @param ub Numeric. Upper bound for the equivalence margin.
#' @param alpha Numeric. Significance level (default = 0.05).
#' @param nsim Integer. Number of simulations (default = 10,000).
#'
#' @return Numeric. Estimated power (probability between 0 and 1).
power_dom <- function(seed, mu_test, mu_control, sigma_test, sigma_control,
                      N_test, N_control, lb, ub, alpha = 0.05, nsim = 10000) {
  set.seed(seed)

  sign <- rep(NA, nsim)
  for (sim in 1:nsim) {
    x_T <- stats::rnorm(n = N_test, mean = mu_test, sd = sigma_test)
    x_R <- stats::rnorm(n = N_control, mean = mu_control, sd = sigma_control)

    xbar_T <- mean(x_T)
    xbar_R <- mean(x_R)

    df <- N_test + N_control - 2
    S <- sqrt(((N_test - 1)*stats::var(x_T) + (N_control - 1)*stats::var(x_R))/df)
    M <- 1/((1/N_test) + (1/N_control))

    # Calculate the test statistics
    tlb <- ((xbar_T - xbar_R) - lb)/(S/sqrt(M))
    tub <- (ub - (xbar_T - xbar_R))/(S/sqrt(M))
    tref <- stats::qt(1 - alpha, df = df)

    # Determine whether to reject H0
    sign[sim] <- (tlb >= tref & tub >= tref)
  }

  # Power
  sum(sign)/nsim
}

#' @title Sample Size Estimation for Multiple Hypothesis Testing Using Mielke's Method
#'
#' @description Estimates the required sample size to achieve a specified power level for multiple hypothesis testing, using the approach described by Mielke et al. (2018). This function is particularly useful for bioequivalence or biosimilar studies with multiple correlated endpoints, where a minimum number of endpoints must meet equivalence criteria.
#'
#' @references
#'
#' Mielke, J., Jones, B., Jilma, B. & König, F. Sample Size for Multiple Hypothesis Testing in Biosimilar Development. Statistics in Biopharmaceutical Research 10, 39–49 (2018).
#'
#'
#' @param power Numeric. Desired statistical power.
#' @param Nmax Integer. Maximum allowable sample size.
#' @param m Integer. Total number of endpoints.
#' @param k Integer. Number of endpoints that must meet the success criteria for overall study success.
#' @param rho Numeric. Constant correlation coefficient among endpoints.
#' @param sigma Numeric or vector. Standard deviation of each endpoint. If a single value is provided, it is assumed to be constant across all endpoints. In a 2x2 crossover design, this is the within-subject standard deviation; in a parallel design, it represents the treatment group’s standard deviation, assumed to be the same for both test and reference.
#' @param true.diff Numeric or vector. Assumed true difference between test and reference for each endpoint. If a single value is provided, it is applied uniformly across all endpoints.
#' @param equi.tol Numeric. Equivalence margin; the equivalence interval is defined as (-equi.tol, +equi.tol).
#' @param design Character. Study design, either "22co" for a 2x2 crossover design or "parallel" for a parallel groups design.
#' @param alpha Numeric. Significance level for the hypothesis test.
#' @param adjust Character. Method for multiplicity adjustment; options are "no" (no adjustment), "bon" (Bonferroni adjustment), and "k" (k-adjustment).
#' @param seed Integer. Random seed for reproducibility.
#' @param nsim Integer. Number of simulations to run for power estimation (default: 10,000).
#'
#' @details This function uses the method proposed by Mielke et al. (2018) to estimate the sample size required to achieve the desired power level in studies with multiple correlated endpoints. The function iteratively increases sample size until the target power is reached or the maximum allowable sample size (Nmax) is exceeded. The approach accounts for endpoint correlation and supports adjustments for multiple testing using various correction methods.
#'
#' @return A named vector containing:
#' \describe{
#'   \item{"power.a"}{Achieved power with the estimated sample size.}
#'   \item{"SS"}{Required sample size per sequence to achieve the target power.}
#' }
#'
#' @examples
#' # Example 1 from Mielke
#' sampleSize_Mielke(power = 0.8, Nmax = 1000, m = 5, k = 5, rho = 0,
#'                   sigma = 0.3, true.diff =  log(1.05), equi.tol = log(1.25),
#'                   design = "parallel", alpha = 0.05, adjust = "no",
#'                   seed = 1234, nsim = 100)
#'
#' @export
sampleSize_Mielke <- function(power, Nmax, m, k, rho, sigma, true.diff,
                              equi.tol, design, alpha, adjust = "no",
                              seed = NULL, nsim = 10000) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  R <- (1 - rho) * diag(m) + rho * rep(1,m) %*% t(rep(1,m))

  #generate covariance matrix
  ##!! for arbitrary correlation matrix replace R with the desired matrix !!
  for (Ni in k:Nmax) {
    # increase one by one until the desired power or the
    # maximal sample size is reached
    power.esti <- power_Mielke(N = Ni, m = m, k = k, R = R, sigma = sigma,
                               true.diff = true.diff, equi.tol = equi.tol,
                               design = design, alpha = alpha, adjust = adjust,
                               nsim = nsim)
    if (power.esti > power) {
      break
    }
  }
  return(c(power.a = power.esti, SS = Ni))
}
