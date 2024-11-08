#' @title Simulated Test Statistic for Given Settings in Noninferiority/Equivalence Trials
#'
#' @description This function simulates test statistics for multiple hypothesis testing in the context of biosimilar development, using the approach described by Mielke et al. (2018). It calculates the necessary sample size for meeting equivalence criteria across multiple endpoints, considering correlation structures, and applying multiplicity adjustments.
#'
#' @details Following the methodology by Mielke et al. (2018), this function is designed for multiple endpoint clinical trials where success is defined as meeting equivalence criteria on at least a subset of tests. The function simulates test statistics based on multivariate normal distribution assumptions and supports k-out-of-m success criteria for regulatory approval. Additionally, Type I error control is achieved through multiplicity adjustments as proposed by Lehmann and Romano (2005) to ensure rigorous error rate management.
#'
#' The approach is particularly relevant for biosimilar studies where sample size estimation must account for multiple comparisons across endpoints, doses, or populations, as regulatory agencies often require equivalence across all relevant comparisons. This method provides a framework for estimating power and sample size, even when equivalence is not required on all endpoints.
#'
#' @param N Integer. The number of subjects per sequence.
#' @param m Integer. The number of endpoints.
#' @param k Integer. The number of endpoints that must be successful to consider the test a success.
#' @param R Matrix. The correlation matrix between the endpoints (e.g., generated with the `variance.const.corr()` function), with dimensions m x m.
#' @param sigma Numeric. The standard deviation of the endpoints. This can either be a vector of length m (one for each endpoint) or a single value. If a single value is provided, it is assumed that the standard deviation is constant across all endpoints. In the case of a 2x2 crossover design, the input is the within-subject variance; for a parallel group design, it represents the standard deviation in the treatment group, assumed identical for both test and reference groups.
#' @param true.diff Numeric. The assumed true difference between the test and reference. This can be a vector of length m (one for each endpoint) or a single value.
#' @param equi.tol Numeric. Equivalence margins, with the interval being (-equi.tol, +equi.tol).
#' @param design Character. The study design, either "22co" for a 2x2 crossover design or "parallel" for a parallel groups design.
#' @param alpha Numeric. The significance level.
#' @param adjust Character. The method for multiplicity adjustment. Options are "no" for no adjustment, "bon" for Bonferroni correction, or "k" for k-adjustment.
#'
#' @references
#'
#' Kong, L., Kohberger, R. C. & Koch, G. G. Type I Error and Power in Noninferiority/Equivalence Trials with Correlated Multiple Endpoints: An Example from Vaccine Development Trials. Journal of Biopharmaceutical Statistics 14, 893–907 (2004).
#'
#' Lehmann, E. L. & Romano, J. P. Generalizations of the Familywise Error Rate. The Annals of Statistics 33, 1138–1154 (2005).
#'
#' Mielke, J., Jones, B., Jilma, B. & König, F. Sample Size for Multiple Hypothesis Testing in Biosimilar Development. Statistics in Biopharmaceutical Research 10, 39–49 (2018).
#'
#' @return A realization of the simulated test statistic for the given setting. If values are provided as a single number (constant for all endpoints), a vector is created.
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

#' @title Power Calculation for Hypothesis Testing Using Mielke's Method
#' @description Calculate power of hypothesis testing using Mielke's
#'
#' @param N Integer. The number of subjects per sequence.
#' @param m Integer. The number of endpoints.
#' @param k Integer. The number of endpoints that must be successful to consider the test a success.
#' @param R Matrix. The correlation matrix between the endpoints (e.g., generated with the `variance.const.corr()` function), with dimensions m x m.
#' @param sigma Numeric. The standard deviation of the endpoints. This can either be a vector of length m (one for each endpoint) or a single value. If a single value is provided, it is assumed that the standard deviation is constant across all endpoints. In the case of a 2x2 crossover design, the input is the within-subject variance; for a parallel group design, it represents the standard deviation in the treatment group, assumed identical for both test and reference groups.
#' @param true.diff Numeric. The assumed true difference between the test and reference. This can be a vector of length m (one for each endpoint) or a single value.
#' @param equi.tol Numeric. Equivalence margins, with the interval being (-equi.tol, +equi.tol). Default is log(1.25).
#' @param design Character. The study design, either "22co" for a 2x2 crossover design or "parallel" for a parallel groups design.
#' @param alpha Numeric. The significance level. Default is 0.05.
#' @param adjust Character. The method for multiplicity adjustment. Options are "no" for no adjustment, "bon" for Bonferroni correction, or "k" for k-adjustment.
#' @param nsim Integer. The number of simulations to perform. Default is 10,000.
#'
#' @return Numeric. The estimated power based on the simulations.
#'
#' @keywords internal
#'
power_Mielke <- function(N, m, k, R, sigma, true.diff, equi.tol = log(1.25),
                         design, alpha=0.05, adjust="no", nsim = 10000) {
  mean(replicate(nsim, sign_Mielke(N = N, m = m, k = k, R = R, sigma = sigma,
                                   true.diff = true.diff, equi.tol = equi.tol,
                                   design = design, alpha = alpha,
                                   adjust = adjust)))
}

#' @title Power Calculation for Hypothesis Testing of Difference of Means (DOM)
#'
#' @description This function calculates the power of hypothesis testing for the difference of means (DOM) between two groups, using simulations to estimate the achieved power based on the provided parameters such as sample sizes, means, standard deviations, and significance level.
#'
#' @param seed Integer. A seed value for reproducibility of the simulations.
#' @param mu_test Numeric. The arithmetic mean of the test group.
#' @param mu_control Numeric. The arithmetic mean of the control group.
#' @param sigma_test Numeric. The standard deviation of the test group.
#' @param sigma_control Numeric. The standard deviation of the control group.
#' @param N_test Integer. The sample size of the test group.
#' @param N_control Integer. The sample size of the control group.
#' @param alpha Numeric. The significance level for the hypothesis test. Default is 0.05.
#' @param nsim Integer. The number of simulations to perform. Default is 10,000.
#'
#' @return Numeric. The estimated power based on the simulations.
#'
#' @keywords internal
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
#' @export
N_Mielke <- function(power, Nmax, m, k, rho, sigma, true.diff, equi.tol,
                     design, alpha, adjust="no", seed = NULL, nsim = 10000) {

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
