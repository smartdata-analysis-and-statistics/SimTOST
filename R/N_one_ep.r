# teststatistic - generates the test statistics as described in the paper
## by Kong, Kohberger, Koch (2004)
## The type-1 error rate adjustment as proposed by
## Lehmann, Roman (2005)

#' @title Simulated Test Statistic for Given Settings Using Mielke's Method
#'
#' @description This function simulates a test statistic for a given setting using Mielke's method, accounting for multiple endpoints, correlation structure, and study design.
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
#' @return A realization of the simulated test statistic for the given setting. If values are provided as a single number (constant for all endpoints), a vector is created.
#'
#' @keywords internal
sign_Mielke <- function(N, m, k, R, sigma, true.diff, equi.tol=log(1.25), design, alpha=0.05, adjust="no") {
  if (length(true.diff)==1) {
    true.diff <- true.diff * rep(1,m)
  }
  if (length(sigma)==1) {
    sigma <- sigma * rep(1,m)
  }
  if (length(equi.tol)==1) {
    equi.tol <- equi.tol * rep(1,m)
  }
  # generate data and calculate cholesky decomposition
  test.raw <- stats::rnorm(m)
  T.matrix <- chol(R)
  # adjustment for type 1-error rate
  if (adjust=="k") {
    alpha <- k*alpha/m
  } else if (adjust=="bon") {
    alpha <- alpha/m
  } else if (adjust=="no") {
    alpha <- alpha
  } # calculate the test statistics (as stated in the paper)
  if (design=="parallel") {
    Z1 <- t(T.matrix) %*% test.raw+sqrt(N/2)/sigma * (equi.tol-true.diff)
    Z2 <- -t(T.matrix) %*% test.raw+sqrt(N/2)/sigma * (equi.tol+true.diff)
  }
  if (design=="22co") {
    Z1 <- t(T.matrix) %*% test.raw+sqrt(2*N)/sigma * (equi.tol-true.diff)
    Z2 <- -t(T.matrix) %*% test.raw+sqrt(2*N)/sigma * (equi.tol+true.diff)
  }
  min.Z <- apply(cbind(Z1,Z2),1,min)
  #decide for each endpoints
  # if Z1 or Z2 is smaller
  critic <- stats::qnorm(1 - alpha)
  # decide how many endpoints are successful and if
  # that number is higher than the requested number k reject H0
  test.dec <- (sum(min.Z > critic)>=k)
  return(dec=test.dec)
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

#' @title N_Mielke
#'
#' @description Sample size estimation for Multiple Hypothesis Testing (Mielke)
#'
#' @param power Desired power
#' @param Nmax Maximal sample size
#' @param m Number of endpoints (integer)
#' @param k Number of endpoints that must be successful (integer)
#' @param rho Correlation (this program assumes a constant correlation between the endpoints)
#' @param sigma Standard deviation of the endpoints (can be a vector of length m or a single value. If it is a single value, it is assumed that the standard deviation is constant for all endpoints). In case of the 2x2 crossover design, it is assumed that the input is the within-subject variance, in case of a parallel groups design, it is the standard deviation of the endpoint in the treatment group and it is assumed that this is identical for test and reference
#' @param true.diff Assumed true difference between test and reference (can be a vector of length m or a single value, see above)
#' @param equi.tol Equivalence margins, positive number, the interval is (-equi.tol, +equi.tol)
#' @param design Study design ("22co" for a 2x2 crossover design and "parallel" for a parallel groups design)
#' @param alpha Significance level
#' @param adjust Should an adjustment for multiplicity be applied ("no": no adjustment. "bon": Bonferroni, "k": k-adjustment)
#' @param seed Random Seed
#' @param nsim Number of iterations for power calculations (default: 10000)
#'
#' @return A vector with the achieved power (first entry) and the required sample size per sequence (second entry)
#' @export
N_Mielke <- function(power, Nmax, m, k, rho, sigma, true.diff, equi.tol, design, alpha, adjust="no", seed = NULL, nsim = 10000) {

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
