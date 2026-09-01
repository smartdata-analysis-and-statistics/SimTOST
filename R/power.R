#' Estimate power at a fixed sample size
#'
#' @description Calculates simulated power for a prespecified sample size.
#' The outcome family is selected through `distribution`, matching the
#' unified `sampleSize()` interface.
#'
#' @param n Integer sample size or vector of sample sizes used for the
#' simulation. For parallel studies this is the base sample size used to
#' derive arm sizes; for 2x2 studies it is the number per sequence. A vector
#' returns a simpower_curve object and enables plotting power across base
#' sample sizes.
#' @param distribution Outcome distribution using R's names: `"norm"`,
#' `"lnorm"`, `"pois"`, or `"nbinom"` (case-insensitive). Longer labels
#' such as `"normal"`, `"lognormal"`, and `"poisson"` are also accepted.
#' @param mu_list Named list of continuous-outcome means per arm.
#' @param varcov_list Optional list of covariance matrices for continuous
#' outcomes.
#' @param sigma_list Optional list of standard-deviation vectors for
#' continuous outcomes.
#' @param cor_mat Optional endpoint correlation matrix. For count outcomes,
#' this is also the endpoint correlation matrix used by the joint count engine.
#' @param sigmaB Between-subject parameter for a 2x2 design.
#' @param rate_list Named arm-rate list for count outcomes.
#' @param exposure Count exposure per subject. This can be a scalar or
#' endpoint vector shared by arms, or a named list of arm-specific values.
#' @param dispersion Negative-binomial dispersion. This can be a scalar or
#' endpoint vector shared by arms, or a named list of arm-specific values.
#' @param Eper Period effects for a 2x2 design.
#' @param Eco Carry-over effects for a 2x2 design.
#' @param rho Common endpoint correlation when `cor_mat` is not supplied.
#' @param TAR Treatment allocation rates for continuous parallel designs.
#' @param arm_names Optional arm names.
#' @param ynames_list Optional endpoint names by arm.
#' @param type_y Endpoint hierarchy for sequential testing for continuous and
#' count outcomes. Use `1` for primary/co-primary and `2` for secondary
#' endpoints.
#' @param list_comparator Named list of treatment-reference comparisons. Each
#'   element must be `c(test, reference)`; the first arm is the test arm and
#'   the second arm is the reference arm.
#' @param list_y_comparator Endpoint selections by comparison.
#'   For count outcomes, the selected endpoints are used to define the count
#'   multiplicity and the effective `k`; in a joint count analysis all
#'   comparisons must currently use the same selected endpoint set.
#' @param alpha One-sided significance level.
#' @param lequi.tol Common lower equivalence bound.
#' @param uequi.tol Common upper equivalence bound.
#' @param list_lequi.tol Comparator-specific lower bounds.
#' @param list_uequi.tol Comparator-specific upper bounds.
#' @param dtype Trial design: `"parallel"` or `"2x2"`.
#' @param ctype Test type. Use `"DOM"` (difference of means) or `"ROM"`
#'   (ratio of means) for Normal and Lognormal outcomes, and `"RR"`
#'   (event-rate ratio) for Poisson and Negative Binomial outcomes.
#' @param vareq Whether variances are assumed equal for continuous outcomes.
#' @param k Number of endpoints required per comparison.
#' @param adjust Multiplicity adjustment: `"no"` (none), `"bon"`
#' (Bonferroni across all selected endpoints), `"sid"` (Sidak), `"k"`
#' (the existing weak K-fold rule), `"t"` (Mielke's strong
#' k-out-of-m adjustment, using `alpha / (m - k + 1)`; legacy `"pc"`
#' aliases are accepted), or
#' `"seq"` (sequential hierarchy).
#' @param dropout Dropout proportions. For count 2x2 studies, supply two
#' sequence-specific values.
#' @param nsim Number of simulated trials.
#' @param keep_sim_data Logical. If `TRUE`, retain model-scale observations
#'   for every simulated trial in `sim_data` for distribution diagnostics.
#'   Defaults to `FALSE` because retained data can be large.
#' @param seed Random seed.
#' @param ncores Number of computation cores. For continuous outcomes this is
#' passed to the compiled simulation backend; for count outcomes it splits
#' Monte Carlo trials into independent seeded chunks whose C++ results are
#' combined.
#' @return An object containing estimated power and its 95% Monte Carlo
#' confidence interval. The unified function returns primary class `simpower`
#' for every distribution. If `keep_sim_data = TRUE`, the object also contains
#' long-format model-scale observations in `sim_data`. Count results additionally inherit from the
#' compatibility class `countpower`.
#' @details
#' For Normal and Log Normal outcomes, supply the continuous-outcome inputs
#' (code{mu_list}, code{sigma_list} or code{varcov_list}, and optionally
#' code{cor_mat} or code{rho}). For Poisson and Negative Binomial outcomes,
#' supply code{rate_list}, code{list_comparator},
#' code{list_lequi.tol}, and code{list_uequi.tol}. Count code{exposure} and
#' code{dispersion} may be scalar, endpoint-specific, or named arm-specific
#' lists. Arguments for the other outcome family are not used.
#' The returned object supports code{summary()}, code{confint()}, and
#' code{plot()}.
#' The effective endpoint count is determined separately for each comparator
#' from `list_y_comparator` (or from the endpoints common to both arms when the
#' argument is omitted). `k` is checked against that count. The function warns
#' when an endpoint-wise adjustment is unnecessary because all selected
#' endpoints are required, and when `adjust = "no"` is used for a `k < m`
#' decision. `type_y` is used with `adjust = "seq"` for both continuous and
#' count outcomes; for other adjustments it is ignored with a warning.
#' For a `k`-of-`m` decision, `adjust = "t"` allocates alpha over the
#' `m - k + 1` boundary endpoints relevant to the strong k-out-of-m null.
#' The legacy `adjust = "pc"` label is accepted as an alias.
#' Comparisons always use the order supplied in `list_comparator`: the first
#' arm is the test and the second arm is the reference. Thus, `c("T", "R")`
#' gives `T - R` for DOM, `T / R` for ROM, and the event-rate ratio
#' `rate_T / rate_R` for RR. Reversing the two names reverses the estimand.
#' For count outcomes, the estimand is the event-rate ratio
#' \\(λ_T / λ_R\\). The equivalence hypotheses are
#' \\(H_0: λ_T/λ_R \\le L \\text{ or } λ_T/λ_R \\ge U\\)
#' versus \\(H_1: L < λ_T/λ_R < U\\), assessed by TOST. The same
#' interval hypotheses apply to continuous mean differences (DOM) or mean
#' ratios (ROM), with the log-normal ROM analysis performed on the log scale.
#' @export
#' @examples
#' simPower(n = 100, distribution = "Poisson",
#'       rate_list = list(TEST = .21, REF = .20),
#'       list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
#'       list_lequi.tol = list(TEST_vs_REF = .80),
#'       list_uequi.tol = list(TEST_vs_REF = 1.25),
#'       exposure = 10, nsim = 100, seed = 1)
simPower <- function(
    n, distribution = c("norm", "lnorm", "pois", "nbinom"),
    mu_list = NULL, varcov_list = NA, sigma_list = NA, cor_mat = NA,
    sigmaB = 0, rate_list = NULL, exposure = 1, dispersion = 0.1,
    Eper = c(0, 0), Eco = c(0, 0), rho = 0,
    TAR = NULL, arm_names = NA, ynames_list = NA, type_y = NA,
    list_comparator = NA, list_y_comparator = NA, alpha = 0.05,
    lequi.tol = NA, uequi.tol = NA, list_lequi.tol = NA,
    list_uequi.tol = NA, dtype = "parallel", ctype = "ROM",
    vareq = TRUE, k = NA, adjust = "no", dropout = NA, nsim = 5000,
    seed = 1234, ncores = 1, keep_sim_data = FALSE,
    .warn_redundant_bon = TRUE) {

  keep_sim_data <- .check_keep_sim_data(keep_sim_data)
  .validate_common_controls(NULL, alpha, nsim, ncores, dropout, rho)

  if (!is.numeric(n) || length(n) < 1L || any(!is.finite(n)) ||
      any(n != as.integer(n)) || any(n < 2L))
    stop("'n' must contain integer sample sizes >= 2.")

  if (length(n) > 1L) {
    original_call <- match.call()
    curve_results <- lapply(seq_along(n), function(i) {
      current_call <- original_call
      current_call$n <- n[[i]]
      eval(current_call, parent.frame())
    })
    curve <- data.frame(
      n = as.integer(n),
      n_total = vapply(curve_results, function(z) {
        if (!is.null(z$n_total)) return(as.numeric(z$n_total))
        if (!is.null(z$result) && !is.null(z$result$response$n_total)) {
          return(as.numeric(z$result$response$n_total))
        }
        as.numeric(z$n)
      }, numeric(1)),
      power = vapply(curve_results, function(z) z$power, numeric(1)),
      power_LCI = vapply(curve_results, function(z) z$power_LCI, numeric(1)),
      power_UCI = vapply(curve_results, function(z) z$power_UCI, numeric(1))
    )
    result <- curve_results[[1L]]
    result$n <- as.integer(n)
    if (!is.null(result$parameters)) result$parameters$n <- as.integer(n)
    result$power_curve <- curve
    result$curve_results <- curve_results
    class(result) <- c("simpower_curve", class(result))
    return(result)
  }
  if (missing(distribution) || is.null(distribution))
    distribution <- "lnorm"
  distribution <- .normalize_distribution(distribution)

  is_count <- distribution %in% c("pois", "nbinom")
  adjust <- .normalize_adjustment(adjust, allow_sequential = TRUE)
  type_y_supplied <- !(is.null(type_y) ||
                       (length(type_y) == 1L && is.na(type_y)))
  if (is_count) {
    if (missing(ctype)) {
      ctype <- "RR"
    } else if (length(ctype) != 1L || !is.character(ctype) ||
               is.na(ctype) || !identical(toupper(ctype), "RR")) {
      warning("For count outcomes, ctype = 'RR' is used. The supplied ctype is not available for count outcomes.")
      ctype <- "RR"
    } else {
      ctype <- "RR"
    }
  } else if (length(ctype) != 1L || !is.character(ctype) ||
             is.na(ctype) || !toupper(ctype) %in% c("DOM", "ROM")) {
    warning("For continuous outcomes, ctype must be 'DOM' or 'ROM'. Default ctype = 'ROM' is used.")
    ctype <- "ROM"
  } else {
    ctype <- toupper(ctype)
  }
  update_parameters <- list(
    n = as.integer(n), distribution = distribution,
    mu_list = mu_list, varcov_list = varcov_list, sigma_list = sigma_list,
    cor_mat = cor_mat, sigmaB = sigmaB, rate_list = rate_list,
    exposure = exposure, dispersion = dispersion, Eper = Eper, Eco = Eco,
    rho = rho, TAR = if (is.null(TAR)) rep(1, length(mu_list)) else TAR,
    arm_names = arm_names, ynames_list = ynames_list, type_y = type_y,
    list_comparator = list_comparator, list_y_comparator = list_y_comparator,
    alpha = alpha, lequi.tol = lequi.tol, uequi.tol = uequi.tol,
    list_lequi.tol = list_lequi.tol, list_uequi.tol = list_uequi.tol,
    dtype = dtype, ctype = ctype, vareq = vareq, k = k, adjust = adjust,
    dropout = dropout, nsim = nsim, seed = seed, ncores = ncores,
    keep_sim_data = keep_sim_data
  )
  if (is_count) {
    count_model <- if (distribution == "pois") "poisson" else
      "negative-binomial"
    count_design <- match.arg(dtype, c("parallel", "2x2"))
    count_endpoint_corr <- if (length(cor_mat) == 1L && is.na(cor_mat)) {
      NULL
    } else {
      cor_mat
    }
    count_k <- if (length(k) == 1L && is.na(k)) NULL else k
    count_adjust <- .normalize_count_adjustment(adjust)
    if (is.null(rate_list) || !is.list(rate_list))
      stop("For count outcomes, supply a named 'rate_list' with one rate vector per arm.")
    if (is.null(names(rate_list)) || any(!nzchar(names(rate_list))))
      stop("'rate_list' must be named so arm-specific parameters can be matched.")
    check_arm_parameter <- function(x, name) {
      if (is.list(x) && (is.null(names(x)) ||
                         !setequal(names(x), names(rate_list))))
        stop(sprintf("Names of arm-specific '%s' must match 'rate_list'.", name))
    }
    check_arm_parameter(exposure, "exposure")
    check_arm_parameter(dispersion, "dispersion")
    if (length(list_comparator) == 1L && is.na(list_comparator))
      stop("For count outcomes, supply 'list_comparator'.")
    if (!is.list(list_lequi.tol) || !is.list(list_uequi.tol))
      stop("For count outcomes, supply comparator-specific 'list_lequi.tol' and 'list_uequi.tol'.")
    comparisons <- list_comparator
    comparison_names <- names(comparisons)
    if (is.null(comparison_names)) {
      comparison_names <- paste0("comparison_", seq_along(comparisons))
      names(comparisons) <- comparison_names
    }
    if (is.null(names(list_lequi.tol)) || is.null(names(list_uequi.tol)) ||
        !setequal(names(list_lequi.tol), comparison_names) ||
        !setequal(names(list_uequi.tol), comparison_names))
      stop("Names of 'list_lequi.tol' and 'list_uequi.tol' must match 'list_comparator'.")

    count_all_endpoints <- names(rate_list[[1L]])
    if (is.null(count_all_endpoints))
      count_all_endpoints <- paste0("endpoint_", seq_along(rate_list[[1L]]))
    count_subset <- .prepare_count_endpoint_subset(
      rates = rate_list, comparators = comparisons,
      list_y_comparator = list_y_comparator, exposure = exposure,
      dispersion = dispersion, cor_mat = count_endpoint_corr,
      list_lequi.tol = list_lequi.tol, list_uequi.tol = list_uequi.tol
    )
    rate_list <- count_subset$rates
    exposure <- count_subset$exposure
    dispersion <- count_subset$dispersion
    count_endpoint_corr <- count_subset$cor_mat
    list_lequi.tol <- count_subset$list_lequi.tol
    list_uequi.tol <- count_subset$list_uequi.tol
    count_m <- length(count_subset$endpoints)
    count_type_info <- .prepare_type_y(
      type_y = type_y, all_endpoints = count_all_endpoints,
      selected_endpoints = count_subset$endpoints,
      adjust = if (adjust == "seq") "seq" else "no"
    )
    count_type_y <- if (count_type_info$active)
      count_type_info$type_y[count_subset$endpoints] else NULL
    count_k <- .normalize_count_k(count_k, count_m, allow_vector = TRUE)
    if (!is.null(count_k) && length(count_k) > 1L) {
      if (length(count_k) != length(comparisons) ||
          length(unique(count_k)) != 1L)
        stop("Count joint outcomes require one common 'k' across comparisons.")
      count_k <- count_k[[1L]]
    }
    if (length(comparisons) == 1L && length(rate_list) == 2L) {
      arms <- comparisons[[1L]]
      if (all(is.na(dropout))) dropout <- c(0, 0)
      out <- power_count(
        n_per_arm = as.integer(n), rate_test = rate_list[[arms[1L]]],
        rate_reference = rate_list[[arms[2L]]], exposure = exposure,
        margin_lower = list_lequi.tol[[comparison_names[1L]]],
        margin_upper = list_uequi.tol[[comparison_names[1L]]],
        model = count_model, dispersion = dispersion, alpha = alpha,
        nsim = nsim, seed = seed, design = count_design,
        k = count_k,
        endpoint_corr = count_endpoint_corr,
        type_y = count_type_y, adjust = count_adjust,
        sigmaB = sigmaB, Eper = Eper, Eco = Eco,
        dropout = dropout, ncores = ncores,
        .warn_redundant_bon = .warn_redundant_bon
      )
      out$distribution <- distribution
      out$endpoint_corr <- count_endpoint_corr
      out$n <- as.integer(n)
      out$parameters <- update_parameters
      if (keep_sim_data) {
        out$sim_data <- .retained_count_data(
          rate_list = rate_list, exposure = exposure, dispersion = dispersion,
          distribution = distribution, design = count_design,
          n = out$n_per_arm, nsim = nsim, seed = seed,
          endpoint_corr = count_endpoint_corr
        )
      }
      class(out) <- c("simpower", "countpower")
      return(out)
    }
    out <- power_count_joint(
      n_per_arm = as.integer(n), rates = rate_list,
      comparisons = comparisons, exposure = exposure,
      margin_lower = NULL, margin_upper = NULL,
      list_margin_lower = if (is.list(list_lequi.tol)) list_lequi.tol else NULL,
      list_margin_upper = if (is.list(list_uequi.tol)) list_uequi.tol else NULL,
      model = count_model, dispersion = dispersion, alpha = alpha,
      endpoint_corr = count_endpoint_corr,
      k = count_k,
      type_y = count_type_y, adjust = count_adjust, nsim = nsim, seed = seed,
      design = count_design, ncores = ncores,
      .warn_redundant_bon = .warn_redundant_bon
    )
    out$distribution <- distribution
    out$endpoint_corr <- count_endpoint_corr
    out$n <- as.integer(n)
    out$parameters <- update_parameters
    if (keep_sim_data && !is.null(out$n_per_arm)) {
      out$sim_data <- .retained_count_data(
        rate_list = rate_list, exposure = exposure, dispersion = dispersion,
        distribution = distribution, design = count_design,
        n = out$n_per_arm, nsim = nsim, seed = seed,
        endpoint_corr = count_endpoint_corr
      )
    }
      class(out) <- c("simpower", "countpower")
    return(out)
  }

  continuous <- sampleSize(
    mu_list = mu_list, varcov_list = varcov_list, sigma_list = sigma_list,
    cor_mat = cor_mat, sigmaB = sigmaB, Eper = Eper, Eco = Eco, rho = rho,
    TAR = if (is.null(TAR)) rep(1, length(mu_list)) else TAR,
    arm_names = arm_names, ynames_list = ynames_list, type_y = type_y,
    list_comparator = list_comparator,
    list_y_comparator = list_y_comparator, power = 0,
    alpha = alpha, lequi.tol = lequi.tol, uequi.tol = uequi.tol,
    list_lequi.tol = list_lequi.tol, list_uequi.tol = list_uequi.tol,
    dtype = dtype, ctype = ctype, vareq = vareq, k = k, adjust = adjust,
    dropout = dropout, nsim = nsim, seed = seed, ncores = ncores,
    optimization_method = "step-by-step", lower = as.integer(n),
    upper = as.integer(n), distribution = distribution,
    .warn_redundant_bon = .warn_redundant_bon
  )
  result <- continuous$response
  result <- as.data.frame(result)
  power_value <- as.numeric(result$power[1L])
  out <- list(power = power_value,
              power_LCI = as.numeric(result$power_LCI[1L]),
              power_UCI = as.numeric(result$power_UCI[1L]),
              n = as.integer(n), n_total = as.numeric(result$n_total[1L]),
              distribution = distribution, nsim = nsim,
              result = continuous)
  out$parameters <- update_parameters
  if (keep_sim_data) {
    out$sim_data <- .retained_continuous_data(
      param = continuous$param, param.d = continuous$param.d,
      n = as.integer(n), nsim = nsim, seed = seed
    )
  }
  class(out) <- "simpower"
  out
}

#' @export
print.simpower <- function(x, ...) {
  if (inherits(x, "simpower_curve")) {
    cat("Fixed-sample-size power curve\n")
    print(x$power_curve, row.names = FALSE)
    return(invisible(x))
  }
  cat("Fixed-sample-size power\n")
  cat("Distribution:", x$distribution, "\n")
  cat("Sample size:", x$n, "\n")
  cat(sprintf("Power: %.4f [%.4f, %.4f]\n",
              x$power, x$power_LCI, x$power_UCI))
  invisible(x)
}
