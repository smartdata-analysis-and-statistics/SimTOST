#' Print count-outcome power results
#' @param x A `countpower` object.
#' @param ... Unused additional arguments.
#' @export
#' @method print countpower
print.countpower <- function(x, ...) {
  cat("Count-rate equivalence power\n")
  cat("Model:", x$model, "\n")
  cat("Design:", x$design, "\n")
  cat("Subjects per arm:", x$n_per_arm, "\n")
  if (!is.null(x$n_endpoints)) {
    cat("Endpoints:", x$n_endpoints, "(required:", x$k, ")\n")
    cat("Alpha adjustment:", x$adjust, "\n")
  }
  cat(sprintf("Power: %.4f [%.4f, %.4f]\n", x$power,
              x$power_LCI, x$power_UCI))
  invisible(x)
}

#' Summarize count-outcome power results
#' @param object An object of class `countpower`.
#' @param ... Unused additional arguments.
#' @return A data frame containing the design, model, sample size, estimated
#' power, and Monte Carlo confidence interval.
#' @export
#' @method summary countpower
summary.countpower <- function(object, ...) {
  data.frame(
    model = object$model,
    ctype = "RR",
    estimand = "event-rate ratio",
    hypothesis = "H0: RR <= L or RR >= U; H1: L < RR < U",
    design = object$design,
    n_per_arm = object$n_per_arm,
    n_total = object$n_total,
    power = object$power,
    power_LCI = object$power_LCI,
    power_UCI = object$power_UCI,
    nsim = object$nsim,
    row.names = NULL
  )
}

#' Print count-outcome sample-size results
#' @param x A `countss` object.
#' @param ... Unused additional arguments.
#' @export
#' @method print countss
print.countss <- function(x, ...) {
  cat("Count-rate equivalence sample size\n")
  cat("Subjects per arm:", x$n_per_arm, "\n")
  cat("Total subjects:", x$n_total, "\n")
  if (!is.null(x$n_endpoints)) {
    cat("Endpoints:", x$n_endpoints, "(required:", x$k, ")\n")
    cat("Alpha adjustment:", x$adjust, "\n")
  }
  cat(sprintf("Achieved power: %.4f [%.4f, %.4f]\n", x$power,
              x$power_LCI, x$power_UCI))
  invisible(x)
}

# Add the same structural result contract used by the continuous sample-size
# engine to a result produced by the unified count entry point. The legacy
# count fields are deliberately retained so existing count-specific code keeps
# working.
.decorate_count_sample_size <- function(result, distribution, rate_list,
                                        exposure, dispersion, comparisons,
                                        selected_endpoints,
                                        requested_endpoints, lower_list,
                                        upper_list, type_y, endpoint_corr,
                                        design, power, alpha, k, adjust,
                                        dropout, sigmaB, Eper, Eco, nsim,
                                        seed, ncores, optimization_method,
                                        lower, upper, step.power, step.up,
                                        pos.side, maxiter, verbose,
                                        keep_sim_data) {
  if (!is.list(result) || !inherits(result, "countss"))
    stop("'result' must be a countss object.")
  if (!is.list(rate_list) || !length(rate_list) ||
      is.null(names(rate_list)))
    stop("'rate_list' must be a named non-empty list.")
  if (!is.list(comparisons) || !length(comparisons))
    stop("'comparisons' must be a non-empty list.")

  arms <- names(rate_list)
  n <- as.numeric(result$n_per_arm)
  if (length(n) != 1L || !is.finite(n))
    stop("The count sample-size result has an invalid selected sample size.")

  dropout_for_arms <- rep(0, length(arms))
  if (is.numeric(dropout) && length(dropout) == length(arms) &&
      all(is.finite(dropout))) {
    dropout_for_arms <- as.numeric(dropout)
    if (!is.null(names(dropout)) && all(arms %in% names(dropout)))
      dropout_for_arms <- as.numeric(dropout[arms])
  }
  observed_per_arm <- ceiling(n * (1 - dropout_for_arms))
  n_drop <- sum(n - observed_per_arm)

  if (identical(design, "2x2")) {
    response <- data.table::as.data.table(list(
      n_iter = as.integer(n),
      n_drop = as.numeric(n_drop),
      n_seq0 = observed_per_arm[[1L]],
      n_seq1 = observed_per_arm[[min(2L, length(observed_per_arm))]],
      n_total = as.numeric(result$n_total),
      power = as.numeric(result$power),
      power_LCI = as.numeric(result$power_LCI),
      power_UCI = as.numeric(result$power_UCI)
    ))
  } else {
    arm_sizes <- stats::setNames(
      as.numeric(observed_per_arm), paste0("n_", arms)
    )
    response <- data.table::as.data.table(c(
      list(n_iter = as.integer(n), n_drop = as.numeric(n_drop)),
      as.list(arm_sizes),
      list(n_total = as.numeric(result$n_total),
           power = as.numeric(result$power),
           power_LCI = as.numeric(result$power_LCI),
           power_UCI = as.numeric(result$power_UCI))
    ))
  }

  param.u <- list(
    distribution = distribution, rate_list = rate_list,
    exposure = exposure, dispersion = dispersion,
    sigmaB = sigmaB, Eper = Eper, Eco = Eco,
    list_comparator = comparisons,
    list_y_comparator = requested_endpoints,
    list_lequi.tol = lower_list, list_uequi.tol = upper_list,
    type_y = type_y, endpoint_corr = endpoint_corr
  )
  param <- list(
    distribution = distribution, rate_list = rate_list,
    exposure = exposure, dispersion = dispersion,
    sigmaB = sigmaB, Eper = Eper, Eco = Eco,
    list_comparator = comparisons,
    list_y_comparator = selected_endpoints,
    list_lequi.tol = lower_list, list_uequi.tol = upper_list,
    type_y = type_y, endpoint_corr = endpoint_corr
  )
  param.d <- list(
    nsim = nsim, power = power, alpha = alpha, dtype = design,
    ctype = "RR", distribution = distribution, vareq = NA,
    k = k, adjust = adjust, dropout = dropout,
    list_lequi.tol = lower_list, list_uequi.tol = upper_list,
    seed = seed, ncores = ncores,
    optimization_method = optimization_method, lower = lower, upper = upper,
    step.power = step.power, step.up = step.up, pos.side = pos.side,
    maxiter = maxiter, verbose = verbose
  )

  # Keep one unambiguous parameter record. In particular, update() must call
  # the unified sampleSize() entry point so comparator-specific endpoint
  # selections are not lost.
  result$response <- response
  result$param.u <- param.u
  result$param <- param
  result$param.d <- param.d
  result$target_power <- power
  result$parameters <- list(
    distribution = distribution, rate_list = rate_list,
    exposure = exposure, dispersion = dispersion,
    sigmaB = sigmaB, Eper = Eper, Eco = Eco,
    cor_mat = endpoint_corr, list_comparator = comparisons,
    list_y_comparator = selected_endpoints, power = power, alpha = alpha,
    list_lequi.tol = lower_list, list_uequi.tol = upper_list,
    dtype = design, ctype = "RR", k = k, adjust = adjust,
    dropout = dropout, nsim = nsim, seed = seed, ncores = ncores,
    optimization_method = optimization_method, lower = lower, upper = upper,
    step.power = step.power, step.up = step.up, pos.side = pos.side,
    maxiter = maxiter, verbose = verbose,
    keep_sim_data = keep_sim_data, .function = "sampleSize"
  )
  result
}

#' @name summary.countss
#' @title Summary for Count Sample-Size Results
#' @param object A `countss` object returned by a count sample-size
#' calculation or by the unified `sampleSize()` function.
#' @param ... Unused additional arguments.
#' @description Prints the same design-oriented sample-size report used for
#' continuous outcomes. Count-specific result fields and search histories are
#' retained in the object, while the invisible return value is a data frame of
#' the selected per-arm (or per-sequence) and total sample sizes.
#' @return Invisibly, a data frame containing the selected sample size for each
#' arm (or sequence) and the total sample size.
#' @export
#' @method summary countss
summary.countss <- function(object, ...) {
  if (!is.list(object) || !inherits(object, "countss"))
    stop("Input must be of class 'countss'.")

  parameters <- object$parameters
  param <- object$param
  param.d <- object$param.d
  unified <- is.list(param) && is.list(param.d)
  first_non_null <- function(...) {
    values <- list(...)
    for (value in values) if (!is.null(value)) return(value)
    NULL
  }

  if (unified) {
    distribution <- param.d$distribution
    design <- param.d$dtype
    alpha <- param.d$alpha
    target_power <- param.d$power
    adjustment <- param.d$adjust
    k <- param.d$k
    comparator_list <- param$list_comparator
    endpoint_list <- param$list_y_comparator
    lower_list <- param.d$list_lequi.tol
    upper_list <- param.d$list_uequi.tol
    response <- object$response
  } else {
    if (!is.list(parameters))
      stop("The countss object does not retain its calculation parameters.")
    distribution <- if (identical(object$model, "poisson")) "pois" else
      "nbinom"
    design <- object$design
    alpha <- first_non_null(parameters$alpha, NA_real_)
    target_power <- first_non_null(object$target_power, parameters$power,
                                   NA_real_)
    adjustment <- first_non_null(object$adjust, parameters$adjust, "none")
    k <- first_non_null(object$k, parameters$k)
    comparator_list <- parameters$comparisons
    if (is.null(comparator_list))
      comparator_list <- list(comparison_1 = c("test", "reference"))
    endpoint_list <- lapply(seq_along(comparator_list), function(i) {
      value <- parameters$rate_test
      if (!is.null(names(value))) names(value) else
        paste0("endpoint_", seq_along(value))
    })
    names(endpoint_list) <- names(comparator_list)
    lower_list <- parameters$list_margin_lower
    upper_list <- parameters$list_margin_upper
    if (is.null(lower_list)) lower_list <- list(parameters$margin_lower)
    if (is.null(upper_list)) upper_list <- list(parameters$margin_upper)
    response <- NULL
  }

  if (is.null(distribution)) distribution <- "pois"
  if (is.null(design)) design <- first_non_null(object$design, "parallel")
  if (is.null(comparator_list) || !length(comparator_list))
    comparator_list <- list(comparison_1 = c("test", "reference"))
  if (is.null(names(comparator_list)))
    names(comparator_list) <- paste0("comparison_", seq_along(comparator_list))
  if (is.null(endpoint_list) || !length(endpoint_list))
    endpoint_list <- lapply(comparator_list, function(x) "endpoint_1")
  if (is.null(names(endpoint_list))) names(endpoint_list) <- names(comparator_list)

  margin_table <- do.call(rbind, lapply(seq_along(comparator_list), function(i) {
    lower <- first_non_null(lower_list[[names(comparator_list)[[i]]]],
                            lower_list[[i]])
    upper <- first_non_null(upper_list[[names(comparator_list)[[i]]]],
                            upper_list[[i]])
    lower <- as.numeric(lower)
    upper <- as.numeric(upper)
    endpoints <- names(lower)
    if (is.null(endpoints)) endpoints <- names(upper)
    if (is.null(endpoints)) endpoints <- endpoint_list[[i]]
    if (is.null(endpoints) || length(endpoints) != length(lower))
      endpoints <- paste0("endpoint_", seq_along(lower))
    data.frame(
      Comparison = paste(comparator_list[[i]], collapse = "_vs_"),
      Endpoint = endpoints, Lower = lower, Upper = upper,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(response)) {
    response <- data.frame(
      n_per_arm = object$n_per_arm, n_total = object$n_total,
      power = object$power, power_LCI = object$power_LCI,
      power_UCI = object$power_UCI, stringsAsFactors = FALSE
    )
  }
  response <- as.data.frame(response)
  sample_columns <- intersect(
    c(setdiff(grep("^n_", names(response), value = TRUE),
              c("n_iter", "n_drop")), "n_per_arm", "n_total"),
    names(response)
  )
  sample_size <- response[, unique(sample_columns), drop = FALSE]
  names(sample_size) <- sub("^n_", "", names(sample_size))
  names(sample_size)[names(sample_size) == "total"] <- "Total"

  cat("Sample Size Summary\n")
  cat(strrep("-", 22), "\n")
  cat("Design type        :", design, "\n")
  cat("Distribution       :", distribution, "\n")
  cat("Comparison type    : RR\n")
  cat("Estimand           : event-rate ratio (lambda_T / lambda_R)\n")
  cat("Hypotheses         : H0: RR <= L or >= U; H1: L < RR < U\n")
  cat("Alpha              :", alpha, "\n")
  cat("Target power       :", sprintf("%.4f", as.numeric(target_power)), "\n")
  cat("Achieved power     :", sprintf("%.4f", as.numeric(object$power)), "\n")
  cat("Power interval     :", sprintf("[%.4f, %.4f]", as.numeric(object$power_LCI),
                                       as.numeric(object$power_UCI)), "\n")
  cat("Required endpoints :", paste(k, collapse = ", "), "\n")
  cat("Alpha adjustment   :", adjustment, "\n")

  cat("\nEquivalence Margins:\n")
  print(margin_table, row.names = FALSE)
  cat("\nEstimated Sample Size:\n")
  print(sample_size, row.names = FALSE)

  invisible(sample_size)
}

#' Extract the Monte Carlo confidence interval from count power results
#'
#' @param object An object returned by a count power calculation.
#' @param parm Unused; included for compatibility with [stats::confint()].
#' @param level Confidence level for the Monte Carlo interval.
#' @param ... Unused additional arguments.
#' @return A named vector containing the lower and upper interval limits.
#' @export
#' @method confint countpower
confint.countpower <- function(object, parm = NULL, level = 0.95, ...) {
  if (!is.null(parm) && missing(level) && is.numeric(parm) &&
      length(parm) == 1L && is.finite(parm) && parm > 0 && parm < 1) {
    level <- parm
    parm <- NULL
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number between 0 and 1.")
  interval <- if (!is.null(object$successes) && !is.null(object$nsim)) {
    stats::binom.test(object$successes, object$nsim,
                      conf.level = level)$conf.int
  } else {
    c(object$power_LCI, object$power_UCI)
  }
  stats::setNames(interval, c("Lower", "Upper"))
}

#' @rdname confint.countpower
#' @export
#' @method confint countss
confint.countss <- function(object, parm = NULL, level = 0.95, ...) {
  confint.countpower(object, parm = parm, level = level, ...)
}
