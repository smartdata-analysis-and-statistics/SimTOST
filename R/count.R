#' Estimate power for count-rate equivalence
#'
#' @param n_per_arm Subjects per arm.
#' @param design Trial design: `"parallel"` or `"2x2"`. For `"2x2"`,
#' `n_per_arm` is interpreted as subjects per sequence.
#' @param rate_test Event rate in the test arm.
#' @param rate_reference Event rate in the reference arm.
#' @param exposure Exposure per subject; a scalar or one value per endpoint.
#' @param margin_lower Lower rate-ratio margin; a scalar or one value per endpoint.
#' @param margin_upper Upper rate-ratio margin; a scalar or one value per endpoint.
#' @param model Count model: `"poisson"` or `"negative-binomial"`.
#' @param dispersion Positive negative-binomial dispersion parameter. The
#' per-subject negative-binomial size is `1 / dispersion`; parallel-arm
#' totals use size `n / dispersion`.
#' @param alpha One-sided significance level.
#' @param nsim Number of simulations.
#' @param seed Optional random seed.
#' @param k Number of endpoints that must demonstrate equivalence. Defaults to
#' all supplied endpoints.
#' @param endpoint_corr Endpoint correlation matrix used by the Gaussian
#' copula for multi-endpoint count simulations. The default is independence.
#' @param type_y Numeric endpoint hierarchy used with `adjust = "seq"`: `1`
#' for primary/co-primary endpoints and `2` for secondary endpoints. Named
#' vectors are recommended when endpoint names are available.
#' @param adjust Multiplicity adjustment for endpoint-wise one-sided alpha:
#' `"none"`, `"bonferroni"`, `"sidak"`, `"t"`, or `"seq"`/`"sequential"`.
#' The `"t"` option uses Mielke's strong k-out-of-m calibration
#' `alpha / (m - k + 1)`; legacy partial-conjunction labels are accepted.
#' Sequential testing applies the same primary-gate/secondary-family rule as
#' the continuous kernels.
#' When `k` equals the number of supplied endpoints, endpoint-wise adjustment
#' is not necessary for the all-endpoints-required intersection-union decision;
#' the requested method is retained but a warning is issued.
#' @param sigmaB Between-subject standard deviation on the log-rate scale for
#' the count `2x2` design.
#' @param Eper Numeric vector of length 2 containing period effects on the
#' log-rate scale.
#' @param Eco Numeric vector of length 2 containing carry-over effects on the
#' log-rate scale, ordered as reference carry-over and treatment carry-over.
#' @param dropout Numeric vector of length 2 containing dropout proportions
#' for the two crossover sequences.
#' @param type_y_active Internal flag indicating whether `type_y` is active.
#' @details For `design = "2x2"`, complete participants contribute one count
#' under each treatment. The kernel analyzes within-participant log-rate
#' contrasts, averages the two sequence-specific estimates to remove period
#' effects, and applies the carry-over correction implied by
#' `Eco = c(reference_carryover, treatment_carryover)`. `exposure` is used as
#' the log-rate offset. `sigmaB` is the standard deviation of a subject
#' random intercept used in the count-generating model; it cancels from the
#' within-participant treatment contrast. The standard error is estimated from
#' the empirical variance of the subject-level contrasts. Participants who
#' drop out before completing both periods do not contribute to this paired
#' analysis.
#' @return An object of class `countpower` containing estimated power and its
#' binomial confidence interval.
#' @examples
#' SimTOST:::power_count(40, 0.20, 0.20, nsim = 100, seed = 1)
.power_count_serial <- function(n_per_arm, rate_test, rate_reference, exposure = 1,
                        margin_lower = 0.80, margin_upper = 1.25,
                        model = c("poisson", "negative-binomial"),
                        dispersion = 0.1, alpha = 0.05, nsim = 5000,
                        seed = NULL, design = c("parallel", "2x2"),
                        k = NULL, endpoint_corr = NULL,
                        type_y = NULL,
                        adjust = c("none", "bonferroni", "sidak", "t", "pc",
                                    "partial-conjunction", "partial_conjunction",
                                    "sequential"),
                        sigmaB = 0, Eper = c(0, 0), Eco = c(0, 0),
                        dropout = c(0, 0), type_y_active = FALSE) {
  model <- match.arg(model)
  design <- match.arg(design)
  adjust <- .normalize_count_adjustment(match.arg(adjust))
  parameter_length <- function(x) {
    if (is.list(x)) max(vapply(x, length, integer(1))) else length(x)
  }
  m <- max(length(rate_test), length(rate_reference), parameter_length(exposure),
           length(margin_lower), length(margin_upper), length(alpha),
           parameter_length(dispersion))
  k <- .normalize_count_k(k, m)
  if (is.null(k)) k <- m
  endpoint_names <- names(rate_test)
  if (is.null(endpoint_names)) endpoint_names <- names(rate_reference)
  if (is.null(endpoint_names)) endpoint_names <- paste0("endpoint_", seq_len(m))
  if (is.null(type_y) || length(type_y) != m)
    type_y <- rep(-1L, m)
  if (isTRUE(type_y_active) && any(!type_y %in% c(1L, 2L)))
    stop("Active count endpoint hierarchies must contain only 1 or 2.")
  if (m > 1L) {
    if (xor(is.null(names(rate_test)), is.null(names(rate_reference))) ||
        (!is.null(names(rate_test)) &&
         !identical(names(rate_test), names(rate_reference))))
      stop("Endpoint names must be identical across test and reference rates.")
    if (!is.null(endpoint_names)) {
      named_parameters <- list(exposure = exposure, dispersion = dispersion,
                               margin_lower = margin_lower,
                               margin_upper = margin_upper, alpha = alpha)
      for (parameter_name in names(named_parameters)) {
        value <- named_parameters[[parameter_name]]
        if (!is.list(value) && length(value) == m &&
            !is.null(names(value)) && !identical(names(value), endpoint_names))
          stop(sprintf("Names of '%s' must match the endpoint names.",
                       parameter_name))
      }
    }
  }
  recycle_endpoint <- function(x, name) {
    if (length(x) == 1L) return(rep(x, m))
    if (length(x) != m)
      stop(sprintf("'%s' must have length 1 or %d.", name, m))
    x
  }
  rate_test <- recycle_endpoint(rate_test, "rate_test")
  rate_reference <- recycle_endpoint(rate_reference, "rate_reference")
  arm_parameter <- function(x, index, name) {
    if (is.list(x)) {
      if (length(x) != 2L)
        stop(sprintf("'%s' as a list must contain test and reference values.", name))
      value <- x[[index]]
    } else {
      value <- x
    }
    recycle_endpoint(value, name)
  }
  exposure_test <- arm_parameter(exposure, 1L, "exposure")
  exposure_reference <- arm_parameter(exposure, 2L, "exposure")
  dispersion_test <- arm_parameter(dispersion, 1L, "dispersion")
  dispersion_reference <- arm_parameter(dispersion, 2L, "dispersion")
  margin_lower <- recycle_endpoint(margin_lower, "margin_lower")
  margin_upper <- recycle_endpoint(margin_upper, "margin_upper")
  alpha <- recycle_endpoint(alpha, "alpha")
  if (length(sigmaB) != 1L || !is.finite(sigmaB) || sigmaB < 0)
    stop("'sigmaB' must be a single non-negative finite value.")
  if (length(Eper) != 2L || any(!is.finite(Eper)) ||
      length(Eco) != 2L || any(!is.finite(Eco)))
    stop("'Eper' and 'Eco' must be finite numeric vectors of length 2.")
  if (length(dropout) != 2L || any(!is.finite(dropout)) ||
      any(dropout < 0) || any(dropout >= 1))
    stop("'dropout' must contain two values in [0, 1).")
  validate_count_inputs(n_per_arm, rate_test, rate_reference,
                        c(exposure_test, exposure_reference),
                        margin_lower, margin_upper, model,
                        c(dispersion_test, dispersion_reference),
                        alpha, nsim)
  if (is.null(endpoint_corr)) endpoint_corr <- diag(m)
  endpoint_corr <- as.matrix(endpoint_corr)
  if (!is.numeric(endpoint_corr) || !all(dim(endpoint_corr) == c(m, m)) ||
      any(!is.finite(endpoint_corr)) ||
      max(abs(endpoint_corr - t(endpoint_corr))) > 1e-10 ||
      max(abs(diag(endpoint_corr) - 1)) > 1e-10 ||
      min(eigen(endpoint_corr, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop(sprintf("'endpoint_corr' must be a positive-definite %d x %d correlation matrix.", m, m))
  alpha_endpoint <- .count_endpoint_alpha(
    alpha, m, k, adjust, type_y, type_y_active
  )
  if (!is.null(seed)) set.seed(seed)
  if (m == 1L) {
    cpp_result <- count_power_cpp(
      n_per_arm, rate_test[1], rate_reference[1], exposure_test[1],
      exposure_reference[1], margin_lower[1], margin_upper[1],
      if (model == "poisson") 0L else 1L, dispersion_test[1],
      dispersion_reference[1],
      alpha_endpoint[1], nsim, if (design == "parallel") 0L else 1L,
      sigmaB, Eper, Eco, dropout
    )
  } else {
    cpp_result <- count_power_multi_cpp(
      n_per_arm, rate_test, rate_reference, exposure_test, exposure_reference,
      margin_lower, margin_upper, if (model == "poisson") 0L else 1L,
      dispersion_test, dispersion_reference,
      alpha_endpoint, nsim, if (design == "parallel") 0L else 1L,
      endpoint_corr, as.integer(type_y), isTRUE(type_y_active),
      as.integer(k), sigmaB, Eper, Eco, dropout
    )
  }
  ci <- stats::prop.test(cpp_result[["successes"]], nsim, correct = TRUE)$conf.int
  endpoint_successes <- cpp_result[grep("^endpoint_success_", names(cpp_result))]
  endpoint_successes <- unname(as.numeric(endpoint_successes))
  names(endpoint_successes) <- endpoint_names
  comparison_successes <- cpp_result[grep("^comparison_success_", names(cpp_result))]
  comparison_successes <- unname(as.numeric(comparison_successes))
  names(comparison_successes) <- "comparison_1"
  out <- list(power = cpp_result[["power"]], power_LCI = ci[1], power_UCI = ci[2],
              successes = unname(cpp_result[["successes"]]),
              endpoint_successes = endpoint_successes,
              comparison_successes = comparison_successes,
              n_per_arm = n_per_arm, n_total = 2 * n_per_arm,
              model = model, design = design, nsim = nsim,
              n_endpoints = m, k = k, endpoint_corr = endpoint_corr,
              type_y = if (type_y_active) type_y else NA,
              adjust = adjust)
  class(out) <- "countpower"
  out
}

# Run independent Monte Carlo chunks on separate R workers.  The C++ kernels
# use R's RNG, so parallelism is deliberately implemented outside the kernel:
# each worker receives its own deterministic seed and performs a serial C++
# simulation.  This is portable to macOS, Linux, and Windows.
.combine_count_power_chunks <- function(serial_fun, serial_name, args, nsim,
                                        ncores, seed) {
  if (ncores <= 1L || nsim < 2L) {
    out <- do.call(serial_fun, args)
    out$ncores <- 1L
    return(out)
  }
  ncores <- min(as.integer(ncores), as.integer(nsim))
  chunks <- rep.int(nsim %/% ncores, ncores)
  chunks[seq_len(nsim %% ncores)] <- chunks[seq_len(nsim %% ncores)] + 1L
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, ncores)

  worker <- function(i, args, chunks, seeds) {
    args$nsim <- chunks[[i]]
    args$seed <- seeds[[i]]
    serial_fun <- get(serial_name, envir = asNamespace("SimTOST"))
    do.call(serial_fun, args)
  }
  if (.Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(ncores), worker,
                                  args = args, chunks = chunks, seeds = seeds,
                                  mc.preschedule = TRUE)
  } else {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      if (!"SimTOST" %in% loadedNamespaces()) loadNamespace("SimTOST")
      NULL
    })
    results <- parallel::parLapply(cl, seq_len(ncores), worker,
                                   args = args, chunks = chunks, seeds = seeds)
  }
  successes <- sum(vapply(results, function(x) x$successes, numeric(1)))
  out <- results[[1L]]
  out$successes <- successes
  if (!is.null(out$endpoint_successes)) {
    out$endpoint_successes <- Reduce(
      `+`, lapply(results, function(x) x$endpoint_successes)
    )
  }
  if (!is.null(out$comparison_successes)) {
    out$comparison_successes <- Reduce(
      `+`, lapply(results, function(x) x$comparison_successes)
    )
  }
  out$nsim <- as.integer(nsim)
  out$power <- successes / nsim
  out$power_LCI <- stats::prop.test(successes, nsim, correct = TRUE)$conf.int[1L]
  out$power_UCI <- stats::prop.test(successes, nsim, correct = TRUE)$conf.int[2L]
  out$ncores <- ncores
  out
}

#' @rdname power_count
#' @param ncores Number of worker processes used for count simulations.
#' Set to 1 for serial execution. Parallel execution splits `nsim` into
#' reproducible independent chunks and combines the resulting successes.
power_count <- function(n_per_arm, rate_test, rate_reference, exposure = 1,
                        margin_lower = 0.80, margin_upper = 1.25,
                        model = c("poisson", "negative-binomial"),
                        dispersion = 0.1, alpha = 0.05, nsim = 5000,
                        seed = NULL, design = c("parallel", "2x2"),
                        k = NULL, endpoint_corr = NULL,
                        type_y = NULL,
                        adjust = c("none", "bonferroni", "sidak", "t", "pc",
                                    "partial-conjunction", "partial_conjunction",
                                    "sequential"),
                        sigmaB = 0, Eper = c(0, 0), Eco = c(0, 0),
                        dropout = c(0, 0), ncores = 1,
                        .warn_redundant_bon = TRUE, type_y_active = FALSE) {
  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores != as.integer(ncores) || ncores < 1L)
    stop("'ncores' must be a positive integer.")
  adjust <- .normalize_count_adjustment(match.arg(adjust))
  parameter_length <- function(x) {
    if (is.list(x)) max(vapply(x, length, integer(1))) else length(x)
  }
  m <- max(length(rate_test), length(rate_reference), parameter_length(exposure),
           length(margin_lower), length(margin_upper), length(alpha),
           parameter_length(dispersion))
  k <- .normalize_count_k(k, m)
  if (is.null(k)) k <- m
  endpoint_names <- names(rate_test)
  if (is.null(endpoint_names)) endpoint_names <- names(rate_reference)
  if (is.null(endpoint_names)) endpoint_names <- paste0("endpoint_", seq_len(m))
  type_info <- .prepare_type_y(
    type_y = type_y, all_endpoints = endpoint_names,
    selected_endpoints = endpoint_names,
    adjust = if (identical(adjust, "sequential")) "seq" else "no"
  )
  if (isTRUE(.warn_redundant_bon))
    .warn_adjustment_configuration(
      k = k, m = m,
      adjust = if (identical(adjust, "sequential")) "seq" else adjust,
      type_y = if (type_info$active) type_info$type_y else NULL,
      type_y_supplied = type_info$active,
      context = "count endpoints"
    )
  args <- list(n_per_arm = n_per_arm, rate_test = rate_test,
               rate_reference = rate_reference, exposure = exposure,
               margin_lower = margin_lower, margin_upper = margin_upper,
               model = model, dispersion = dispersion, alpha = alpha,
               nsim = nsim, seed = seed, design = design, k = k,
               endpoint_corr = endpoint_corr, type_y = type_info$type_y,
               type_y_active = type_info$active, adjust = adjust,
               sigmaB = sigmaB, Eper = Eper, Eco = Eco, dropout = dropout)
  out <- .combine_count_power_chunks(.power_count_serial, ".power_count_serial", args, nsim,
                                     as.integer(ncores), seed)
  out$parameters <- c(args, list(ncores = ncores, .function = "power_count"))
  out
}

#' Estimate sample size for count-rate equivalence
#'
#' @description Searches for the smallest number of subjects per arm whose
#' simulated power reaches the target for a rate-ratio equivalence test.
#' @inheritParams power_count
#' @param power Target power.
#' @param rate_test Event rate in the test arm.
#' @param rate_reference Event rate in the reference arm.
#' @param exposure Exposure per subject; a scalar or one value per endpoint.
#' @param margin_lower Lower rate-ratio equivalence margin.
#' @param margin_upper Upper rate-ratio equivalence margin.
#' @param model Count model: `"poisson"` or `"negative-binomial"`.
#' @param dispersion Positive negative-binomial dispersion parameter. The
#' per-subject negative-binomial size is `1 / dispersion`; parallel-arm
#' totals use size `n / dispersion`.
#' @param alpha One-sided significance level.
#' @param nsim Number of simulated trials.
#' @param seed Optional random seed.
#' @param design Trial design: `"parallel"` or `"2x2"`. In a crossover
#' design, the returned sample size is per sequence.
#' @param lower Minimum subjects per arm.
#' @param upper Maximum subjects per arm.
#' @param k Number of endpoints that must demonstrate equivalence. Defaults to
#' all supplied endpoints.
#' @param endpoint_corr Endpoint correlation matrix used by the Gaussian
#' copula for multi-endpoint count simulations. The default is independence.
#' @param type_y Numeric endpoint hierarchy used with `adjust = "seq"`: `1`
#' for primary/co-primary endpoints and `2` for secondary endpoints.
#' @param adjust Multiplicity adjustment for endpoint-wise one-sided alpha:
#' `"none"`, `"bonferroni"`, `"sidak"`, `"t"`, or `"seq"`/`"sequential"`.
#' The `"t"` option uses Mielke's strong k-out-of-m calibration
#' `alpha / (m - k + 1)`; legacy partial-conjunction labels are accepted.
#' @param optimization_method Search method. `"fast"` brackets the power
#' crossing and uses integer bisection; `"step-by-step"` evaluates every
#' candidate sample size.
#' @param step.power Initial power-of-two jump used by the fast search.
#' @param step.up Direction of the initial bracketing search.
#' @param pos.side Retained for compatibility with `sampleSize()`; count
#' searches return the smallest candidate reaching the target.
#' @param maxiter Maximum number of power evaluations.
#' @param sigmaB Between-subject standard deviation for the count 2x2 design.
#' @param Eper Period effects for the count 2x2 design.
#' @param Eco Carry-over effects for the count 2x2 design.
#' @param dropout Dropout proportions for the count 2x2 design.
#' @param .warn_redundant_bon Logical. If `TRUE`, warn about redundant or
#'   uncalibrated multiplicity configurations.
#' @return An object of class `countss` containing the selected sample size,
#' achieved power, confidence interval, input parameters, and the search
#' history in `table.iter` and `table.test`. For count outcomes, `table.iter`
#' has one row per evaluated candidate sample size. `table.test` contains
#' complete-trial, comparator, and endpoint decision indicators for each
#' simulated trial and candidate. The count kernel returns aggregate decision
#' counts rather than raw endpoint-level test statistics, so component columns
#' preserve the simulated marginal success counts.
#' @examples
#' SimTOST:::sampleSize_count(0.80, 0.20, 0.20, lower = 100, upper = 2000,
#'                  nsim = 100, seed = 1)
sampleSize_count <- function(power = 0.80, rate_test, rate_reference,
                              exposure = 1, margin_lower = 0.80,
                              margin_upper = 1.25,
                              model = c("poisson", "negative-binomial"),
                              dispersion = 0.1, alpha = 0.05, nsim = 5000,
                              seed = NULL, lower = 2, upper = 500,
                              design = c("parallel", "2x2"), k = NULL,
                              endpoint_corr = NULL,
                              type_y = NULL,
                              adjust = c("none", "bonferroni", "sidak", "t", "pc",
                                          "partial-conjunction", "partial_conjunction",
                                          "sequential"),
                              sigmaB = 0, Eper = c(0, 0), Eco = c(0, 0),
                              dropout = c(0, 0),
                              optimization_method = c("fast", "step-by-step"),
                              step.power = 6, step.up = TRUE,
                              pos.side = FALSE, maxiter = 1000, ncores = 1,
                              .warn_redundant_bon = TRUE) {
  model <- match.arg(model)
  design <- match.arg(design)
  adjust <- .normalize_count_adjustment(match.arg(adjust))
  optimization_method <- match.arg(optimization_method)
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("'power' must be a single number between 0 and 1.")
  if (upper < lower) stop("'upper' must be greater than or equal to 'lower'.")
  if (lower != as.integer(lower) || upper != as.integer(upper) || lower < 2)
    stop("'lower' and 'upper' must be integers with lower >= 2.")
  parameter_length <- function(x) {
    if (is.list(x)) max(vapply(x, length, integer(1))) else length(x)
  }
  m <- max(length(rate_test), length(rate_reference), parameter_length(exposure),
           length(margin_lower), length(margin_upper), length(alpha),
           parameter_length(dispersion))
  k <- .normalize_count_k(k, m)
  if (is.null(k)) k <- m
  endpoint_names <- names(rate_test)
  if (is.null(endpoint_names)) endpoint_names <- names(rate_reference)
  if (is.null(endpoint_names)) endpoint_names <- paste0("endpoint_", seq_len(m))
  type_info <- .prepare_type_y(
    type_y = type_y, all_endpoints = endpoint_names,
    selected_endpoints = endpoint_names,
    adjust = if (identical(adjust, "sequential")) "seq" else "no"
  )
  if (isTRUE(.warn_redundant_bon))
    .warn_adjustment_configuration(
      k = k, m = m,
      adjust = if (identical(adjust, "sequential")) "seq" else adjust,
      type_y = if (type_info$active) type_info$type_y else NULL,
      type_y_supplied = type_info$active,
      context = "count endpoints"
    )
  power_fun <- function(n) {
    power_count(
      n_per_arm = n, rate_test = rate_test,
      rate_reference = rate_reference, exposure = exposure,
      margin_lower = margin_lower, margin_upper = margin_upper,
      model = model, dispersion = dispersion, alpha = alpha, nsim = nsim,
      seed = seed, design = design, k = k, endpoint_corr = endpoint_corr,
      type_y = type_y, adjust = adjust,
      .warn_redundant_bon = FALSE,
      sigmaB = sigmaB, Eper = Eper, Eco = Eco, dropout = dropout,
      ncores = ncores
    )
  }
  search <- find_count_sample_size(
    power_fun = power_fun, target_power = power, lower = lower, upper = upper,
    optimization_method = optimization_method, step.power = step.power,
    step.up = step.up, pos.side = pos.side, maxiter = maxiter,
    return_history = TRUE
  )
  selected <- search$selected
  achieved <- search$results[[match(selected, search$n_per_arm)]]
  history <- .count_search_tables(search)
  out <- c(achieved, list(target_power = power, parameters = list(
    rate_test = rate_test, rate_reference = rate_reference,
    exposure = exposure, margin_lower = margin_lower,
    margin_upper = margin_upper, model = model, dispersion = dispersion,
    alpha = alpha, lower = lower, upper = upper, design = design,
    type_y = type_y,
    k = achieved$k, endpoint_corr = endpoint_corr, adjust = adjust,
    sigmaB = sigmaB, Eper = Eper,
    Eco = Eco, dropout = dropout, optimization_method = optimization_method,
    step.power = step.power, step.up = step.up, pos.side = pos.side,
    maxiter = maxiter)))
  class(out) <- "countss"
  out$parameters <- list(
    rate_test = rate_test, rate_reference = rate_reference,
    exposure = exposure, margin_lower = margin_lower,
    margin_upper = margin_upper, model = model, dispersion = dispersion,
    alpha = alpha, nsim = nsim, seed = seed, design = design, k = out$k,
    endpoint_corr = endpoint_corr, type_y = type_y, adjust = adjust, sigmaB = sigmaB,
    Eper = Eper, Eco = Eco, dropout = dropout, ncores = ncores,
    power = power, lower = lower, upper = upper,
    optimization_method = optimization_method, step.power = step.power,
    step.up = step.up, pos.side = pos.side, maxiter = maxiter,
    .function = "sampleSize_count"
  )
  out$table.iter <- history$table.iter
  out$table.test <- history$table.test
  out
}

find_count_sample_size <- function(power_fun, target_power, lower, upper,
                                   optimization_method = c("fast", "step-by-step"),
                                   step.power = 6, step.up = TRUE,
                                   pos.side = FALSE, maxiter = 1000,
                                   return_history = FALSE) {
  optimization_method <- match.arg(optimization_method)
  if (!is.numeric(step.power) || length(step.power) != 1L ||
      !is.finite(step.power) || step.power < 0)
    stop("'step.power' must be a single non-negative finite number.")
  if (!is.logical(step.up) || length(step.up) != 1L || is.na(step.up) ||
      !is.logical(pos.side) || length(pos.side) != 1L || is.na(pos.side))
    stop("'step.up' and 'pos.side' must be single logical values.")
  if (!is.numeric(maxiter) || length(maxiter) != 1L ||
      !is.finite(maxiter) || maxiter < 1 || maxiter != as.integer(maxiter))
    stop("'maxiter' must be a positive integer.")
  if (!is.logical(return_history) || length(return_history) != 1L ||
      is.na(return_history))
    stop("'return_history' must be a single logical value.")

  cache <- new.env(parent = emptyenv())
  result_cache <- new.env(parent = emptyenv())
  evaluations <- 0L
  evaluate <- function(n) {
    key <- as.character(as.integer(n))
    if (exists(key, envir = cache, inherits = FALSE))
      return(get(key, envir = cache, inherits = FALSE))
    if (evaluations >= maxiter)
      stop("reached 'maxiter' without finding a sample size.")
    result <- power_fun(as.integer(n))
    value <- result$power
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
      stop("The count power calculation returned an invalid power value.")
    evaluations <<- evaluations + 1L
    assign(key, value, envir = cache)
    assign(key, result, envir = result_cache)
    value
  }

  finish <- function(selected) {
    selected <- as.integer(selected)
    if (!isTRUE(return_history)) return(selected)
    n_evaluated <- sort(as.integer(ls(envir = cache, all.names = TRUE)))
    list(
      selected = selected,
      evaluations = evaluations,
      n_per_arm = n_evaluated,
      results = lapply(as.character(n_evaluated), function(key) {
        get(key, envir = result_cache, inherits = FALSE)
      })
    )
  }

  p_lower <- evaluate(lower)
  if (p_lower >= target_power) return(finish(lower))
  p_upper <- evaluate(upper)
  if (p_upper < target_power)
    stop(sprintf("Target power %.3f was not reached by upper = %d; increase 'upper'.",
                 target_power, upper))
  if (optimization_method == "step-by-step") {
    for (n in seq.int(lower + 1L, upper)) {
      if (evaluate(n) >= target_power) return(finish(n))
    }
  }

  # The fast search first brackets the crossing using a power-of-two jump,
  # then uses integer bisection. Common random numbers (the same seed at each
  # candidate n) make this substantially more stable for simulated power.
  if (step.up) {
    left <- as.integer(lower)
    right <- min(as.integer(upper),
                 left + max(1L, as.integer(2^step.power)))
    while (evaluate(right) < target_power && right < upper) {
      jump <- max(1L, right - left)
      left <- right
      right <- min(as.integer(upper), right + max(1L, 2L * jump))
    }
  } else {
    right <- as.integer(upper)
    jump <- max(1L, as.integer(2^step.power))
    left <- max(as.integer(lower), right - jump)
    while (evaluate(left) >= target_power && left > lower) {
      right <- left
      jump <- max(1L, 2L * jump)
      left <- max(as.integer(lower), right - jump)
    }
  }

  while (right - left > 1L) {
    middle <- left + (right - left) %/% 2L
    if (evaluate(middle) >= target_power) right <- middle else left <- middle
  }
  # pos.side is retained for compatibility with sampleSize(); count planning
  # always returns the smallest integer whose estimated power reaches target.
  finish(right)
}

# The count C++ kernels return aggregate component and complete-trial successes
# rather than the raw endpoint-level statistics retained by the continuous
# simulation. Reconstruct decision indicators with the same marginal success
# counts so count sample-size objects expose compatible search-history fields.
# A common deterministic permutation is applied to each candidate's indicators:
# the kernel does not return trial-level decisions, so placing all successes
# before all failures would create an artificial block-shaped stability plot.
.count_search_tables <- function(search) {
  if (!is.list(search) || !length(search$results) ||
      length(search$n_per_arm) != length(search$results))
    stop("The count sample-size search did not retain a valid history.")

  rows_iter <- lapply(seq_along(search$results), function(i) {
    result <- search$results[[i]]
    data.frame(
      n_iter = as.integer(search$n_per_arm[[i]]),
      n_total = as.numeric(result$n_total),
      power = as.numeric(result$power),
      power_LCI = as.numeric(result$power_LCI),
      power_UCI = as.numeric(result$power_UCI),
      stringsAsFactors = FALSE
    )
  })
  table.iter <- do.call(rbind, rows_iter)
  rownames(table.iter) <- NULL

  indicator_permutation <- function(nsim, candidate) {
    if (nsim < 2L) return(seq_len(nsim))
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(104729L + as.integer(candidate))
    sample.int(nsim)
  }

  rows_test <- lapply(seq_along(search$results), function(i) {
    result <- search$results[[i]]
    nsim <- as.integer(result$nsim)
    successes <- as.integer(result$successes)
    endpoint_successes <- result$endpoint_successes
    endpoint_successes_numeric <- as.numeric(endpoint_successes)
    comparison_names <- names(result$comparison_successes)
    comparison_successes <- as.numeric(result$comparison_successes)
    names(comparison_successes) <- comparison_names
    if (length(nsim) != 1L || length(successes) != 1L ||
        successes < 0L || successes > nsim ||
        length(endpoint_successes_numeric) == 0L ||
        any(!is.finite(endpoint_successes_numeric)) ||
        any(endpoint_successes_numeric < 0 |
            endpoint_successes_numeric > nsim) ||
        length(comparison_successes) == 0L ||
        any(!is.finite(comparison_successes)) ||
        any(comparison_successes < successes | comparison_successes > nsim))
      stop("The count power result contains invalid simulation counts.")
    data <- data.frame(
      n_iter = rep.int(as.integer(search$n_per_arm[[i]]), nsim),
      n_total = rep.int(as.numeric(result$n_total), nsim),
      stringsAsFactors = FALSE
    )
    permutation <- indicator_permutation(nsim, i)
    indicator_values <- function(success_count) {
      c(rep.int(1L, success_count), rep.int(0L, nsim - success_count))[permutation]
    }
    endpoint_is_matrix <- is.matrix(endpoint_successes)
    if (endpoint_is_matrix &&
        !identical(dim(endpoint_successes),
                   c(length(comparison_successes),
                     ncol(endpoint_successes))))
      stop("The count power result has incompatible component dimensions.")
    comparison_names <- names(comparison_successes)
    if (is.null(comparison_names))
      comparison_names <- paste0("comparison_", seq_along(comparison_successes))
    endpoint_names <- if (endpoint_is_matrix) colnames(endpoint_successes) else
      names(endpoint_successes)
    if (is.null(endpoint_names))
      endpoint_names <- paste0("endpoint_", seq_len(if (endpoint_is_matrix)
        ncol(endpoint_successes) else length(endpoint_successes)))
    for (comparison in seq_along(comparison_successes)) {
      comparison_label <- comparison_names[[comparison]]
      data[[paste0("totalyComp:", comparison_label)]] <-
        indicator_values(comparison_successes[[comparison]])
      for (endpoint in seq_along(endpoint_names)) {
        component_count <- if (endpoint_is_matrix)
          endpoint_successes[comparison, endpoint] else
          endpoint_successes[[endpoint]]
        data[[paste0(endpoint_names[[endpoint]], "Comp:", comparison_label)]] <-
          indicator_values(component_count)
      }
    }
    data$totaly <- indicator_values(successes)
    data
  })
  table.test <- do.call(rbind, rows_test)
  rownames(table.test) <- NULL

  list(table.iter = table.iter, table.test = table.test)
}

validate_count_inputs <- function(n_per_arm, rate_test, rate_reference,
                                  exposure, margin_lower, margin_upper,
                                  model, dispersion, alpha, nsim) {
  if (!is.numeric(n_per_arm) || length(n_per_arm) != 1 || n_per_arm < 2 ||
      n_per_arm != as.integer(n_per_arm)) stop("'n_per_arm' must be an integer >= 2.")
  if (any(!is.finite(c(rate_test, rate_reference, exposure))) ||
      any(c(rate_test, rate_reference, exposure) <= 0))
    stop("Rates and exposure must be positive finite numbers.")
  if (length(margin_lower) < 1 || length(margin_upper) < 1 ||
      any(!is.finite(c(margin_lower, margin_upper))) ||
      any(margin_lower <= 0) || any(margin_lower >= margin_upper))
    stop("Margins must satisfy 0 < margin_lower < margin_upper.")
  if (model == "negative-binomial" &&
      (any(!is.finite(dispersion)) || any(dispersion <= 0)))
    stop("'dispersion' must be positive for the negative-binomial model.")
  if (any(!is.finite(alpha)) || any(alpha <= 0) || any(alpha >= 1) ||
      length(nsim) != 1 || nsim < 1 || nsim != as.integer(nsim))
    stop("'alpha' must be in (0, 1) and 'nsim' must be a positive integer.")
  invisible(TRUE)
}
