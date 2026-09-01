#' Prepare inputs for a joint multi-arm count simulation
#'
#' @keywords internal
prepare_joint_count_inputs <- function(rates, comparisons, exposure,
                                        margin_lower, margin_upper, alpha,
                                        endpoint_corr,
                                        list_margin_lower = NULL,
                                        list_margin_upper = NULL,
                                        dispersion = 0.1) {
  if (!is.list(rates) || length(rates) < 2L || is.null(names(rates)) ||
      any(!nzchar(names(rates))))
    stop("'rates' must be a named list with at least two treatment arms.")
  if (any(vapply(rates, function(x) !is.numeric(x) || length(x) < 1L,
                logical(1))))
    stop("Each element of 'rates' must be a non-empty numeric vector.")
  m <- length(rates[[1L]])
  if (any(vapply(rates, length, integer(1)) != m))
    stop("All arms in 'rates' must have the same number of endpoints.")
  endpoint_names <- names(rates[[1L]])
  rate_names <- lapply(rates, names)
  if (any(vapply(rate_names, is.null, logical(1)))) {
    if (any(!vapply(rate_names, is.null, logical(1))))
      stop("Either all rate vectors must be named or none may be named.")
  } else if (any(vapply(rate_names[-1L], function(x)
                         !identical(x, endpoint_names), logical(1)))) {
    stop("Endpoint names must be identical and in the same order across arms.")
  }
  rates_mat <- do.call(rbind, lapply(rates, function(x) unname(as.numeric(x))))
  rownames(rates_mat) <- names(rates)
  if (is.null(endpoint_names)) endpoint_names <- paste0("endpoint_", seq_len(m))
  colnames(rates_mat) <- endpoint_names

  if (!is.list(comparisons) || length(comparisons) < 1L)
    stop("'comparisons' must be a non-empty list of two arm names.")
  if (is.null(names(comparisons)))
    names(comparisons) <- paste0("comparison_", seq_along(comparisons))
  if (anyDuplicated(names(comparisons)))
    stop("Comparison names must be unique.")
  if (any(vapply(comparisons, length, integer(1)) != 2L))
    stop("Each comparison must contain exactly two arm names.")
  if (any(!vapply(comparisons, function(x) all(x %in% names(rates)),
                  logical(1))))
    stop("Every comparison arm must be present in 'rates'.")
  comparison_matrix <- matrix(
    match(unlist(comparisons), names(rates)), ncol = 2L, byrow = TRUE
  ) - 1L

  recycle_endpoint <- function(x, name) {
    if (!is.numeric(x) || !(length(x) %in% c(1L, m)))
      stop(sprintf("'%s' must have length 1 or %d.", name, m))
    if (length(x) == 1L) rep(x, m) else {
      if (!is.null(names(x)) && !identical(names(x), endpoint_names))
        stop(sprintf("Names of '%s' must match the endpoint names.", name))
      as.numeric(x)
    }
  }

  arm_endpoint_matrix <- function(x, name) {
    if (is.list(x)) {
      if (length(x) != length(rates) || is.null(names(x)) ||
          !setequal(names(x), names(rates)))
        stop(sprintf("'%s' as a list must contain one named vector per arm.", name))
      out <- matrix(NA_real_, nrow = length(rates), ncol = m,
                    dimnames = list(names(rates), endpoint_names))
      for (i in seq_along(rates)) {
        value <- x[[names(rates)[i]]]
        if (!is.numeric(value) || !(length(value) %in% c(1L, m)))
          stop(sprintf("Each '%s' arm value must have length 1 or %d.", name, m))
        if (length(value) > 1L && !is.null(names(value))) {
          if (!identical(names(value), endpoint_names))
            stop(sprintf("Names of '%s' must match the endpoint names.", name))
          value <- value[endpoint_names]
        }
        out[i, ] <- if (length(value) == 1L) rep(value, m) else as.numeric(value)
      }
      return(out)
    }
    if (!is.numeric(x) || !(length(x) %in% c(1L, m)))
      stop(sprintf("'%s' must be scalar, endpoint-specific, or a named arm list.", name))
    if (length(x) > 1L && !is.null(names(x))) {
      if (!identical(names(x), endpoint_names))
        stop(sprintf("Names of '%s' must match the endpoint names.", name))
      x <- x[endpoint_names]
    }
    matrix(if (length(x) == 1L) rep(x, m) else as.numeric(x),
           nrow = length(rates), ncol = m, byrow = TRUE,
           dimnames = list(names(rates), endpoint_names))
  }

  bound_matrix <- function(x, list_x, name) {
    if (!is.null(list_x)) {
      if (!is.list(list_x) || length(list_x) != length(comparisons))
        stop(sprintf("'%s' must contain one vector for each comparison.", name))
      if (is.null(names(list_x))) names(list_x) <- names(comparisons)
      if (!setequal(names(list_x), names(comparisons)))
        stop(sprintf("Names of '%s' must match the comparison names.", name))
      out <- matrix(NA_real_, nrow = length(comparisons), ncol = m,
                    dimnames = list(names(comparisons), endpoint_names))
      for (i in seq_along(comparisons)) {
        value <- list_x[[names(comparisons)[i]]]
        if (!is.numeric(value) || !(length(value) %in% c(1L, m)))
          stop(sprintf("Each '%s' element must have length 1 or %d.", name, m))
        if (length(value) > 1L && !is.null(names(value))) {
          if (!identical(names(value), endpoint_names))
            stop(sprintf("Names of '%s' must match the endpoint names.", name))
          value <- value[endpoint_names]
        }
        out[i, ] <- if (length(value) == 1L) rep(value, m) else value
      }
      return(out)
    }
    if (!is.numeric(x) || !(length(x) %in% c(1L, m)))
      stop(sprintf("'%s' must have length 1 or %d.", name, m))
    if (length(x) > 1L && !is.null(names(x))) {
      if (!identical(names(x), endpoint_names))
        stop(sprintf("Names of '%s' must match the endpoint names.", name))
      x <- x[endpoint_names]
    }
    matrix(if (length(x) == 1L) rep(x, m) else as.numeric(x),
           nrow = length(comparisons), ncol = m, byrow = TRUE,
           dimnames = list(names(comparisons), endpoint_names))
  }
  exposure <- arm_endpoint_matrix(exposure, "exposure")
  dispersion <- arm_endpoint_matrix(dispersion, "dispersion")
  margin_lower <- bound_matrix(margin_lower, list_margin_lower, "margin_lower")
  margin_upper <- bound_matrix(margin_upper, list_margin_upper, "margin_upper")
  alpha <- recycle_endpoint(alpha, "alpha")
  if (is.null(endpoint_corr)) endpoint_corr <- diag(m)
  endpoint_corr <- as.matrix(endpoint_corr)
  if (!is.numeric(endpoint_corr) || !all(dim(endpoint_corr) == c(m, m)))
    stop(sprintf("'endpoint_corr' must be a %d x %d numeric matrix.", m, m))
  if (any(!is.finite(endpoint_corr)) ||
      max(abs(endpoint_corr - t(endpoint_corr))) > 1e-10 ||
      max(abs(diag(endpoint_corr) - 1)) > 1e-10)
    stop("'endpoint_corr' must be finite, symmetric, and have unit diagonal.")
  eigen_values <- eigen(endpoint_corr, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigen_values) <= 0)
    stop("'endpoint_corr' must be positive definite.")

  list(rates = rates_mat, comparisons = comparisons,
       comparison_matrix = comparison_matrix, exposure = exposure,
       dispersion = dispersion,
       margin_lower = margin_lower, margin_upper = margin_upper,
       alpha = alpha, endpoint_corr = endpoint_corr,
       n_arms = nrow(rates_mat), n_endpoints = m,
       n_comparisons = length(comparisons))
}

#' Estimate joint power for correlated count endpoints and multiple comparisons
#'
#' Each simulated trial contains all arms and endpoints. Endpoint counts are
#' generated with their requested marginal Poisson or negative-binomial
#' distributions and a Gaussian-copula dependence structure. Every comparison
#' must pass at least `k` endpoints for the trial to count as a success.
#'
#' @param n_per_arm Subjects in each arm. This joint implementation supports
#' parallel-group designs.
#' @param rates Named list of equal-length endpoint-rate vectors, one per arm.
#' @param comparisons Named list of length-two character vectors. The first arm
#' is the test arm and the second is the reference arm.
#' @param exposure Exposure per subject, scalar or one value per endpoint,
#' or a named list with one scalar/vector per arm.
#' @param margin_lower Lower rate-ratio equivalence margin.
#' @param margin_upper Upper rate-ratio equivalence margin.
#' @param list_margin_lower Optional named list of lower margins, one vector
#' per comparison. Each vector is scalar or has one value per endpoint.
#' @param list_margin_upper Optional named list of upper margins, one vector
#' per comparison. Each vector is scalar or has one value per endpoint.
#' @param model Count model: `"poisson"` or `"negative-binomial"`.
#' @param dispersion Positive negative-binomial dispersion parameter, scalar
#' or a named list with one scalar/vector per arm.
#' @param alpha One-sided significance level, scalar or one value per endpoint.
#' @param endpoint_corr Positive-definite latent Gaussian correlation matrix
#' across endpoints. The default is independence.
#' @param type_y Numeric endpoint hierarchy used with `adjust = "seq"`: `1`
#' for primary/co-primary endpoints and `2` for secondary endpoints.
#' @param k Number of endpoints that must pass within every comparison.
#' @param adjust Multiplicity adjustment within each comparison's selected
#' endpoint family: `"none"`, `"bonferroni"`, `"sidak"`, `"t"`, or
#' `"seq"`/`"sequential"`. The `"t"` option uses Mielke's strong
#' k-out-of-m calibration `alpha / (m - k + 1)`; legacy partial-conjunction
#' labels are accepted. Sequential testing uses the primary gate and
#' secondary-family rule used by the continuous kernels.
#' When all supplied endpoints are required (`k` equals the endpoint count),
#' endpoint-wise adjustment is not necessary for the intersection-union
#' decision; a requested adjustment remains available with a warning.
#' @param nsim Number of simulated trials.
#' @param seed Optional random seed.
#' @param design Joint multi-arm design; currently only `"parallel"` is
#' supported because multiple reference arms do not define a single standard
#' 2x2 crossover design.
#' @return An object of class `countpower` containing joint power and a
#' binomial confidence interval.
#' @examples
#' rates <- list(TEST = c(.20, .20), REF = c(.20, .20), ALT = c(.20, .20))
#' SimTOST:::power_count_joint(100, rates, list(REF = c("TEST", "REF"),
#'                      ALT = c("TEST", "ALT")), nsim = 100, seed = 1)
.power_count_joint_serial <- function(n_per_arm, rates, comparisons,
                              exposure = 1, margin_lower = 0.80,
                              margin_upper = 1.25,
                              model = c("poisson", "negative-binomial"),
                              dispersion = 0.1, alpha = 0.05,
                              endpoint_corr = NULL, k = NULL,
                              type_y = NULL,
                              adjust = c("none", "bonferroni", "sidak", "t", "pc",
                                          "partial-conjunction", "partial_conjunction",
                                          "sequential"),
                              nsim = 5000, seed = NULL,
                              design = c("parallel"),
                              list_margin_lower = NULL,
                              list_margin_upper = NULL,
                              type_y_active = FALSE) {
  model <- match.arg(model)
  adjust <- .normalize_count_adjustment(match.arg(adjust))
  design <- match.arg(design)
  inputs <- prepare_joint_count_inputs(
    rates, comparisons, exposure, margin_lower, margin_upper, alpha,
    endpoint_corr, list_margin_lower, list_margin_upper, dispersion
  )
  k <- .normalize_count_k(k, inputs$n_endpoints)
  if (is.null(k)) k <- inputs$n_endpoints
  if (is.null(type_y) || length(type_y) != inputs$n_endpoints)
    type_y <- rep(-1L, inputs$n_endpoints)
  if (isTRUE(type_y_active) && any(!type_y %in% c(1L, 2L)))
    stop("Active count endpoint hierarchies must contain only 1 or 2.")
  validate_count_inputs(n_per_arm, as.vector(inputs$rates),
                        as.vector(inputs$rates), as.vector(inputs$exposure),
                        inputs$margin_lower, inputs$margin_upper, model,
                        as.vector(inputs$dispersion), inputs$alpha, nsim)
  alpha_endpoint <- .count_endpoint_alpha(
    inputs$alpha, inputs$n_endpoints, k, adjust, type_y, type_y_active
  )
  if (!is.null(seed)) set.seed(seed)
  cpp_result <- count_power_joint_cpp(
    as.integer(n_per_arm), inputs$rates, inputs$exposure,
    inputs$margin_lower, inputs$margin_upper,
    if (model == "poisson") 0L else 1L, inputs$dispersion, alpha_endpoint,
    inputs$endpoint_corr, inputs$comparison_matrix, as.integer(type_y),
    isTRUE(type_y_active), as.integer(k),
    as.integer(nsim)
  )
  ci <- stats::prop.test(cpp_result[["successes"]], nsim, correct = TRUE)$conf.int
  endpoint_successes <- matrix(
    0, nrow = inputs$n_comparisons, ncol = inputs$n_endpoints,
    dimnames = list(names(inputs$comparisons), colnames(inputs$rates))
  )
  for (comparison in seq_len(inputs$n_comparisons)) {
    for (endpoint in seq_len(inputs$n_endpoints)) {
      field <- paste0("comparison_endpoint_success_", comparison, "_", endpoint)
      endpoint_successes[comparison, endpoint] <- cpp_result[[field]]
    }
  }
  comparison_successes <- cpp_result[
    paste0("comparison_success_", seq_len(inputs$n_comparisons))
  ]
  comparison_successes <- unname(as.numeric(comparison_successes))
  names(comparison_successes) <- names(inputs$comparisons)
  out <- list(
    power = cpp_result[["power"]], power_LCI = ci[1], power_UCI = ci[2],
    successes = unname(cpp_result[["successes"]]),
    endpoint_successes = endpoint_successes,
    comparison_successes = comparison_successes,
    n_per_arm = n_per_arm, n_total = n_per_arm * inputs$n_arms,
    n_arms = inputs$n_arms, n_endpoints = inputs$n_endpoints,
    n_comparisons = inputs$n_comparisons, comparisons = inputs$comparisons,
    k = k, model = model, dispersion = inputs$dispersion,
    type_y = if (type_y_active) type_y else NA, adjust = adjust,
    endpoint_corr = inputs$endpoint_corr, nsim = nsim, design = design
  )
  class(out) <- "countpower"
  out
}

#' @rdname power_count_joint
#' @param ncores Number of worker processes used for count simulations.
#' Set to 1 for serial execution. Parallel execution splits `nsim` into
#' reproducible independent chunks and combines the resulting successes.
power_count_joint <- function(n_per_arm, rates, comparisons,
                              exposure = 1, margin_lower = 0.80,
                              margin_upper = 1.25,
                              model = c("poisson", "negative-binomial"),
                              dispersion = 0.1, alpha = 0.05,
                              endpoint_corr = NULL, k = NULL,
                              type_y = NULL,
                              adjust = c("none", "bonferroni", "sidak", "t", "pc",
                                          "partial-conjunction", "partial_conjunction",
                                          "sequential"),
                              nsim = 5000, seed = NULL,
                              design = c("parallel"),
                              list_margin_lower = NULL,
                              list_margin_upper = NULL, ncores = 1,
                              .warn_redundant_bon = TRUE,
                              type_y_active = FALSE) {
  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores != as.integer(ncores) || ncores < 1L)
    stop("'ncores' must be a positive integer.")
  adjust <- .normalize_count_adjustment(match.arg(adjust))
  inputs <- prepare_joint_count_inputs(
    rates, comparisons, exposure, margin_lower, margin_upper, alpha,
    endpoint_corr, list_margin_lower, list_margin_upper, dispersion
  )
  k <- .normalize_count_k(k, inputs$n_endpoints)
  if (is.null(k)) k <- inputs$n_endpoints
  type_info <- .prepare_type_y(
    type_y = type_y, all_endpoints = colnames(inputs$rates),
    selected_endpoints = colnames(inputs$rates),
    adjust = if (identical(adjust, "sequential")) "seq" else "no"
  )
  if (isTRUE(.warn_redundant_bon))
    .warn_adjustment_configuration(
      k = k, m = inputs$n_endpoints,
      adjust = if (identical(adjust, "sequential")) "seq" else adjust,
      type_y = if (type_info$active) type_info$type_y else NULL,
      type_y_supplied = type_info$active,
      n_comparators = inputs$n_comparisons, context = "count endpoints"
    )
  args <- list(n_per_arm = n_per_arm, rates = rates,
               comparisons = comparisons, exposure = exposure,
               margin_lower = margin_lower, margin_upper = margin_upper,
               model = model, dispersion = dispersion, alpha = alpha,
               endpoint_corr = endpoint_corr, k = k,
               type_y = type_info$type_y, type_y_active = type_info$active,
               adjust = adjust,
               nsim = nsim, seed = seed, design = design,
               list_margin_lower = list_margin_lower,
               list_margin_upper = list_margin_upper)
  out <- .combine_count_power_chunks(.power_count_joint_serial,
                                     ".power_count_joint_serial", args, nsim,
                                     as.integer(ncores), seed)
  out$parameters <- c(args, list(ncores = ncores,
                                 .function = "power_count_joint"))
  out
}

#' Estimate sample size for joint correlated count equivalence
#'
#' @inheritParams power_count_joint
#' @param power Target joint power.
#' @param type_y Numeric endpoint hierarchy used with `adjust = "seq"`: `1`
#' for primary/co-primary endpoints and `2` for secondary endpoints.
#' @param lower Minimum subjects per arm.
#' @param upper Maximum subjects per arm.
#' @param optimization_method Search method: `"fast"` uses bracketing and
#' integer bisection; `"step-by-step"` evaluates every candidate.
#' @param step.power Initial power-of-two jump for the fast search.
#' @param step.up Direction of the initial fast-search bracketing.
#' @param pos.side Retained for compatibility with `sampleSize()`; count
#' searches always return the smallest candidate reaching the target.
#' @param maxiter Maximum number of power evaluations.
#' @return An object of class `countss` containing the selected sample
#' size, achieved joint power, confidence interval, input parameters, and the
#' search history in `table.iter` and `table.test`. For count outcomes,
#' `table.iter` has one row per evaluated candidate sample size and
#' `table.test` contains complete-trial, comparator, and endpoint decision
#' indicators for each simulated trial and candidate. The count kernel returns
#' aggregate decision counts rather than raw endpoint-level test statistics, so
#' component columns preserve the simulated marginal success counts.
sampleSize_count_joint <- function(power = 0.80, rates, comparisons,
                                   exposure = 1, margin_lower = 0.80,
                                   margin_upper = 1.25,
                                   model = c("poisson", "negative-binomial"),
                                   dispersion = 0.1, alpha = 0.05,
                                   endpoint_corr = NULL, k = NULL,
                                   type_y = NULL,
                                   adjust = c("none", "bonferroni", "sidak", "t", "pc",
                                               "partial-conjunction", "partial_conjunction",
                                               "sequential"),
                                   nsim = 5000, seed = NULL, lower = 2,
                                   upper = 500, design = c("parallel"),
                                   list_margin_lower = NULL,
                                   list_margin_upper = NULL,
                                   optimization_method = c("fast", "step-by-step"),
                                   step.power = 6, step.up = TRUE,
                                   pos.side = FALSE, maxiter = 1000, ncores = 1,
                                   .warn_redundant_bon = TRUE) {
  design <- match.arg(design)
  optimization_method <- match.arg(optimization_method)
  if (!is.numeric(power) || length(power) != 1L || power <= 0 || power >= 1)
    stop("'power' must be a single number between 0 and 1.")
  if (!is.numeric(lower) || !is.numeric(upper) || lower > upper ||
      lower != as.integer(lower) || upper != as.integer(upper) || lower < 2)
    stop("'lower' and 'upper' must be integers with lower >= 2.")
  adjust <- .normalize_count_adjustment(match.arg(adjust))
  if (!is.list(rates) || !length(rates))
    stop("'rates' must be a non-empty list.")
  m <- max(vapply(rates, length, integer(1)))
  joint_endpoint_names <- names(rates[[1L]])
  if (is.null(joint_endpoint_names))
    joint_endpoint_names <- paste0("endpoint_", seq_len(m))
  k <- .normalize_count_k(k, m)
  if (is.null(k)) k <- m
  type_info <- .prepare_type_y(
    type_y = type_y, all_endpoints = joint_endpoint_names,
    selected_endpoints = joint_endpoint_names,
    adjust = if (identical(adjust, "sequential")) "seq" else "no"
  )
  if (isTRUE(.warn_redundant_bon))
    .warn_adjustment_configuration(
      k = k, m = m,
      adjust = if (identical(adjust, "sequential")) "seq" else adjust,
      type_y = if (type_info$active) type_info$type_y else NULL,
      type_y_supplied = type_info$active,
      n_comparators = length(comparisons), context = "count endpoints"
    )
  power_fun <- function(n) {
    power_count_joint(
      n_per_arm = n, rates = rates, comparisons = comparisons,
      exposure = exposure, margin_lower = margin_lower,
      margin_upper = margin_upper, model = model, dispersion = dispersion,
      alpha = alpha, endpoint_corr = endpoint_corr, k = k, adjust = adjust,
      type_y = type_y,
      nsim = nsim, seed = seed, design = design,
      list_margin_lower = list_margin_lower,
      list_margin_upper = list_margin_upper, ncores = ncores,
      .warn_redundant_bon = FALSE
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
  out <- c(achieved, list(
    target_power = power,
    parameters = list(rates = rates, comparisons = comparisons,
                      exposure = exposure, margin_lower = margin_lower,
                      margin_upper = margin_upper, model = model,
                      dispersion = dispersion, alpha = alpha,
                      endpoint_corr = endpoint_corr, k = achieved$k,
                      type_y = type_y,
                      adjust = adjust, lower = lower, upper = upper,
                      list_margin_lower = list_margin_lower,
                      list_margin_upper = list_margin_upper,
                      design = design,
                      optimization_method = optimization_method,
                      step.power = step.power, step.up = step.up,
                      pos.side = pos.side, maxiter = maxiter)
  ))
  class(out) <- "countss"
  out$parameters <- list(
    rates = rates, comparisons = comparisons, exposure = exposure,
    margin_lower = margin_lower, margin_upper = margin_upper,
    model = model, dispersion = dispersion, alpha = alpha, nsim = nsim,
    seed = seed, design = design, k = out$k,
    endpoint_corr = endpoint_corr, type_y = type_y, adjust = adjust,
    list_margin_lower = list_margin_lower,
    list_margin_upper = list_margin_upper, power = power,
    lower = lower, upper = upper, ncores = ncores,
    optimization_method = optimization_method, step.power = step.power,
    step.up = step.up, pos.side = pos.side, maxiter = maxiter,
    .function = "sampleSize_count_joint"
  )
  out$table.iter <- history$table.iter
  out$table.test <- history$table.test
  out
}
