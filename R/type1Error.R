#' Empirical Type I Error at the Least-Favorable Null
#'
#' Generates a boundary-null configuration and estimates its empirical
#' rejection probability using [simPower()]. The comparator convention is
#' `c(test, reference)`. For normal ROM and count RR, a lower-bound null
#' sets `mean(test) / mean(reference) = L`. For log-normal ROM, the boundary
#' is imposed on the log-analysis scale used by the test; this coincides with
#' the arithmetic mean ratio when the two arms have the same coefficient of
#' variation.
#'
#' @param null Boundary to evaluate: `"lower"`, `"upper"`, or `"both"`.
#' @param x Optional existing `simss` or `simpower` object. When supplied, its
#'   design and outcome parameters are reused; arguments in `...` override
#'   stored settings such as `n` or `nsim`.
#' @param comparator Optional comparator name identifying the boundary
#'   scenario when `joint = FALSE`. It must be one of the names in
#'   `list_comparator`; by default, the first comparator is used. It is
#'   ignored when `joint = TRUE`, because joint mode evaluates all comparators.
#' @param endpoint Optional endpoint identifying the boundary component when
#'   `joint = FALSE`. It is available for the all-endpoints-required case
#'   (`k` equal to the number of endpoints), where one endpoint is placed on
#'   the boundary. For `k` smaller than the number of endpoints, use
#'   `joint = TRUE` so that all composite boundary configurations are evaluated.
#' @param joint Logical. If `TRUE`, evaluate every comparator and every
#'   endpoint partial-null configuration and return the joint Type I-error
#'   result for each complete-trial decision. For a comparator with `m`
#'   endpoints and an at-least-`k` rule, all boundary counts from
#'   `m - k + 1` through `m` are evaluated. Defaults to `FALSE` for backward
#'   compatibility.
#' @param conf.level Confidence level for the simultaneous one-sided Monte
#'   Carlo upper bound across the evaluated joint scenarios. Defaults to
#'   `0.95`.
#' @details Comparators use the convention `c(test, reference)`. For ROM,
#'   the tested estimand is `test / reference`; for DOM it is
#'   `test - reference`. With log-normal ROM, the test is performed after
#'   converting the supplied original-scale means and variances to the
#'   log-analysis scale, so the boundary scenario is calibrated on that
#'   scale. For count rate ratios, the analogous midpoint is `sqrt(L * U)` on
#'   the rate-ratio scale because the test is performed on the log-rate-ratio
#'   scale. Absolute rates and dispersion remain nuisance parameters, so this
#'   midpoint should be supplemented by a grid or optimization when a global
#'   supremum is required. If a comparator has `m` endpoints and requires `k`
#'   endpoints to pass, a composite null configuration has at least `m - k + 1`
#'   non-equivalent endpoints. With `joint = TRUE`, all endpoint subsets with
#'   boundary counts from `m - k + 1` through `m`, and all lower/upper
#'   direction combinations, are evaluated, while the complete decision still
#'   requires the `k`-of-`m` rule for every comparator in the same simulated
#'   trial. The returned joint object includes a Bonferroni simultaneous
#'   one-sided Monte Carlo upper bound for the maximum scenario probability.
#' @param ... Arguments passed to [simPower()]. The call must include the
#'   outcome parameters, comparator definitions, and equivalence margins.
#' @return For `null = "lower"` or `"upper"`, an object of class
#'   `type1error` containing the empirical Type I error in `type1_error` and
#'   its Monte Carlo interval. For `null = "both"`, a named list containing
#'   the lower- and upper-bound results. With `joint = TRUE`, a
#'   `type1error_joint` object containing one joint result for every valid
#'   comparator, endpoint-subset, boundary-direction combination is returned.
#' @export
type1Error <- function(null = c("lower", "upper", "both"), x = NULL,
                       comparator = NULL, endpoint = NULL, joint = FALSE,
                       conf.level = 0.95, ...) {
  null_missing <- missing(null)
  null <- match.arg(null)
  if (!is.logical(joint) || length(joint) != 1L || is.na(joint))
    stop("'joint' must be TRUE or FALSE.")
  if (joint && null_missing) null <- "both"
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number strictly between 0 and 1.")
  args <- list(...)
  if (!is.null(x)) {
    if (!(inherits(x, "simss") || inherits(x, "simpower")))
      stop("'x' must be a simss or simpower object.")
    args <- utils::modifyList(.type1_args_from_object(x), args)
  }
  required <- c("n", "list_comparator")
  missing_args <- required[vapply(required, function(z) is.null(args[[z]]), logical(1))]
  if (length(missing_args))
    stop("type1Error() requires: ", paste(missing_args, collapse = ", "))

  distribution <- .normalize_distribution(
    if (is.null(args$distribution)) "norm" else args$distribution
  )
  is_count <- distribution %in% c("pois", "nbinom")
  ctype <- toupper(if (is.null(args$ctype)) if (is_count) "RR" else "ROM" else args$ctype)
  if ((!is_count && !ctype %in% c("DOM", "ROM")) ||
      (is_count && ctype != "RR"))
    stop("The supplied ctype is not available for this distribution.")

  comparators <- args$list_comparator
  if (!is.list(comparators) || !length(comparators))
    stop("'list_comparator' must be a non-empty list.")
  comparator_names <- names(comparators)
  if (is.null(comparator_names)) {
    comparator_names <- paste0("comparison_", seq_along(comparators))
    names(comparators) <- comparator_names
  }
  if (any(vapply(comparators, length, integer(1)) != 2L))
    stop("Each comparator must contain exactly two arms: c(test, reference).")
  if (joint || is.null(comparator)) {
    comparator_index <- 1L
  } else {
    if (!is.character(comparator) || length(comparator) != 1L ||
        is.na(comparator) || !comparator %in% comparator_names)
      stop("'comparator' must be one of the names in 'list_comparator'.")
    comparator_index <- match(comparator, comparator_names)
  }

  arm_parameters <- if (is_count) args$rate_list else args$mu_list
  arm_endpoints <- lapply(arm_parameters, function(value) {
    if (is.matrix(value) || is.data.frame(value)) colnames(value) else
      names(value)
  })
  arm_endpoints <- lapply(seq_along(arm_parameters), function(i) {
    value <- arm_endpoints[[i]]
    if (!is.null(value)) value else {
      n <- length(arm_parameters[[i]])
      prefix <- if (is_count) "endpoint_" else "y"
      paste0(prefix, seq_len(n))
    }
  })
  names(arm_endpoints) <- names(arm_parameters)
  requested_list_y_comparator <- args$list_y_comparator
  endpoint_sets <- .resolve_comparator_endpoints(
    comparators = comparators, arm_endpoints = arm_endpoints,
    requested = args$list_y_comparator
  )
  .warn_inferred_endpoint_reduction(
    comparators = comparators, arm_endpoints = arm_endpoints,
    endpoint_sets = endpoint_sets, requested = requested_list_y_comparator,
    context = "list_y_comparator"
  )
  endpoint_names <- function(i, arms, lower_i) endpoint_sets[[i]]
  first_endpoints <- endpoint_names(comparator_index,
                                    comparators[[comparator_index]], NULL)
  if (is.null(endpoint)) endpoint <- first_endpoints[[1L]]
  if (length(endpoint) != 1L || !is.character(endpoint) || is.na(endpoint) ||
      !endpoint %in% first_endpoints)
    stop("'endpoint' must identify one endpoint of the selected comparator.")

  lower <- args$list_lequi.tol
  upper <- args$list_uequi.tol
  if (!is.list(lower)) {
    if (is.null(args$lequi.tol) || all(is.na(args$lequi.tol)))
      stop("Supply 'list_lequi.tol' or 'lequi.tol'.")
    lower <- rep(list(args$lequi.tol), length(comparators))
  }
  if (!is.list(upper)) {
    if (is.null(args$uequi.tol) || all(is.na(args$uequi.tol)))
      stop("Supply 'list_uequi.tol' or 'uequi.tol'.")
    upper <- rep(list(args$uequi.tol), length(comparators))
  }
  if (is.null(names(lower))) names(lower) <- comparator_names
  if (is.null(names(upper))) names(upper) <- comparator_names
  lower <- lower[comparator_names]
  upper <- upper[comparator_names]

  k_values <- args$k
  k_scalar_input <- !is.list(k_values) && length(k_values) == 1L &&
    !is.null(k_values) && !is.na(k_values)
  if (is.list(k_values)) k_values <- unlist(k_values, use.names = FALSE)
  if (is.null(k_values) || !length(k_values) || all(is.na(k_values))) {
    k_values <- vapply(endpoint_sets, length, integer(1))
  } else if (!is.numeric(k_values)) {
    stop("'k' must be numeric.")
  } else if (length(k_values) == 1L) {
    k_values <- rep(k_values, length(comparators))
  } else if (length(k_values) != length(comparators)) {
    stop("'k' must have one value per comparator or be a single value.")
  }
  if (any(!is.finite(k_values)) || any(k_values != as.integer(k_values)) ||
      any(k_values < 1L))
    stop("'k' must contain positive integers.")
  endpoint_counts <- lengths(endpoint_sets)
  oversized <- which(k_values > endpoint_counts)
  if (length(oversized)) {
    warning("'k' is larger than the number of selected endpoints for comparator(s) ",
            paste(oversized, collapse = ", "), "; setting it to the maximum possible value.",
            call. = FALSE)
    k_values[oversized] <- endpoint_counts[oversized]
  }
  k_values <- as.integer(k_values)
  if (k_scalar_input && length(unique(endpoint_counts)) > 1L) {
    warning("The selected endpoint counts differ by comparator; the supplied scalar 'k' is interpreted separately for each comparator and capped at that comparator's m.",
            call. = FALSE)
  }
  .warn_adjustment_configuration(
    k = k_values, m = endpoint_counts,
    adjust = if (is.null(args$adjust)) "no" else args$adjust,
    n_comparators = length(comparators), context = "selected endpoints"
  )
  boundary_counts <- vapply(seq_along(endpoint_sets), function(i) {
    length(endpoint_sets[[i]]) - k_values[[i]] + 1L
  }, integer(1))

  if (!joint && any(boundary_counts != 1L)) {
    stop("For k smaller than the number of endpoints, use joint = TRUE to evaluate all composite boundary configurations.")
  }

  boundary_label <- function(directions) {
    if (all(directions == "lower")) return("lower")
    if (all(directions == "upper")) return("upper")
    paste(directions, collapse = "/")
  }

  make_scenarios <- function(indices) {
    scenarios <- list()
    index <- 0L
    for (i in indices) {
      boundary_sizes <- if (joint) {
        seq.int(boundary_counts[[i]], length(endpoint_sets[[i]]))
      } else {
        1L
      }
      for (boundary_size in boundary_sizes) {
        subsets <- utils::combn(endpoint_sets[[i]], boundary_size,
                                simplify = FALSE)
        direction_sets <- if (null == "both") {
          grid <- expand.grid(
            rep(list(c("lower", "upper")), boundary_size),
            stringsAsFactors = FALSE
          )
          lapply(seq_len(nrow(grid)), function(j) as.character(grid[j, ]))
        } else {
          list(rep(null, boundary_size))
        }
        for (selected_endpoints in subsets) {
          for (directions in direction_sets) {
            if (!joint && !identical(selected_endpoints[[1L]], endpoint)) next
            index <- index + 1L
            scenarios[[index]] <- list(
              comparator_index = i,
              endpoints = selected_endpoints,
              directions = directions,
              Boundary = boundary_label(directions),
              NullCount = boundary_size,
              Endpoint = paste(selected_endpoints, collapse = "+")
            )
          }
        }
      }
    }
    scenarios
  }

  # The log-normal ROM implementation converts the supplied original-scale
  # means and variances to log-scale means before applying the DOM kernel. A
  # boundary therefore has to be imposed on that analysis scale. Simply
  # setting arithmetic_mean(test) = arithmetic_mean(reference) * margin is
  # only equivalent when the two arms have the same coefficient of variation.
  endpoint_variance <- function(arm, ep) {
    vc <- args$varcov_list
    if (is.list(vc) && !is.null(vc[[arm]])) {
      value <- vc[[arm]]
      if (is.matrix(value) || is.data.frame(value)) {
        if (!is.null(rownames(value)) && ep %in% rownames(value) &&
            !is.null(colnames(value)) && ep %in% colnames(value))
          return(as.numeric(value[ep, ep]))
        j <- match(ep, colnames(value))
        if (!is.na(j)) return(as.numeric(value[j, j]))
      }
    }
    sig <- args$sigma_list
    if (is.list(sig) && !is.null(sig[[arm]])) {
      value <- .endpoint_parameter(sig[[arm]], ep)
      if (length(value) == 1L && is.finite(value)) return(value^2)
    }
    stop("For log-normal ROM boundaries, supply arm-specific variances " ,
         "through 'varcov_list' or standard deviations through 'sigma_list'.")
  }

  lognormal_mean <- function(mean, variance) {
    if (length(mean) != 1L || !is.finite(mean) || mean <= 0 ||
        length(variance) != 1L || !is.finite(variance) || variance < 0)
      stop("Log-normal ROM means must be positive and variances non-negative.")
    log(mean) - 0.5 * log1p(variance / mean^2)
  }

  arithmetic_mean_at_lognormal_boundary <- function(reference, margin,
                                                     reference_variance,
                                                     target_variance) {
    if (!is.finite(margin) || margin <= 0)
      stop("Log-normal ROM equivalence margins must be positive.")
    target_log_mean <- lognormal_mean(reference, reference_variance) +
      log(margin)
    q <- exp(2 * target_log_mean)
    sqrt((q + sqrt(q^2 + 4 * q * target_variance)) / 2)
  }

  analysis_estimand <- function(target, reference, target_arm, reference_arm,
                                ep) {
    if (distribution == "lnorm" && ctype == "ROM") {
      exp(lognormal_mean(target, endpoint_variance(target_arm, ep)) -
            lognormal_mean(reference, endpoint_variance(reference_arm, ep)))
    } else if (ctype == "DOM") {
      target - reference
    } else {
      target / reference
    }
  }

  make_boundary_call <- function(target_comparator, target_endpoints,
                                 target_directions) {
    boundary_args <- args
    if (is_count) {
      if (is.null(boundary_args$rate_list))
        stop("'rate_list' is required for count outcomes.")
      parameters <- boundary_args$rate_list
      parameter_name <- "rate_list"
    } else {
      if (is.null(boundary_args$mu_list))
        stop("'mu_list' is required for continuous outcomes.")
      parameters <- boundary_args$mu_list
      parameter_name <- "mu_list"
    }

    set_endpoint <- function(values, arm, ep, target, position) {
      value <- values[[arm]]
      if (is.matrix(value) || is.data.frame(value)) {
        col <- match(ep, colnames(value))
        if (is.na(col)) stop("Endpoint names do not match '", parameter_name, "'.")
        value[, col] <- target
      } else if (!is.null(names(value))) {
        value[[ep]] <- target
      } else {
        value[[position]] <- target
      }
      values[[arm]] <- value
      values
    }

    margin_for <- function(margins, ep, position) {
      if (length(margins) == 1L) margins[[1L]] else
        if (!is.null(names(margins))) margins[[ep]] else margins[[position]]
    }

    analysis_midpoint <- function(lower_margin, upper_margin) {
      if (ctype == "DOM") (lower_margin + upper_margin) / 2 else
        sqrt(lower_margin * upper_margin)
    }

    target_from_reference <- function(reference, margin, test_arm,
                                      reference_arm, ep) {
      if (ctype == "DOM") {
        reference + margin
      } else if (distribution == "lnorm") {
        arithmetic_mean_at_lognormal_boundary(
          reference = reference, margin = margin,
          reference_variance = endpoint_variance(reference_arm, ep),
          target_variance = endpoint_variance(test_arm, ep)
        )
      } else {
        reference * margin
      }
    }

    reference_from_target <- function(target, margin, test_arm,
                                      reference_arm, ep) {
      if (ctype == "DOM") {
        target - margin
      } else if (distribution == "lnorm") {
        target_log_mean <- lognormal_mean(
          target, endpoint_variance(test_arm, ep)
        ) - log(margin)
        q <- exp(2 * target_log_mean)
        sqrt((q + sqrt(q^2 + 4 * q * endpoint_variance(reference_arm, ep))) / 2)
      } else {
        target / margin
      }
    }

    # First place every non-boundary component at the midpoint of its
    # user-supplied equivalence interval on the analysis scale. For log-normal
    # ROM, target_from_reference() converts this midpoint back to an
    # arithmetic mean with the same variance calibration used for the selected
    # boundary. For count models, the midpoint fixes the rate ratio but does
    # not fix the nuisance-rate or dispersion configuration.
    for (i in seq_along(comparators)) {
      arms <- comparators[[i]]
      endpoints <- endpoint_names(i, arms, lower[[i]])
      lower_i <- lower[[i]]
      upper_i <- upper[[i]]
      for (j in seq_along(endpoints)) {
        ep <- endpoints[[j]]
        center <- analysis_midpoint(
          margin_for(lower_i, ep, j), margin_for(upper_i, ep, j)
        )
        reference <- .endpoint_parameter(boundary_args[[parameter_name]][[arms[[2L]]]], ep)
        target <- target_from_reference(reference, center, arms[[1L]], arms[[2L]], ep)
        parameters <- set_endpoint(parameters, arms[[1L]], ep, target, j)
      }
    }

    selected_arms <- comparators[[target_comparator]]
    selected_all_endpoints <- endpoint_names(
      target_comparator, selected_arms, lower[[target_comparator]]
    )
    selected_positions <- match(target_endpoints, selected_all_endpoints)
    fixed_keys <- paste(selected_arms[[1L]], target_endpoints, sep = "::")
    for (j in seq_along(target_endpoints)) {
      target_ep <- target_endpoints[[j]]
      target_position <- selected_positions[[j]]
      selected_margins <- if (target_directions[[j]] == "lower")
        lower[[target_comparator]] else upper[[target_comparator]]
      selected_margin <- margin_for(selected_margins, target_ep, target_position)
      selected_reference <- .endpoint_parameter(
        boundary_args[[parameter_name]][[selected_arms[[2L]]]], target_ep
      )
      selected_target <- target_from_reference(
        selected_reference, selected_margin, selected_arms[[1L]],
        selected_arms[[2L]], target_ep
      )
      parameters <- set_endpoint(
        parameters, selected_arms[[1L]], target_ep,
        selected_target, target_position
      )
    }

    # Propagate the selected boundary through shared arms while keeping every
    # other comparator--endpoint component at the interval centre.
    for (i in seq_along(comparators)) {
      arms <- comparators[[i]]
      endpoints <- endpoint_names(i, arms, lower[[i]])
      lower_i <- lower[[i]]
      upper_i <- upper[[i]]
      for (j in seq_along(endpoints)) {
        ep <- endpoints[[j]]
        key <- paste(arms[[1L]], ep, sep = "::")
        if (i == target_comparator && ep %in% target_endpoints) next
        center <- analysis_midpoint(
          margin_for(lower_i, ep, j), margin_for(upper_i, ep, j)
        )
        reference <- .endpoint_parameter(parameters[[arms[[2L]]]], ep)
        if (key %in% fixed_keys) {
          reference <- reference_from_target(
            .endpoint_parameter(parameters[[arms[[1L]]]], ep), center,
            arms[[1L]], arms[[2L]], ep
          )
          parameters <- set_endpoint(parameters, arms[[2L]], ep,
                                      reference, j)
        } else {
          target <- target_from_reference(
            reference, center, arms[[1L]], arms[[2L]], ep
          )
          parameters <- set_endpoint(parameters, arms[[1L]], ep,
                                      target, j)
        }
      }
    }
    boundary_args[[parameter_name]] <- parameters
    boundary_args$ctype <- ctype
    boundary_args$distribution <- distribution
    boundary_args
  }

  boundary_estimands <- function(boundary_args, target_comparator,
                                 target_endpoints) {
    parameter_values <- if (is_count) boundary_args$rate_list else
      boundary_args$mu_list
    rows <- list()
    index <- 0L
    for (i in seq_along(comparators)) {
      arms <- comparators[[i]]
      endpoints <- endpoint_names(i, arms, lower[[i]])
      for (ep in endpoints) {
        reference <- .endpoint_parameter(parameter_values[[arms[[2L]]]], ep)
        target <- .endpoint_parameter(parameter_values[[arms[[1L]]]], ep)
        if (i != target_comparator || !ep %in% target_endpoints) next
        index <- index + 1L
        rows[[index]] <- data.frame(
          Comparator = comparator_names[[i]], Endpoint = ep,
          Estimand = analysis_estimand(
            target = target, reference = reference,
            target_arm = arms[[1L]], reference_arm = arms[[2L]], ep = ep
          ),
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  }

  run_one <- function(scenario) {
    target_comparator <- scenario$comparator_index
    target_endpoints <- scenario$endpoints
    boundary_args <- make_boundary_call(
      target_comparator = target_comparator,
      target_endpoints = target_endpoints,
      target_directions = scenario$directions
    )
    # The public type1Error() call warns once for the requested adjustment;
    # boundary scenarios should not repeat that same warning.
    boundary_args$.warn_redundant_bon <- FALSE
    result <- do.call(simPower, boundary_args)
    result$null_boundary <- scenario$Boundary
    result$boundary_comparator <- comparator_names[[target_comparator]]
    result$boundary_endpoints <- target_endpoints
    result$boundary_endpoint <- scenario$Endpoint
    result$boundary_directions <- scenario$directions
    result$boundary_null_count <- scenario$NullCount
    result$boundary_estimands <- boundary_estimands(
      boundary_args, target_comparator, target_endpoints
    )
    test_table <- if (!is.null(result$result)) result$result$table.test else NULL
    result$mc_nsim <- as.integer(if (!is.null(result$nsim)) result$nsim else
      if (!is.null(args$nsim)) args$nsim else 5000L)
    total_col <- if (!is.null(test_table)) {
      if ("totaly" %in% names(test_table)) "totaly" else
        grep("^totaly", names(test_table), value = TRUE)[1L]
    } else {
      NA_character_
    }
    result$mc_successes <- if (!is.na(total_col)) {
      sum(test_table[[total_col]], na.rm = TRUE)
    } else if (!is.null(result$successes)) {
      as.integer(result$successes)
    } else {
      as.integer(round(result$power * result$mc_nsim))
    }
    # simPower() receives k from the source simss object through
    # .type1_args_from_object(). Its decision engine applies that k value for
    # every comparator, and its returned power and Monte Carlo interval are
    # therefore the Type I-error estimate and interval for this scenario.
    result$type1_error <- result$power
    comparator_label <- paste(comparators[[target_comparator]], collapse = " vs ")
    component_col <- paste0(target_endpoints[[1L]], "Comp:", comparator_label)
    if (!is.null(test_table) && component_col %in% names(test_table)) {
      successes <- sum(test_table[[component_col]], na.rm = TRUE)
      trials <- sum(!is.na(test_table[[component_col]]))
      result$component_type1_error <- successes / trials
      component_ci <- stats::prop.test(successes, trials, correct = TRUE)$conf.int
      result$component_LCI <- unname(component_ci[[1L]])
      result$component_UCI <- unname(component_ci[[2L]])
    } else {
      result$component_type1_error <- NA_real_
      result$component_LCI <- NA_real_
      result$component_UCI <- NA_real_
    }
    class(result) <- c("type1error", class(result))
    result
  }
  if (joint) {
    scenarios <- make_scenarios(seq_along(comparators))
    results <- list()
    rows <- list()
    index <- 0L
    for (scenario in scenarios) {
      result <- run_one(scenario)
      key <- paste(scenario$Boundary,
                   comparator_names[[scenario$comparator_index]],
                   scenario$Endpoint, sep = "__")
      results[[key]] <- result
      index <- index + 1L
      rows[[index]] <- data.frame(
        Boundary = scenario$Boundary,
        Comparator = comparator_names[[scenario$comparator_index]],
        Endpoint = scenario$Endpoint,
        NullCount = scenario$NullCount,
        Type1_Error = result$type1_error,
        Lower = result$power_LCI,
        Upper = result$power_UCI,
        Simultaneous_Upper = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    table <- do.call(rbind, rows)
    n_scenarios <- nrow(table)
    table$Simultaneous_Upper <- vapply(seq_len(n_scenarios), function(i) {
      result <- results[[paste(
        table$Boundary[[i]], table$Comparator[[i]], table$Endpoint[[i]],
        sep = "__"
      )]]
      .type1_simultaneous_upper(
        result$mc_successes, result$mc_nsim, conf.level, n_scenarios
      )
    }, numeric(1))
    global_row <- table[which.max(table$Type1_Error), , drop = FALSE]
    out <- list(
      results = results, table = table, global = global_row,
      alpha = if (is.null(args$alpha)) 0.05 else as.numeric(args$alpha[[1L]]),
      null = null,
      nsim = if (is.null(args$nsim)) NA_real_ else as.numeric(args$nsim[[1L]]),
      conf.level = conf.level,
      global_upper = max(table$Simultaneous_Upper, na.rm = TRUE),
      simultaneous_method = "Bonferroni one-sided Clopper-Pearson upper bound"
    )
    for (i in seq_along(results)) {
      key <- names(results)[[i]]
      results[[i]]$simultaneous_upper <-
        table$Simultaneous_Upper[[match(key, names(results))]]
    }
    out$results <- results
    class(out) <- "type1error_joint"
    return(out)
  }
  if (null == "both") {
    scenarios <- make_scenarios(comparator_index)
    out <- lapply(scenarios, run_one)
    n_scenarios <- length(out)
    for (i in seq_along(out)) {
      out[[i]]$simultaneous_upper <- .type1_simultaneous_upper(
        out[[i]]$mc_successes, out[[i]]$mc_nsim, conf.level, n_scenarios
      )
    }
    names(out) <- vapply(scenarios, function(scenario) scenario$Boundary,
                         character(1))
    attr(out, "conf.level") <- conf.level
    attr(out, "global_upper") <- max(vapply(
      out, function(z) z$simultaneous_upper, numeric(1)
    ), na.rm = TRUE)
    class(out) <- "type1error_set"
    return(out)
  }
  scenarios <- make_scenarios(comparator_index)
  result <- run_one(scenarios[[1L]])
  result$simultaneous_upper <- .type1_simultaneous_upper(
    result$mc_successes, result$mc_nsim, conf.level, length(scenarios)
  )
  result
}

.type1_simultaneous_upper <- function(successes, nsim, conf.level,
                                      n_scenarios) {
  if (!is.finite(successes) || !is.finite(nsim) || nsim < 1L ||
      successes < 0 || successes > nsim || n_scenarios < 1L)
    return(NA_real_)
  tail_alpha <- (1 - conf.level) / n_scenarios
  if (successes >= nsim) return(1)
  stats::qbeta(1 - tail_alpha, successes + 1, nsim - successes)
}

.type1_args_from_object <- function(x) {
  # Count objects retain their original count inputs explicitly. This includes
  # unified simPower() objects, which inherit from countpower but not countss.
  if ((inherits(x, "countss") || inherits(x, "countpower")) &&
      !is.null(x$parameters)) {
    p <- x$parameters
    if (!is.null(p$rate_list)) {
      return(list(
        n = if (!is.null(x$n_per_arm)) x$n_per_arm else x$n,
        distribution = p$distribution, rate_list = p$rate_list,
        list_comparator = p$list_comparator,
        list_y_comparator = p$list_y_comparator,
        list_lequi.tol = p$list_lequi.tol,
        list_uequi.tol = p$list_uequi.tol,
        exposure = p$exposure, dispersion = p$dispersion,
        nsim = p$nsim, seed = NULL, dtype = p$dtype, ctype = "RR",
        type_y = if (!is.null(p$type_y)) p$type_y else NA,
        k = if (!is.null(x$k)) x$k else p$k,
        adjust = if (!is.null(p$adjust)) p$adjust else "no"
      ))
    }
    distribution <- if (identical(p$model, "poisson")) "pois" else "nbinom"
    comparisons <- if (!is.null(p$comparisons)) p$comparisons else
      list(comparison_1 = c("test", "reference"))
    rates <- if (!is.null(p$rates)) p$rates else
      list(test = p$rate_test, reference = p$rate_reference)
    return(list(
      n = x$n_per_arm, distribution = distribution, rate_list = rates,
      list_comparator = comparisons,
      list_lequi.tol = if (!is.null(p$list_margin_lower)) p$list_margin_lower else
        list(comparison_1 = p$margin_lower),
      list_uequi.tol = if (!is.null(p$list_margin_upper)) p$list_margin_upper else
        list(comparison_1 = p$margin_upper),
      exposure = p$exposure, dispersion = p$dispersion,
      nsim = x$nsim, seed = NULL, dtype = p$design, ctype = "RR",
      type_y = if (!is.null(p$type_y)) p$type_y else NA,
      k = if (!is.null(x$k)) x$k else p$k,
      adjust = if (!is.null(p$adjust)) p$adjust else "none"
    ))
  }

  result <- if (inherits(x, "simpower")) x$result else x
  param <- result$param
  param.d <- result$param.d
  if (is.null(param) || is.null(param.d))
    stop("The supplied object does not retain the parameters needed for type1Error().")
  mu_list <- lapply(seq_along(param$mu), function(i) {
    value <- as.numeric(param$mu[[i]])
    endpoint_names <- param$ynames_list[[param$arm_names[[i]]]]
    if (!is.null(endpoint_names)) names(value) <- endpoint_names
    value
  })
  names(mu_list) <- names(param$mu)
  list(
    n = if (inherits(x, "simpower")) x$n else as.integer(result$response$n_iter[[1L]]),
    distribution = param.d$distribution,
    mu_list = mu_list,
    varcov_list = param$varcov,
    sigmaB = param$sigmaB,
    Eper = param$Eper, Eco = param$Eco,
    TAR = unlist(param$TAR_list), arm_names = param$arm_names,
    ynames_list = param$ynames_list, type_y = param$type_y,
    list_comparator = param$list_comparator,
    list_y_comparator = param$list_y_comparator,
    list_lequi.tol = param.d$list_lequi.tol,
    list_uequi.tol = param.d$list_uequi.tol,
    alpha = param.d$alpha, dtype = param.d$dtype, ctype = param.d$ctype,
    vareq = param.d$vareq, k = param.d$k, adjust = param.d$adjust,
    dropout = param.d$dropout, nsim = param.d$nsim, seed = 1234,
    ncores = 1
  )
}

#' @export
print.type1error <- function(x, ...) {
  if (inherits(x, "simpower_curve")) {
    cat("Empirical Type I error curve (", x$null_boundary, " boundary)\n", sep = "")
    print(x$power_curve, row.names = FALSE)
  } else {
    cat("Empirical Type I error (", x$null_boundary, " boundary)\n", sep = "")
    if (!is.null(x$boundary_estimands)) {
      values <- x$boundary_estimands$Estimand
      cat(sprintf("Boundary comparator: %s\n", x$boundary_comparator))
      cat(sprintf("Boundary endpoint(s): %s\n",
                  paste(x$boundary_endpoints, collapse = ", ")))
      cat(sprintf("Null estimand: %s\n",
                  paste(format(values, digits = 4), collapse = ", ")))
    }
    cat(sprintf("Global/joint Type I error: %.4f [%.4f, %.4f]\n",
                x$type1_error, x$power_LCI, x$power_UCI))
  }
  invisible(x)
}

#' @export
print.type1error_set <- function(x, ...) {
  cat("Empirical Type I error at both equivalence boundaries\n")
  lapply(x, print)
  invisible(x)
}
