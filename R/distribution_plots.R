.distribution_sim_data <- function(x, display = "all", endpoint = NULL,
                                    arm = NULL, max_points = 100000L) {
  data <- if (inherits(x, "simpower_curve") && !is.null(x$curve_results)) {
    retained <- lapply(x$curve_results, function(z) {
      d <- z$sim_data
      if (!is.null(d) && !"n_total" %in% names(d))
        d$n_total <- if (!is.null(z$n_total)) z$n_total else z$n
      d
    })
    retained <- retained[!vapply(retained, is.null, logical(1L))]
    if (length(retained)) do.call(rbind, retained) else NULL
  } else {
    x$sim_data
  }
  if (is.null(data))
    stop("Raw simulated data are not available. Re-run with keep_sim_data = TRUE.")
  data <- as.data.frame(data)
  if (!all(c("trial", "arm", "endpoint", "value") %in% names(data)))
    stop("'sim_data' does not contain the required trial, arm, endpoint, and value columns.")

  param <- if (inherits(x, "simpower") && !is.null(x$result)) {
    x$result$param
  } else if (!is.null(x$param)) {
    x$param
  } else {
    NULL
  }
  if (!is.null(param$list_comparator) && !any(display == "all")) {
    display <- .diagnostic_display(display)
    selected <- param$list_comparator[
      vapply(param$list_comparator, function(z) {
        gsub(" vs ", "_vs_", paste(z, collapse = "_vs_"), fixed = TRUE) %in% display
      }, logical(1))
    ]
    if (!length(selected))
      stop("Unknown comparator(s) in 'display'.")
    keep_arms <- unique(unlist(selected, use.names = FALSE))
    data <- data[data$arm %in% keep_arms, , drop = FALSE]
  }
  if (!is.null(endpoint)) data <- data[data$endpoint %in% endpoint, , drop = FALSE]
  if (!is.null(arm)) data <- data[data$arm %in% arm, , drop = FALSE]
  if (!nrow(data)) stop("No retained simulated observations match the requested filters.")
  if (nrow(data) > max_points) {
    set.seed(1L)
    data <- data[sample.int(nrow(data), max_points), , drop = FALSE]
  }
  data
}

.rate_distribution_data <- function(x, display = "all", endpoint = NULL,
                                    arm = NULL, max_points = 100000L) {
  data <- .distribution_sim_data(x, display, endpoint, arm, Inf)
  distribution <- unique(as.character(data$distribution))
  if (length(distribution) != 1L ||
      !distribution %in% c("pois", "nbinom"))
    stop("quantity = \"rate\" is available only for Poisson or negative binomial outcomes.")
  if (!all(c("exposure", "true_rate") %in% names(data)))
    stop("Retained count data do not contain exposure and true_rate columns.")
  split_key <- interaction(data$trial, data$arm, data$endpoint,
                           drop = TRUE, lex.order = TRUE)
  groups <- split(data, split_key)
  result <- lapply(groups, function(z) {
    data.frame(trial = z$trial[[1L]], arm = z$arm[[1L]],
               endpoint = z$endpoint[[1L]],
               value = sum(z$value) / sum(z$exposure),
               truth = z$true_rate[[1L]], stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, result)
  if (!is.null(arm) && !any(arm == "all")) {
    result <- result[order(match(as.character(result$arm), arm)), , drop = FALSE]
    result$arm <- factor(result$arm, levels = arm)
    rownames(result) <- NULL
  }
  if (nrow(result) > max_points) {
    set.seed(1L)
    result <- result[sample.int(nrow(result), max_points), , drop = FALSE]
  }
  result
}

.parameter_distribution_data <- function(x, endpoint = NULL, arm = NULL,
                                         parameter = "mean",
                                         max_points = 100000L) {
  if (!is.null(arm) && (!is.character(arm) || anyNA(arm)))
    stop("'arm' must be NULL, 'all', or a character vector of arm names.")
  arm_filter <- if (is.null(arm) || any(arm == "all")) NULL else arm
  data <- .distribution_sim_data(x, display = "all", max_points = Inf)
  available_endpoints <- unique(as.character(data$endpoint))
  if (!is.null(endpoint)) {
    unknown <- setdiff(endpoint, available_endpoints)
    if (length(unknown))
      stop("Unknown endpoint(s) in 'endpoints': ", paste(unknown, collapse = ", "))
    data <- data[data$endpoint %in% endpoint, , drop = FALSE]
  }
  if (!is.null(arm_filter)) {
    unknown <- setdiff(arm_filter, unique(as.character(data$arm)))
    if (length(unknown))
      stop("Unknown arm(s) in 'arms': ", paste(unknown, collapse = ", "))
    data <- data[data$arm %in% arm_filter, , drop = FALSE]
  }
  parameter <- match.arg(parameter, c("mean", "mu", "sigma", "dispersion", "rate"))
  if (parameter == "mu") parameter <- "mean"
  groups <- split(data, interaction(data$trial, data$arm, data$endpoint,
                                    drop = TRUE, lex.order = TRUE))
  context <- .distribution_parameters(x)
  param <- context$param
  param.d <- context$param.d
  rows <- lapply(groups, function(z) {
    values <- z$value
    mean_value <- mean(values, na.rm = TRUE)
    variance_value <- stats::var(values, na.rm = TRUE)
    distribution <- tolower(as.character(z$distribution[[1L]]))
    exposure <- if ("exposure" %in% names(z)) mean(z$exposure) else 1
    true_rate <- if ("true_rate" %in% names(z)) z$true_rate[[1L]] else NA_real_
    true_dispersion <- if ("true_dispersion" %in% names(z))
      z$true_dispersion[[1L]] else NA_real_
    truth <- NA_real_
    estimate <- switch(parameter,
      mean = mean_value,
      sigma = sqrt(variance_value),
    dispersion = if (distribution %in% c("pois", "nbinom"))
        (variance_value - mean_value) / mean_value^2 else
        variance_value / mean_value^2,
      rate = if (!is.na(true_rate)) mean_value / exposure else NA_real_
    )
    current_arm <- as.character(z$arm[[1L]])
    if (!is.null(param) && !is.null(param$mu) && !is.null(param$varcov)) {
      mu <- .endpoint_parameter(param$mu[[current_arm]], z$endpoint[[1L]])
      Sigma <- as.matrix(param$varcov[[current_arm]])
      endpoints <- param$ynames_list[[current_arm]]
      j <- match(z$endpoint[[1L]], endpoints)
      if (!is.na(j)) {
        sigma <- sqrt(Sigma[j, j])
        truth <- switch(parameter,
          mean = mu,
          sigma = sigma,
          dispersion = sigma^2 / mu^2,
          rate = NA_real_
        )
      }
    }
    if (is.na(truth) && !is.na(true_rate)) {
      truth <- switch(parameter,
        mean = true_rate * exposure,
        sigma = sqrt(true_rate * exposure + true_dispersion *
                       (true_rate * exposure)^2),
        dispersion = true_dispersion,
        rate = true_rate
      )
    }
    data.frame(trial = z$trial[[1L]], arm = z$arm[[1L]],
               endpoint = z$endpoint[[1L]], value = estimate,
               truth = truth, parameter = parameter,
               stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result <- result[is.finite(result$value), , drop = FALSE]
  if (nrow(result) > max_points) {
    set.seed(1L)
    result <- result[sample.int(nrow(result), max_points), , drop = FALSE]
  }
  if (!nrow(result)) stop("No parameter estimates could be calculated.")
  result
}

.target_endpoint_correlation <- function(x, arm, endpoint1, endpoint2) {
  context <- .distribution_parameters(x)
  param <- context$param
  corr <- if (!is.null(context$endpoint_corr)) context$endpoint_corr else
    if (!is.null(context$parameters)) context$parameters$endpoint_corr else NULL
  if (!is.null(corr)) {
    corr <- as.matrix(corr)
    if (!is.null(rownames(corr))) {
      i <- match(endpoint1, rownames(corr)); j <- match(endpoint2, colnames(corr))
    } else {
      i <- match(endpoint1, names(param$mu[[arm]])); j <- match(endpoint2, names(param$mu[[arm]]))
    }
    if (!is.na(i) && !is.na(j)) return(corr[i, j])
  }
  if (!is.null(param$varcov[[arm]])) {
    endpoints <- param$ynames_list[[arm]]
    i <- match(endpoint1, endpoints); j <- match(endpoint2, endpoints)
    Sigma <- as.matrix(param$varcov[[arm]])
    if (!is.na(i) && !is.na(j))
      return(Sigma[i, j] / sqrt(Sigma[i, i] * Sigma[j, j]))
  }
  NA_real_
}

.correlation_distribution_data <- function(x, display = "all", endpoint = NULL,
                                           arm = NULL, max_points = 100000L) {
  data <- .distribution_sim_data(x, display, endpoint, arm, Inf)
  if (!"subject" %in% names(data))
    stop("Retained simulated data do not contain subject identifiers; re-run with keep_sim_data = TRUE.")
  endpoints <- unique(as.character(data$endpoint))
  if (length(endpoints) < 2L)
    stop("At least two endpoints are required to calculate endpoint correlations.")
  grouping <- if ("n_total" %in% names(data))
    interaction(data$trial, data$arm, data$n_total, drop = TRUE,
                lex.order = TRUE) else
    interaction(data$trial, data$arm, drop = TRUE, lex.order = TRUE)
  groups <- split(data, grouping)
  rows <- list(); index <- 0L
  for (z in groups) {
    wide <- stats::reshape(z[, c("trial", "subject", "endpoint", "value")],
                    idvar = c("trial", "subject"), timevar = "endpoint",
                    direction = "wide")
    for (i in seq_len(length(endpoints) - 1L)) for (j in (i + 1L):length(endpoints)) {
      n1 <- paste0("value.", endpoints[[i]]); n2 <- paste0("value.", endpoints[[j]])
      if (!all(c(n1, n2) %in% names(wide))) next
      keep <- stats::complete.cases(wide[, c(n1, n2)])
      if (sum(keep) < 3L) next
      index <- index + 1L
      rows[[index]] <- data.frame(
        trial = z$trial[[1L]],
        n_total = if ("n_total" %in% names(z)) z$n_total[[1L]] else NA_real_,
        arm = z$arm[[1L]],
        pair = paste(endpoints[[i]], endpoints[[j]], sep = " vs "),
        value = stats::cor(wide[[n1]][keep], wide[[n2]][keep]),
        truth = .target_endpoint_correlation(x, z$arm[[1L]], endpoints[[i]], endpoints[[j]]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) stop("No endpoint correlations could be calculated.")
  result <- do.call(rbind, rows)
  if (!is.null(arm) && !any(arm == "all")) {
    result <- result[order(match(as.character(result$arm), arm)), , drop = FALSE]
    result$arm <- factor(result$arm, levels = arm)
    rownames(result) <- NULL
  }
  if (nrow(result) > max_points) {
    set.seed(1L)
    result <- result[sample.int(nrow(result), max_points), , drop = FALSE]
  }
  result
}

.distribution_parameters <- function(x) {
  if (inherits(x, "simpower") && !is.null(x$result))
    return(x$result)
  if (!is.null(x$param)) return(x)
  x
}

.distribution_comparators <- function(x, data) {
  context <- .distribution_parameters(x)
  param <- context$param
  comps <- if (!is.null(param) && is.list(param$list_comparator))
    param$list_comparator else if (is.list(x$comparisons)) x$comparisons else NULL
  if (!is.null(comps)) {
    if (is.null(names(comps))) names(comps) <- paste0("comparison_", seq_along(comps))
    names(comps) <- gsub(" vs ", "_vs_", names(comps), fixed = TRUE)
    return(comps)
  }
  arms <- unique(as.character(data$arm))
  if (length(arms) < 2L)
    stop("At least two arms are required to construct a comparator.")
  setNames(list(arms[seq_len(2L)]), paste(arms[seq_len(2L)], collapse = "_vs_"))
}

.comparator_display_from_arms <- function(x, data, arms = NULL) {
  comps <- .distribution_comparators(x, data)
  if (is.null(arms) || any(arms == "all")) return("all")
  if (!is.character(arms) || !length(arms) || anyNA(arms))
    stop("'arms' must be NULL, 'all', or a non-empty character vector of arm names.")
  available_arms <- unique(unlist(comps, use.names = FALSE))
  unknown <- setdiff(arms, available_arms)
  if (length(unknown))
    stop("Unknown arm(s) in 'arms': ", paste(unknown, collapse = ", "))
  selected <- names(comps)[vapply(comps, function(comp) {
    any(arms %in% comp)
  }, logical(1L))]
  if (!length(selected))
    stop("'arms' did not select any comparator panels.")
  gsub(" vs ", "_vs_", selected, fixed = TRUE)
}

.rate_ratio_data <- function(x, endpoint = NULL, arms = NULL, display = "all",
                             max_points = 100000L) {
  raw <- .distribution_sim_data(x, display = "all", endpoint = endpoint,
                                arm = NULL, max_points = Inf)
  distribution <- unique(tolower(as.character(raw$distribution)))
  if (length(distribution) != 1L || !distribution %in% c("pois", "nbinom"))
    stop("estimand = \"RR\" is available only for Poisson or negative-binomial outcomes.")
  comps <- .distribution_comparators(x, raw)
  selected <- .comparator_display_from_arms(x, raw, arms)
  if (any(selected == "all")) selected <- names(comps)
  selected <- if (any(display == "all")) selected else {
    display <- .diagnostic_display(display)
    unknown <- setdiff(display, gsub(" vs ", "_vs_", names(comps), fixed = TRUE))
    if (length(unknown))
      stop("Unknown comparator(s) in 'display': ", paste(unknown, collapse = ", "))
    intersect(selected, display)
  }
  if (!length(selected)) stop("No comparator panels remain after filtering.")
  rate_data <- .rate_distribution_data(x, display = "all", endpoint = endpoint,
                                       arm = unique(unlist(comps[selected], use.names = FALSE)),
                                       max_points = Inf)
  rows <- list(); index <- 0L
  for (comp_name in selected) {
    comp_index <- match(comp_name, gsub(" vs ", "_vs_", names(comps), fixed = TRUE))
    comp <- comps[[comp_index]]
    endpoints <- unique(rate_data$endpoint[
      rate_data$arm %in% comp
    ])
    for (ep in endpoints) {
      test <- rate_data[rate_data$arm == comp[[1L]] & rate_data$endpoint == ep,
                        c("trial", "value", "truth"), drop = FALSE]
      reference <- rate_data[rate_data$arm == comp[[2L]] & rate_data$endpoint == ep,
                             c("trial", "value", "truth"), drop = FALSE]
      names(test)[-1L] <- c("test", "test_truth")
      names(reference)[-1L] <- c("reference", "reference_truth")
      merged <- merge(test, reference, by = "trial")
      if (!nrow(merged)) next
      index <- index + 1L
      rows[[index]] <- data.frame(
        trial = merged$trial,
        value = merged$test / merged$reference,
        truth = merged$test_truth / merged$reference_truth,
        Comparator = comp_name, endpoint = ep, estimand = "RR",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) stop("No trial-level rate-ratio estimates were found in 'x'.")
  result <- do.call(rbind, rows)
  result <- result[is.finite(result$value), , drop = FALSE]
  if (nrow(result) > max_points) {
    set.seed(1L)
    result <- result[sample.int(nrow(result), max_points), , drop = FALSE]
  }
  if (!nrow(result)) stop("No finite rate-ratio estimates were available to plot.")
  result
}

.select_distribution_trials <- function(data, n_trials, seed) {
  if (is.null(n_trials)) return(data)
  if (!is.numeric(n_trials) || length(n_trials) != 1L ||
      !is.finite(n_trials) || n_trials < 1 || n_trials != as.integer(n_trials))
    stop("'n_trials' must be NULL or a positive integer.")
  trials <- unique(data$trial)
  if (length(trials) <= n_trials) return(data)
  if (!is.null(seed)) set.seed(as.integer(seed))
  selected <- sample(trials, n_trials)
  data[data$trial %in% selected, , drop = FALSE]
}

.endpoint_parameter <- function(value, endpoint) {
  if (is.matrix(value) || is.data.frame(value)) {
    if (!is.null(colnames(value)) && endpoint %in% colnames(value))
      return(as.numeric(value[, endpoint]))
    if (is.null(colnames(value)) && length(value) == 1L)
      return(as.numeric(value[[1L]]))
    j <- match(endpoint, colnames(value))
    if (is.na(j)) return(NA_real_)
    return(as.numeric(value[[j]]))
  }
  if (!is.null(names(value)) && endpoint %in% names(value))
    return(as.numeric(value[[endpoint]]))
  if (length(value) == 1L) return(as.numeric(value[[1L]]))
  if (is.null(names(value)) && length(endpoint) == 1L) {
    index <- if (grepl("^endpoint_[0-9]+$", endpoint)) {
      as.integer(sub("^endpoint_", "", endpoint))
    } else if (grepl("^y[0-9]+$", endpoint)) {
      as.integer(sub("^y", "", endpoint))
    } else {
      NA_integer_
    }
    if (!is.na(index) && index >= 1L && index <= length(value))
      return(as.numeric(value[[index]]))
  }
  j <- match(endpoint, names(value))
  if (is.na(j)) return(NA_real_)
  as.numeric(value[[j]])
}

.distribution_reference <- function(data, x, scale) {
  context <- .distribution_parameters(x)
  param <- context$param
  param.d <- context$param.d
  distribution <- if (!is.null(param.d)) param.d$distribution else
    unique(as.character(data$distribution))[[1L]]
  if (is.null(distribution) || is.na(distribution)) return(NULL)

  reference <- list()
  index <- 0L
  for (arm in unique(data$arm)) {
    for (endpoint in unique(data$endpoint)) {
      observed <- data$value[data$arm == arm & data$endpoint == endpoint]
      observed <- observed[is.finite(observed)]
      if (!length(observed)) next
      if (scale == "log" && distribution == "norm") next
      observed_plot <- if (scale == "log") log(observed) else observed
      if (distribution %in% c("norm", "lnorm")) {
        if (is.null(param) || is.null(param$mu) || is.null(param$varcov)) next
        mu <- as.numeric(param$mu[[as.character(arm)]])
        endpoints <- param$ynames_list[[as.character(arm)]]
        if (is.null(endpoints)) endpoints <- paste0("y", seq_along(mu))
        j <- match(as.character(endpoint), endpoints)
        if (is.na(j)) next
        Sigma <- as.matrix(param$varcov[[as.character(arm)]])
        sd_raw <- sqrt(Sigma[j, j])
      }
      if (distribution == "lnorm") {
        sd_log <- sqrt(log(Sigma[j, j] / mu[j]^2 + 1))
        mean_log <- log(mu[j]) - 0.5 * sd_log^2
        if (scale == "log") {
          grid <- seq(min(observed_plot), max(observed_plot), length.out = 200)
          density <- stats::dnorm(grid, mean_log, sd_log)
          cdf <- stats::pnorm(grid, mean_log, sd_log)
        } else {
          grid <- seq(min(observed_plot), max(observed_plot), length.out = 200)
          density <- stats::dlnorm(grid, mean_log, sd_log)
          cdf <- stats::plnorm(grid, mean_log, sdlog = sd_log)
        }
      } else if (distribution == "norm") {
        grid <- seq(min(observed_plot), max(observed_plot), length.out = 200)
        density <- stats::dnorm(grid, mu[j], sd_raw)
        cdf <- stats::pnorm(grid, mu[j], sd_raw)
      } else if (distribution %in% c("pois", "nbinom") &&
                 all(c("true_rate", "exposure", "true_dispersion") %in% names(data))) {
        rate <- unique(data$true_rate[data$arm == arm & data$endpoint == endpoint])[[1L]]
        exposure <- unique(data$exposure[data$arm == arm & data$endpoint == endpoint])[[1L]]
        dispersion <- unique(data$true_dispersion[data$arm == arm & data$endpoint == endpoint])[[1L]]
        upper <- max(stats::qpois(0.999, rate * exposure), max(observed), 1)
        grid <- 0:min(upper, 200)
        if (distribution == "pois") {
          density <- stats::dpois(grid, rate * exposure)
          cdf <- stats::ppois(grid, rate * exposure)
        } else {
          density <- stats::dnbinom(grid, size = 1 / dispersion,
                                    mu = rate * exposure)
          cdf <- stats::pnbinom(grid, size = 1 / dispersion,
                                mu = rate * exposure)
        }
      }
      index <- index + 1L
      reference[[index]] <- data.frame(
        arm = as.character(arm), endpoint = as.character(endpoint),
        x = grid, density = density, cdf = cdf,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(reference)) do.call(rbind, reference) else NULL
}

.estimand_data <- function(x, display = "all", endpoint = NULL,
                           max_points = 100000L) {
  context <- .distribution_parameters(x)
  param <- context$param
  param.d <- context$param.d
  if (is.null(param) || is.null(param.d))
    stop("Estimand distributions require continuous simulation parameters.")
  tables <- if (inherits(x, "simpower_curve")) {
    lapply(x$curve_results, function(z) z$result$table.test)
  } else if (inherits(x, "simpower")) {
    list(x$result$table.test)
  } else {
    list(x$table.test)
  }
  comps <- param$list_comparator
  display <- .diagnostic_display(display)
  rows <- list()
  index <- 0L
  for (tab in tables) {
    for (i in seq_along(comps)) {
      comp <- comps[[i]]
      comp_name <- gsub(" vs ", "_vs_", paste(comp, collapse = "_vs_"), fixed = TRUE)
      if (!any(display == "all") && !comp_name %in% display) next
      endpoints <- if (!is.null(param$list_y_comparator))
        param$list_y_comparator[[i]] else names(param$mu[[comp[[1L]]]])
      for (ep in endpoints) {
        if (!is.null(endpoint) && !ep %in% endpoint) next
        prefix <- paste0("mu_", ep, "_")
        col1 <- paste0(prefix, comp[[1L]], "Comp:", paste(comp, collapse = " vs "))
        col2 <- paste0(prefix, comp[[2L]], "Comp:", paste(comp, collapse = " vs "))
        if (!all(c(col1, col2) %in% names(tab))) next
        estimate1 <- tab[[col1]]
        estimate2 <- tab[[col2]]
        if (param.d$ctype == "DOM") {
          estimate <- estimate1 - estimate2
          truth <- .endpoint_parameter(param$mu[[comp[[1L]]]], ep) -
            .endpoint_parameter(param$mu[[comp[[2L]]]], ep)
          estimand <- "DOM"
        } else {
          estimate <- if (param.d$distribution == "lnorm")
            exp(estimate1 - estimate2) else estimate1 / estimate2
          truth <- .endpoint_parameter(param$mu[[comp[[1L]]]], ep) /
            .endpoint_parameter(param$mu[[comp[[2L]]]], ep)
          estimand <- "ROM"
        }
        index <- index + 1L
        rows[[index]] <- data.frame(
          trial = seq_along(estimate), value = as.numeric(estimate), truth = truth,
          Comparator = comp_name, endpoint = ep, estimand = estimand,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) stop("No trial-level estimand estimates were found in 'x'.")
  data <- do.call(rbind, rows)
  data <- data[is.finite(data$value), , drop = FALSE]
  if (nrow(data) > max_points) {
    set.seed(1L)
    data <- data[sample.int(nrow(data), max_points), , drop = FALSE]
  }
  data
}

.t_value_distribution_data <- function(x, display = "all", endpoint = NULL,
                                       max_points = 100000L) {
  context <- .distribution_parameters(x)
  param <- context$param
  param.d <- context$param.d
  if (is.null(param) || is.null(param.d))
    stop("Test-statistic distributions require continuous simulation parameters.")
  if (!identical(param.d$dtype, "parallel"))
    stop("quantity = \"t_value\" is currently available for parallel designs only.")
  if (!param.d$distribution %in% c("norm", "lnorm"))
    stop("quantity = \"t_value\" is available for continuous outcomes only.")
  if (!param.d$ctype %in% c("DOM", "ROM"))
    stop("quantity = \"t_value\" requires ctype = \"DOM\" or \"ROM\".")

  data <- .distribution_sim_data(x, display = display, endpoint = endpoint,
                                 arm = NULL, max_points = Inf)
  comps <- param$list_comparator
  display <- .diagnostic_display(display)
  rows <- list(); index <- 0L
  value_at <- function(z, ep, j) {
    if (!is.null(names(z)) && ep %in% names(z)) return(as.numeric(z[[ep]]))
    as.numeric(z[[j]])
  }
  alpha_at <- function(i, ep, j) {
    alpha <- as.numeric(param.d$alpha[[1L]])
    m <- length(param$list_y_comparator[[i]])
    adjust <- as.character(param.d$adjust)
    if (adjust == "bon") return(alpha / m)
    if (adjust == "sid") return(1 - (1 - alpha)^(1 / m))
    if (adjust == "k") return(as.numeric(param.d$k[[i]]) * alpha / m)
    if (adjust == "t") {
      k <- as.numeric(param.d$k[[i]])
      return(alpha / (m - k + 1))
    }
    if (adjust == "seq" && !is.null(param$weight_seq)) {
      weights <- param$weight_seq
      if (!is.null(names(weights)) && ep %in% names(weights))
        return(alpha * as.numeric(weights[[ep]]))
      return(alpha * as.numeric(weights[[j]]))
    }
    alpha
  }
  for (i in seq_along(comps)) {
    comp <- comps[[i]]
    comp_name <- gsub(" vs ", "_vs_", paste(comp, collapse = "_vs_"), fixed = TRUE)
    if (!any(display == "all") && !comp_name %in% display) next
    endpoints <- param$list_y_comparator[[i]]
    lower <- param.d$list_lequi.tol[[i]]
    upper <- param.d$list_uequi.tol[[i]]
    for (j in seq_along(endpoints)) {
      ep <- endpoints[[j]]
      lower_j <- value_at(lower, ep, j)
      upper_j <- value_at(upper, ep, j)
      alpha_j <- alpha_at(i, ep, j)
      trials <- sort(unique(data$trial))
      for (trial in trials) {
        test_values <- data$value[data$trial == trial &
                                  data$arm == comp[[1L]] & data$endpoint == ep]
        reference_values <- data$value[data$trial == trial &
                                       data$arm == comp[[2L]] & data$endpoint == ep]
        if (length(test_values) < 2L || length(reference_values) < 2L) next
        mean_test <- mean(test_values)
        mean_reference <- mean(reference_values)
        sd_test <- stats::sd(test_values)
        sd_reference <- stats::sd(reference_values)
        n_test <- length(test_values)
        n_reference <- length(reference_values)
        lower_work <- lower_j
        upper_work <- upper_j
        if (param.d$distribution == "lnorm") {
          test_values <- log(test_values)
          reference_values <- log(reference_values)
          mean_test <- mean(test_values)
          mean_reference <- mean(reference_values)
          sd_test <- stats::sd(test_values)
          sd_reference <- stats::sd(reference_values)
          lower_work <- log(lower_work)
          upper_work <- log(upper_work)
        }
        if (param.d$ctype == "DOM") {
          pooled <- ((n_test - 1) * sd_test^2 +
                     (n_reference - 1) * sd_reference^2) /
            (n_test + n_reference - 2)
          se <- if (isTRUE(param.d$vareq)) {
            sqrt(pooled * (1 / n_test + 1 / n_reference))
          } else {
            sqrt(sd_test^2 / n_test + sd_reference^2 / n_reference)
          }
          t_lower <- (mean_test - mean_reference - lower_work) / se
          t_upper <- (mean_test - mean_reference - upper_work) / se
          df <- if (isTRUE(param.d$vareq)) n_test + n_reference - 2 else
            (sd_test^2 / n_test + sd_reference^2 / n_reference)^2 /
              ((sd_test^2 / n_test)^2 / (n_test - 1) +
               (sd_reference^2 / n_reference)^2 / (n_reference - 1))
        } else {
          pooled <- ((n_test - 1) * sd_test^2 +
                     (n_reference - 1) * sd_reference^2) /
            (n_test + n_reference - 2)
          se_lower <- if (isTRUE(param.d$vareq)) {
            sqrt(pooled * (1 / n_test + lower_work^2 / n_reference))
          } else {
            sqrt(sd_test^2 / n_test + lower_work^2 * sd_reference^2 / n_reference)
          }
          se_upper <- if (isTRUE(param.d$vareq)) {
            sqrt(pooled * (1 / n_test + upper_work^2 / n_reference))
          } else {
            sqrt(sd_test^2 / n_test + upper_work^2 * sd_reference^2 / n_reference)
          }
          t_lower <- (mean_test - mean_reference * lower_work) / se_lower
          t_upper <- (mean_test - mean_reference * upper_work) / se_upper
          df <- n_test + n_reference - 2
        }
        values <- c(lower = t_lower, upper = t_upper)
        critical <- c(lower = stats::qt(1 - alpha_j, df),
                      upper = stats::qt(alpha_j, df))
        for (boundary in names(values)) {
          index <- index + 1L
          rows[[index]] <- data.frame(
            trial = trial, value = as.numeric(values[[boundary]]),
            critical = as.numeric(critical[[boundary]]),
            boundary = boundary, Comparator = comp_name, endpoint = ep,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(rows))
    stop("No test statistics could be calculated from the retained simulated data.")
  result <- do.call(rbind, rows)
  result <- result[is.finite(result$value), , drop = FALSE]
  if (nrow(result) > max_points) {
    set.seed(1L)
    result <- result[sample.int(nrow(result), max_points), , drop = FALSE]
  }
  if (!nrow(result)) stop("No finite test statistics were available to plot.")
  result
}

#' Plot distributions of retained simulated observations
#'
#' @description Creates density, histogram, ECDF, or Q--Q plots for a selected
#' trial-level quantity. The `estimand` argument selects an arm parameter, a
#' reconstructed TOST statistic, a continuous comparison, a count rate ratio,
#' or endpoint correlations.
#' @param x A result returned by `sampleSize()` or `simPower()` with
#'   `keep_sim_data = TRUE`.
#' @param estimand Quantity to display: `"mu"` (or `"mean"`), `"sigma"`,
#'   `"dispersion"`, `"rate"`, `"t_value"`,
#'   `"DOM"`, `"ROM"`, `"RR"`, or `"correlation"`. `"mean"` and `"sigma"`
#'   display arm-specific simulation parameters; `"DOM"` and `"ROM"` display
#'   continuous comparison estimates; `"RR"` displays count rate ratios; and
#'   `"correlation"` displays endpoint correlations.
#' @param arms Optional character vector of arm names. For `mean` and `sigma`,
#'   these are the arms to facet. For comparison estimands, these select
#'   comparator panels containing at least one supplied arm. If `NULL`, all
#'   are shown.
#' @param endpoints Optional character vector of endpoint names. If `NULL`,
#'   all endpoints are displayed.
#' @param ... Optional plotting controls, including `type`, `show_reference`,
#'   `n_trials`, `seed`, and `max_points`. Legacy selector names are accepted
#'   for backward compatibility.
#' @details Every comparator follows the package
#'   convention `c(test, reference)`. For `t_value`, the lower
#'   and upper TOST statistics are reconstructed from the retained outcomes
#'   using the same continuous parallel formulas as the test. The dashed lines
#'   are the corresponding one-sided critical values. The plotted DOM is
#'   `test - reference`
#'   and the plotted ROM is `test / reference`; reversing the comparator reverses
#'   the displayed estimand. The plotted RR is `rate_test / rate_reference`.
#' @return A `ggplot` object.
#' @export
plot_distribution <- function(x, estimand = "mu", arms = NULL,
                              endpoints = NULL, ...) {
  dots <- list(...)
  type <- if (is.null(dots$type)) "density" else dots$type
  show_reference <- if (is.null(dots$show_reference)) TRUE else dots$show_reference
  max_points <- if (is.null(dots$max_points)) 100000L else dots$max_points
  n_trials <- dots$n_trials
  seed <- if (is.null(dots$seed)) 1 else dots$seed
  if (any(c("quantity", "parameter", "endpoint", "arm", "display", "scale") %in%
          names(dots)))
    stop("Use plot_distribution(x, estimand = ..., arms = ..., endpoints = ...).")
  if ("overall" %in% names(dots))
    stop("'overall' is not available for plot_distribution(); use plot_stability() or plot_mc_error().")
  type <- match.arg(type, c("density", "histogram", "ecdf", "qq"))
  estimand <- match.arg(toupper(as.character(estimand)),
                        c("MEAN", "MU", "SIGMA", "DISPERSION", "T_VALUE", "DOM", "ROM", "RR",
                          "CORRELATION", "OUTCOME", "RATE"))
  estimand <- tolower(estimand)
  if (estimand == "mu") estimand <- "mean"
  if (estimand == "correlation") {
    if (type == "qq") stop("'type = \"qq\"' is not available for correlation plots.")
    return(.plot_correlation(
      x, type = type, endpoint = endpoints, arm = arms,
      show_reference = show_reference, max_points = as.integer(max_points)
    ))
  }
  quantity <- switch(estimand,
    mean = "parameter", sigma = "parameter", dispersion = "parameter",
    rate = "parameter",
    t_value = "t_value",
    dom = "estimand", rom = "estimand", rr = "rate_ratio",
    outcome = "outcome"
  )
  parameter <- if (estimand %in% c("mean", "sigma", "dispersion", "rate"))
    estimand else "mean"
  endpoint <- endpoints
  arm <- arms
  display <- "all"
  scale <- "original"
  if (!is.logical(show_reference) || length(show_reference) != 1L ||
      is.na(show_reference))
    stop("'show_reference' must be TRUE or FALSE.")
  if (!is.numeric(max_points) || length(max_points) != 1L ||
      max_points < 1 || max_points != as.integer(max_points))
    stop("'max_points' must be a positive integer.")
  if (quantity == "outcome") {
    data <- .distribution_sim_data(
      x, display = display, endpoint = endpoint, arm = arm,
      max_points = as.integer(max_points)
    )
    data <- .select_distribution_trials(data, n_trials, seed)
    reference <- if (show_reference)
      .distribution_reference(data, x, scale) else NULL
    if (scale == "log") {
      if (any(data$value <= 0, na.rm = TRUE))
        stop("'scale = \"log\"' requires strictly positive simulated values.")
      data$value <- log(data$value)
    }
    arm_levels <- if (!is.null(arms) && !any(arms == "all")) {
      arms[arms %in% unique(as.character(data$arm))]
    } else {
      if (is.factor(data$arm)) levels(data$arm) else unique(as.character(data$arm))
    }
    data$arm <- factor(as.character(data$arm), levels = arm_levels)
    data$endpoint <- factor(data$endpoint, levels = unique(data$endpoint))
    facet_formula <- stats::as.formula("endpoint ~ arm")
    x_label <- if (scale == "log") "log(simulated outcome)" else "Simulated outcome"
  } else if (quantity == "rate") {
    if (scale == "log")
      stop("scale = \"log\" is not supported for quantity = \"rate\".")
    data <- .rate_distribution_data(x, display, endpoint, arm, as.integer(max_points))
    data <- .select_distribution_trials(data, n_trials, seed)
    reference <- if (show_reference) unique(data[, c("arm", "endpoint", "truth")]) else NULL
    arm_levels <- if (!is.null(arms) && !any(arms == "all")) {
      arms[arms %in% unique(as.character(data$arm))]
    } else if (is.factor(data$arm)) {
      levels(data$arm)
    } else {
      unique(as.character(data$arm))
    }
    data$arm <- factor(as.character(data$arm), levels = arm_levels)
    data$endpoint <- factor(data$endpoint, levels = unique(data$endpoint))
    facet_formula <- stats::as.formula("endpoint ~ arm")
    x_label <- "Estimated event rate"
  } else if (quantity == "t_value") {
    if (!is.null(arm)) {
      filter_data <- .distribution_sim_data(x, display = "all",
                                            endpoint = endpoint, max_points = Inf)
      display <- .comparator_display_from_arms(x, filter_data, arm)
    }
    data <- .t_value_distribution_data(x, display = display, endpoint = endpoint,
                                       max_points = as.integer(max_points))
    data <- .select_distribution_trials(data, n_trials, seed)
    reference <- if (show_reference)
      unique(data[, c("Comparator", "endpoint", "boundary", "critical")]) else NULL
    data$Comparator <- factor(data$Comparator, levels = unique(data$Comparator))
    data$endpoint <- factor(data$endpoint, levels = unique(data$endpoint))
    data$boundary <- factor(data$boundary, levels = c("lower", "upper"))
    facet_formula <- stats::as.formula("endpoint ~ Comparator")
    x_label <- "TOST t-statistic"
  } else if (quantity == "parameter") {
    data <- .parameter_distribution_data(x, endpoint = endpoint, arm = arm,
                                         parameter = parameter,
                                         max_points = as.integer(max_points))
    data <- .select_distribution_trials(data, n_trials, seed)
    reference <- if (show_reference) unique(data[, c("arm", "endpoint", "truth")]) else NULL
    arm_levels <- if (!is.null(arms) && !any(arms == "all")) {
      arms[arms %in% unique(as.character(data$arm))]
    } else if (is.factor(data$arm)) {
      levels(data$arm)
    } else {
      unique(as.character(data$arm))
    }
    data$arm <- factor(as.character(data$arm), levels = arm_levels)
    data$endpoint <- factor(data$endpoint, levels = unique(data$endpoint))
    facet_formula <- stats::as.formula("endpoint ~ arm")
    x_label <- paste("Simulated", parameter)
  } else if (quantity == "rate_ratio") {
    data <- .rate_ratio_data(x, endpoint = endpoint, arms = arm,
                             display = display,
                             max_points = as.integer(max_points))
    data <- .select_distribution_trials(data, n_trials, seed)
    reference <- if (show_reference)
      unique(data[, c("Comparator", "endpoint", "truth"), drop = FALSE]) else NULL
    data$Comparator <- factor(data$Comparator, levels = unique(data$Comparator))
    data$endpoint <- factor(data$endpoint, levels = unique(data$endpoint))
    facet_formula <- stats::as.formula("endpoint ~ Comparator")
    x_label <- "Estimated event-rate ratio"
  } else {
    if (!is.null(arm)) {
      filter_data <- .distribution_sim_data(x, display = "all",
                                            endpoint = endpoint, max_points = Inf)
      display <- .comparator_display_from_arms(x, filter_data, arm)
    }
    data <- .estimand_data(x, display = display, endpoint = endpoint,
                           max_points = as.integer(max_points))
    data <- .select_distribution_trials(data, n_trials, seed)
    data_estimand <- tolower(unique(data$estimand))
    if (length(data_estimand) != 1L || data_estimand != estimand) {
      context <- .distribution_parameters(x)
      object_ctype <- if (!is.null(context$param.d$ctype))
        toupper(as.character(context$param.d$ctype[[1L]])) else "unknown"
      actual <- toupper(data_estimand[[1L]])
      warning(
        "plot_distribution(): requested estimand '", toupper(estimand),
        "' is incompatible with this object's ctype = '", object_ctype,
        "'. The object contains '", actual, "' estimates, so '", actual,
        "' will be plotted instead. Use plot_distribution(x, estimand = \"",
        tolower(data_estimand[[1L]]), "\") or create the object with ctype = \"",
        actual, "\".",
        call. = FALSE
      )
      estimand <- data_estimand[[1L]]
    }
    reference <- if (show_reference)
      unique(data[, c("Comparator", "endpoint", "truth"), drop = FALSE]) else NULL
    data$Comparator <- factor(data$Comparator, levels = unique(data$Comparator))
    data$endpoint <- factor(data$endpoint, levels = unique(data$endpoint))
    facet_formula <- stats::as.formula("endpoint ~ Comparator")
    x_label <- if (estimand == "dom")
      "Estimated mean difference" else "Estimated mean ratio"
  }
  p <- if (quantity == "t_value") {
    ggplot2::ggplot(data, ggplot2::aes(x = value, color = boundary,
                                       fill = boundary, group = boundary))
  } else {
    ggplot2::ggplot(data, ggplot2::aes(x = value))
  }
  if (type == "density") {
    p <- p + if (quantity == "t_value") {
      ggplot2::geom_density(alpha = 0.25, na.rm = TRUE)
    } else {
      ggplot2::geom_density(fill = "#0072B2", alpha = 0.35,
                            color = "#0072B2", na.rm = TRUE)
    }
  } else if (type == "histogram") {
    p <- p + ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30,
      position = if (quantity == "t_value") "identity" else "stack",
      fill = if (quantity == "t_value") NULL else "#0072B2",
      color = "white",
      alpha = if (quantity == "t_value") 0.35 else 0.8,
      na.rm = TRUE)
  } else if (type == "ecdf") {
    p <- p + ggplot2::stat_ecdf(
      color = if (quantity == "t_value") NULL else "#0072B2",
      linewidth = 0.7, na.rm = TRUE) +
      ggplot2::scale_y_continuous(labels = scales::label_percent())
  } else {
    p <- p + ggplot2::stat_qq(
      data = data,
      ggplot2::aes(sample = value),
      inherit.aes = quantity == "t_value",
      color = if (quantity == "t_value") NULL else "#0072B2",
      alpha = 0.55, na.rm = TRUE
    ) +
      ggplot2::stat_qq_line(
        data = data,
        ggplot2::aes(sample = value),
        inherit.aes = quantity == "t_value",
        color = if (quantity == "t_value") NULL else "#D55E00",
        linewidth = 0.6, na.rm = TRUE
      )
  }
  if (show_reference && !is.null(reference) && type != "qq") {
    if (quantity == "outcome") {
      reference$arm <- factor(reference$arm, levels = levels(data$arm))
      reference$endpoint <- factor(reference$endpoint, levels = levels(data$endpoint))
      if (type == "ecdf") {
        p <- p + ggplot2::geom_line(
          data = reference,
          ggplot2::aes(x = x, y = cdf, group = interaction(arm, endpoint)),
          inherit.aes = FALSE, color = "#D55E00", linewidth = 0.7,
          linetype = "dashed"
        )
      } else {
        p <- p + ggplot2::geom_line(
          data = reference,
          ggplot2::aes(x = x, y = density, group = interaction(arm, endpoint)),
          inherit.aes = FALSE, color = "#D55E00", linewidth = 0.7,
          linetype = "dashed"
        )
      }
    } else if (quantity %in% c("rate", "parameter")) {
      reference$arm <- factor(reference$arm, levels = levels(data$arm))
      reference$endpoint <- factor(reference$endpoint, levels = levels(data$endpoint))
      p <- p + ggplot2::geom_vline(
        data = reference, ggplot2::aes(xintercept = truth), inherit.aes = FALSE,
        color = "#D55E00", linewidth = 0.7, linetype = "dashed"
      )
    } else if (quantity == "t_value") {
      reference$Comparator <- factor(reference$Comparator,
                                     levels = levels(data$Comparator))
      reference$endpoint <- factor(reference$endpoint,
                                   levels = levels(data$endpoint))
      reference$boundary <- factor(reference$boundary,
                                   levels = levels(data$boundary))
      p <- p + ggplot2::geom_vline(
        data = reference,
        ggplot2::aes(xintercept = critical, color = boundary),
        inherit.aes = FALSE, linewidth = 0.7, linetype = "dashed"
      )
    } else {
      reference$Comparator <- factor(reference$Comparator,
                                     levels = levels(data$Comparator))
      reference$endpoint <- factor(reference$endpoint,
                                   levels = levels(data$endpoint))
      if (type != "qq") {
        p <- p + ggplot2::geom_vline(
          data = reference,
          ggplot2::aes(xintercept = truth), inherit.aes = FALSE,
          color = "#D55E00", linewidth = 0.7, linetype = "dashed"
        )
      }
    }
  }
  if (quantity == "t_value") {
    p <- p + ggplot2::scale_color_manual(
      values = c(lower = "#0072B2", upper = "#D55E00")
    ) + ggplot2::scale_fill_manual(
      values = c(lower = "#0072B2", upper = "#D55E00")
    )
  }
  p + ggplot2::facet_grid(facet_formula, scales = "free") +
    ggplot2::labs(x = x_label, y = NULL,
                  color = if (quantity == "t_value") "TOST boundary" else NULL,
                  fill = if (quantity == "t_value") "TOST boundary" else NULL) +
    .diagnostic_theme()
}

#' Plot simulated endpoint correlations
#'
#' @description Shows the distribution of within-trial, within-arm endpoint
#' correlations. The dashed line is the user-specified target correlation.
#' @param x A result returned by `sampleSize()` or `simPower()` with
#'   `keep_sim_data = TRUE`.
#' @param type Plot type: `"density"`, `"histogram"`, or `"ecdf"`.
#' @param display Comparator names to display, or `"all"`.
#' @param endpoint Optional endpoint names to display.
#' @param arm Optional arm names to display.
#' @param show_reference Logical; overlay the target correlation.
#' @param max_points Maximum number of retained observations used.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object.
#' @export
.plot_correlation <- function(x, type = c("density", "histogram", "ecdf"),
                             display = "all", endpoint = NULL, arm = NULL,
                             show_reference = TRUE, max_points = 100000L,
                             ...) {
  type <- match.arg(type)
  data <- .correlation_distribution_data(x, display, endpoint, arm,
                                         as.integer(max_points))
  arm_levels <- if (!is.null(arm) && !any(arm == "all")) {
    arm[arm %in% unique(as.character(data$arm))]
  } else {
    if (is.factor(data$arm)) levels(data$arm) else unique(as.character(data$arm))
  }
  data$arm <- factor(as.character(data$arm), levels = arm_levels)
  data$pair <- factor(data$pair, levels = unique(data$pair))
  p <- ggplot2::ggplot(data, ggplot2::aes(x = value))
  if (type == "density") {
    p <- p + ggplot2::geom_density(fill = "#0072B2", alpha = 0.35,
                                   color = "#0072B2", na.rm = TRUE)
  } else if (type == "histogram") {
    p <- p + ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30,
      fill = "#0072B2", color = "white", alpha = 0.8, na.rm = TRUE)
  } else {
    p <- p + ggplot2::stat_ecdf(color = "#0072B2", linewidth = 0.7,
                                na.rm = TRUE) +
      ggplot2::scale_y_continuous(labels = scales::label_percent())
  }
  if (show_reference) {
    p <- p + ggplot2::geom_vline(
      data = unique(data[, c("arm", "pair", "truth")]),
      ggplot2::aes(xintercept = truth), inherit.aes = FALSE,
      color = "#D55E00", linewidth = 0.7, linetype = "dashed"
    )
  }
  facet_formula <- if (any(!is.na(data$n_total)))
    stats::as.formula("n_total + pair ~ arm") else
    stats::as.formula("pair ~ arm")
  p + ggplot2::facet_grid(facet_formula, scales = "free") +
    ggplot2::labs(x = "Within-trial endpoint correlation", y = NULL) +
    .diagnostic_theme()
}
