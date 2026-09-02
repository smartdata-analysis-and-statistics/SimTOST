# Helpers shared by the simulation-diagnostic plots.
.diagnostic_tables <- function(x) {
  if (is.list(x) && !inherits(x, c("simss", "simpower", "simpower_curve"))) {
    if (!length(x) || !all(vapply(x, function(z) {
      inherits(z, c("simss", "simpower", "simpower_curve"))
    }, logical(1)))) {
      stop("A diagnostic list must contain simss, simpower, or simpower_curve objects.")
    }
    info <- lapply(x, .diagnostic_tables)
    return(list(
      tables = unlist(lapply(info, `[[`, "tables"), recursive = FALSE),
      n_total = unlist(lapply(info, `[[`, "n_total"), use.names = FALSE),
      nsim = unlist(lapply(info, `[[`, "nsim"), use.names = FALSE)
    ))
  }
  if (inherits(x, "simpower_curve")) {
    tables <- lapply(x$curve_results, function(z) z$result$table.test)
    n_total <- vapply(x$curve_results, function(z) {
      if (!is.null(z$n_total)) return(as.numeric(z$n_total))
      if (!is.null(z$result$response$n_total))
        return(as.numeric(z$result$response$n_total[[1L]]))
      as.numeric(z$n)
      }, numeric(1))
    nsim <- vapply(x$curve_results, function(z) {
      if (!is.null(z$nsim)) return(as.numeric(z$nsim))
      if (!is.null(z$result$param.d$nsim)) return(as.numeric(z$result$param.d$nsim))
      NA_real_
    }, numeric(1))
    return(list(tables = tables, n_total = n_total, nsim = nsim))
  }
  if (inherits(x, "simpower") && !is.null(x$result$table.test))
    return(list(tables = list(x$result$table.test),
                n_total = if (is.null(x$n_total)) NA_real_ else x$n_total,
                nsim = if (is.null(x$nsim)) NA_real_ else x$nsim))
  if (inherits(x, "simss") && !is.null(x$table.test)) {
    suggested_n_total <- NA_real_
    if (is.data.frame(x$response) && nrow(x$response) > 0L &&
        "n_total" %in% names(x$response) &&
        is.finite(as.numeric(x$response$n_total[[1L]]))) {
      suggested_n_total <- as.numeric(x$response$n_total[[1L]])
    } else if (is.data.frame(x$table.iter) && nrow(x$table.iter) > 0L &&
               "n_total" %in% names(x$table.iter)) {
      suggested_n_total <- as.numeric(utils::tail(x$table.iter$n_total, 1L))
    }
    table <- x$table.test
    if (is.finite(suggested_n_total) && "n_total" %in% names(table)) {
      table <- table[as.numeric(table$n_total) == suggested_n_total, , drop = FALSE]
    }
    return(list(
      tables = list(table), n_total = suggested_n_total,
      nsim = if (!is.null(x$param.d$nsim)) x$param.d$nsim else NA_real_
    ))
  }
  stop("'x' must be a simss, simpower, or simpower_curve object.")
}

.diagnostic_display <- function(display) {
  if (!is.character(display) || length(display) == 0L || anyNA(display))
    stop("'display' must be 'all' or a non-empty character vector of comparator names.")
  unique(gsub(" vs ", "_vs_", display, fixed = TRUE))
}

.diagnostic_long <- function(x) {
  info <- .diagnostic_tables(x)
  out <- vector("list", length(info$tables))
  for (i in seq_along(info$tables)) {
    tab <- as.data.frame(info$tables[[i]])
    if (!nrow(tab)) next
    decision_cols <- grep("Comp:", names(tab), value = TRUE)
    decision_cols <- decision_cols[
      !grepl("^(mu_|sd_|eql_|equ_)", decision_cols)
    ]
    # The stored `totaly` column is the all-comparators decision. Include it
    # only when there is more than one comparator, so a one-comparison plot
    # does not duplicate the individual comparator panel.
    comparator_names <- sub(".*Comp:", "", decision_cols)
    if ("totaly" %in% names(tab) && length(unique(comparator_names)) > 1L)
      decision_cols <- c(decision_cols, "totaly")
    if (!length(decision_cols) && "totaly" %in% names(tab))
      decision_cols <- "totaly"
    if (!length(decision_cols)) next

    n_total <- if ("n_total" %in% names(tab)) tab$n_total else
      rep(info$n_total[[min(i, length(info$n_total))]], nrow(tab))
    nsim <- rep(info$nsim[[min(i, length(info$nsim))]], nrow(tab))
    trial <- seq_len(nrow(tab))
    out[[i]] <- do.call(rbind, lapply(decision_cols, function(col) {
      comparator <- sub(".*Comp:", "", col)
      endpoint <- sub("Comp:.*", "", col)
      if (col == "totaly") {
        comparator <- "All comparators"
        endpoint <- "Total"
      }
      if (comparator == "totaly") comparator <- "All comparators"
      if (endpoint == "totaly") endpoint <- "Total"
        data.frame(
        trial = trial, n_total = n_total, nsim = nsim,
        value = as.numeric(tab[[col]]),
        Comparator = gsub(" vs ", "_vs_", comparator, fixed = TRUE),
        Endpoint = endpoint, stringsAsFactors = FALSE
      )
    }))
  }
  out <- out[!vapply(out, is.null, logical(1L))]
  if (!length(out)) stop("No trial-level decision columns were found in 'x'.")
  do.call(rbind, out)
}

.filter_diagnostic_display <- function(data, display) {
  display <- .diagnostic_display(display)
  data$Comparator <- gsub(" vs ", "_vs_", data$Comparator, fixed = TRUE)
  if (any(display == "all")) return(data)
  available <- unique(data$Comparator)
  unknown <- setdiff(display, available)
  if (length(unknown))
    stop("Unknown comparator(s) in 'display': ", paste(unknown, collapse = ", "))
  data[data$Comparator %in% display, , drop = FALSE]
}

.filter_diagnostic_endpoint <- function(data, endpoint) {
  if (!is.character(endpoint) || length(endpoint) == 0L || anyNA(endpoint))
    stop("'endpoint' must be 'all' or a non-empty character vector of endpoint names.")
  endpoint <- unique(endpoint)
  if (any(endpoint == "all")) return(data)
  available <- unique(data$Endpoint)
  unknown <- setdiff(endpoint, available)
  if (length(unknown))
    stop("Unknown endpoint(s) in 'endpoint': ", paste(unknown, collapse = ", "))
  data[data$Endpoint %in% endpoint, , drop = FALSE]
}

.diagnostic_theme <- function() {
  ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#EAF1F8", color = NA),
      strip.text = ggplot2::element_text(face = "bold", color = "#20364D"),
      axis.title = ggplot2::element_text(face = "bold", color = "#20364D"),
      legend.position = "bottom"
    )
}

.diagnostic_sample_size_label <- function(data) {
  values <- sort(unique(data$n_total[is.finite(data$n_total)]))
  if (!length(values)) return("Suggested total sample size is unavailable")
  if (length(values) == 1L)
    return(paste0("Suggested total sample size = ", values[[1L]]))
  paste0("Compared total sample sizes = ", paste(values, collapse = ", "))
}

.diagnostic_endpoint_scale <- function(data) {
  ggplot2::scale_color_manual(
    values = c(stats::setNames(grDevices::hcl.colors(
      length(setdiff(levels(data$Endpoint), "Total")), "Dark 3"),
      setdiff(levels(data$Endpoint), "Total")), Total = "black"),
    drop = FALSE
  )
}

#' Plot simulation stability
#'
#' Displays cumulative achieved power over simulated trials. A stable curve
#' approaches a horizontal plateau, while large late movements indicate that
#' more simulations may be needed.
#' @param x A `simss`, `simpower`, or `simpower_curve` object, or a list of
#'   such objects. A list can compare simulations with different `n` and
#'   `nsim` values.
#' @param target_power Target power shown as a horizontal reference line.
#' @param display Comparator names to display, or `"all"`.
#' @param endpoint Endpoint names to display. The default, `"Total"`, shows
#'   only the aggregate endpoint. Use `"all"` to show all retained endpoints.
#' @param overall Logical. If `TRUE`, display only the `All comparators`
#'   aggregate result. This is useful for assessing the overall decision when
#'   several comparisons are required; `display` is ignored in this mode.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object.
#' @export
plot_stability <- function(x, target_power = 0.80, display = "all",
                           overall = FALSE, endpoint = "Total", ...) {
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1)
    stop("'target_power' must be a single number between 0 and 1.")
  if (!is.logical(overall) || length(overall) != 1L || is.na(overall))
    stop("'overall' must be TRUE or FALSE.")
  raw_data <- .diagnostic_long(x)
  if (isTRUE(overall)) {
    if (!any(raw_data$Comparator == "All comparators"))
      stop("'overall = TRUE' requires an 'All comparators' result; use a simulation with multiple comparisons.")
    data <- raw_data[raw_data$Comparator == "All comparators", , drop = FALSE]
  } else {
    data <- .filter_diagnostic_display(raw_data, display)
  }
  data <- .filter_diagnostic_endpoint(data, endpoint)
  data <- data[!is.na(data$value), , drop = FALSE]
  data$trial <- stats::ave(data$value, data$n_total, data$nsim,
                    data$Comparator, data$Endpoint,
                    FUN = seq_along)
  data$cumulative_power <- stats::ave(data$value, data$n_total, data$nsim,
                               data$Comparator, data$Endpoint,
                               FUN = function(z) cumsum(z) / seq_along(z))
  data$Endpoint <- factor(data$Endpoint,
                          levels = c(setdiff(unique(data$Endpoint), "Total"), "Total"))
  ggplot2::ggplot(data, ggplot2::aes(x = trial, y = cumulative_power,
                                     color = Endpoint, group = Endpoint)) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                        color = "#4B5563", linewidth = 0.5) +
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.7) +
    ggplot2::facet_grid(
      n_total + nsim ~ Comparator, scales = "free_x",
      labeller = ggplot2::labeller(
        n_total = function(x) paste0("Total sample size = ", x),
        nsim = function(x) paste0("Simulation trials = ", x)
      )
    ) +
    .diagnostic_endpoint_scale(data) +
    ggplot2::scale_y_continuous(name = "Cumulative achieved power",
                                limits = c(0, 1),
                                labels = scales::label_percent()) +
    ggplot2::labs(
      x = "Simulation trial", color = "Endpoint",
      subtitle = paste0(.diagnostic_sample_size_label(data),
                         "; target power = ", scales::label_percent()(target_power))
    ) +
    .diagnostic_theme()
}

#' Plot Monte Carlo error
#'
#' Displays the 95\% Monte Carlo confidence-interval half-width as the number
#' of simulated trials increases. Decreasing curves indicate improving
#' precision; a curve that has not flattened may require a larger `nsim`.
#' @param x A `simss`, `simpower`, or `simpower_curve` object, or a list of
#'   such objects. A list can compare simulations with different `n` and
#'   `nsim` values.
#' @param display Comparator names to display, or `"all"`.
#' @param endpoint Endpoint names to display. The default, `"Total"`, shows
#'   only the aggregate endpoint. Use `"all"` to show all retained endpoints.
#' @param overall Logical. If `TRUE`, display only the `All comparators`
#'   aggregate result. This is useful for assessing the overall decision when
#'   several comparisons are required; `display` is ignored in this mode.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object.
#' @export
plot_mc_error <- function(x, display = "all", overall = FALSE,
                          endpoint = "Total", ...) {
  dots <- list(...)
  if ("target_error" %in% names(dots))
    stop("'target_error' is no longer available; inspect the achieved final Monte Carlo error instead.")
  if (!is.logical(overall) || length(overall) != 1L || is.na(overall))
    stop("'overall' must be TRUE or FALSE.")
  raw_data <- .diagnostic_long(x)
  if (isTRUE(overall)) {
    if (!any(raw_data$Comparator == "All comparators"))
      stop("'overall = TRUE' requires an 'All comparators' result; use a simulation with multiple comparisons.")
    data <- raw_data[raw_data$Comparator == "All comparators", , drop = FALSE]
  } else {
    data <- .filter_diagnostic_display(raw_data, display)
  }
  data <- .filter_diagnostic_endpoint(data, endpoint)
  data <- data[!is.na(data$value), , drop = FALSE]
  data$trial <- stats::ave(data$value, data$n_total, data$nsim,
                    data$Comparator, data$Endpoint,
                    FUN = seq_along)
  data$mc_error <- stats::ave(data$value, data$n_total, data$nsim,
                       data$Comparator, data$Endpoint, FUN = function(z) {
                         vapply(seq_along(z), function(k) {
                           interval <- stats::binom.test(
                             sum(z[seq_len(k)]), k, conf.level = 0.95
                           )$conf.int
                           (interval[[2L]] - interval[[1L]]) / 2
                         }, numeric(1))
                       })
  group_id <- interaction(data$n_total, data$nsim, data$Comparator,
                          data$Endpoint, drop = TRUE)
  data$last_trial <- data$trial == stats::ave(data$trial, group_id, FUN = max)
  final_data <- data[data$last_trial, , drop = FALSE]
  final_data$label <- scales::label_percent(accuracy = 0.1)(final_data$mc_error)
  data$label <- NA_character_
  data$label[data$last_trial] <- final_data$label
  achieved_title <- paste0(
    "Achieved Monte Carlo error = ",
    paste(unique(final_data$label), collapse = ", "),
    " (at final simulation trial)"
  )
  data$Endpoint <- factor(data$Endpoint,
                          levels = c(setdiff(unique(data$Endpoint), "Total"), "Total"))
  ggplot2::ggplot(data, ggplot2::aes(x = trial, y = mc_error,
                                     color = Endpoint, group = Endpoint)) +
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.7) +
    ggplot2::geom_text(
      data = final_data,
      ggplot2::aes(label = label),
      hjust = -0.1, vjust = -0.5, size = 3, show.legend = FALSE
    ) +
    ggplot2::facet_grid(
      n_total + nsim ~ Comparator, scales = "free_x",
      labeller = ggplot2::labeller(
        n_total = function(x) paste0("Total sample size = ", x),
        nsim = function(x) paste0("Simulation trials = ", x)
      )
    ) +
    ggplot2::scale_y_continuous(name = "Monte Carlo error (95% CI half-width)",
                                labels = scales::label_percent(),
                                expand = ggplot2::expansion(mult = c(0.02, 0.15))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.15))) +
    .diagnostic_endpoint_scale(data) +
    ggplot2::labs(
      title = achieved_title,
      x = "Simulation trial", color = "Endpoint",
      subtitle = .diagnostic_sample_size_label(data)
    ) +
    .diagnostic_theme() +
    ggplot2::coord_cartesian(clip = "off")
}

#' Plot simulated decision heatmaps
#'
#' Shows the equivalence decision (pass or fail) for every endpoint and
#' comparator in every simulated trial. Each tile is one binary decision;
#' green means that the criterion passed and orange means that it failed. The
#' `Total` row is the combined decision for the comparator, and an
#' `All comparators` panel, when present, shows the overall decision. This
#' makes isolated failures, unstable endpoints, and systematic comparator
#' differences easy to identify.
#' @param x A `simss`, `simpower`, or `simpower_curve` object.
#' @param display Comparator names to display, or `"all"`.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object.
#' @details The x-axis is the simulation-trial number and has no scientific
#' meaning beyond showing the sequence of simulated datasets. Endpoint rows
#' show individual decisions; `Total` shows the combined endpoint decision
#' according to `k` and the specified decision rule. The heatmap is a
#' diagnostic of trial-level decisions, not a replacement for the numerical
#' power estimate.
#' @export
plot_decision_heatmap <- function(x, display = "all", ...) {
  data <- .filter_diagnostic_display(.diagnostic_long(x), display)
  data$Decision <- factor(data$value, levels = c(0, 1),
                          labels = c("Fail", "Pass"))
  data$Endpoint <- factor(data$Endpoint,
                          levels = c(setdiff(unique(data$Endpoint), "Total"), "Total"))
  ggplot2::ggplot(data, ggplot2::aes(x = trial, y = Endpoint, fill = Decision)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(
      n_total ~ Comparator, scales = "free_x",
      labeller = ggplot2::labeller(
        n_total = function(x) paste0("Total sample size = ", x)
      )
    ) +
    ggplot2::scale_fill_manual(values = c(Fail = "#D55E00", Pass = "#009E73"),
                               drop = FALSE) +
    ggplot2::labs(x = "Simulation trial", y = "Endpoint", fill = "Decision") +
    .diagnostic_theme()
}
