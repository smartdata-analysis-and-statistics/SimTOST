#' Summarize fixed-sample-size power results
#'
#' @param object An object returned by [simPower()].
#' @param ... Unused additional arguments.
#' @return A data frame containing the distribution, design, comparison,
#' estimand, hypotheses, sample size, estimated power, and Monte Carlo
#' confidence interval.
#' @export
#' @method summary simpower
summary.simpower <- function(object, ...) {
  if (inherits(object, "countpower")) {
    out <- summary.countpower(object, ...)
    out$distribution <- object$distribution
    return(out[, c("distribution", setdiff(names(out), "distribution")),
               drop = FALSE])
  }
  curve <- inherits(object, "simpower_curve")
  results <- if (curve) object$curve_results else list(object)
  first <- results[[1L]]
  continuous_result <- first$result
  param_d <- if (!is.null(continuous_result)) continuous_result$param.d else NULL
  comparisons <- if (!is.null(continuous_result)) {
    vapply(continuous_result$param$list_comparator,
           paste, collapse = "_vs_", FUN.VALUE = character(1))
  } else {
    NA_character_
  }
  estimand <- if (identical(object$distribution, "norm") &&
                 identical(param_d$ctype, "DOM")) {
    "mean difference"
  } else if (identical(object$distribution, "lnorm")) {
    "arithmetic mean ratio"
  } else {
    "mean ratio"
  }
  power_rows <- lapply(results, function(result) {
    response <- result$result$response
    row <- list(
      n = if (!is.null(result$n)) result$n[[1L]] else
        if (!is.null(result$n_per_arm)) result$n_per_arm else NA_real_
    )
    row$n_total <- if (!is.null(result$n_total)) result$n_total else
      if (!is.null(response)) response$n_total[[1L]] else NA_real_
    row[["Achieved power"]] <- result$power
    row$power_LCI <- result$power_LCI
    row$power_UCI <- result$power_UCI
    as.data.frame(row, check.names = FALSE)
  })
  out <- as.data.frame(data.table::rbindlist(power_rows, fill = TRUE))
  out <- out[, c("n", "n_total", "Achieved power", "power_LCI", "power_UCI"),
             drop = FALSE]

  cat("Fixed Sample Power Summary\n")
  cat(strrep("-", 26), "\n")
  cat("Design type        :", if (is.null(param_d)) first$design else param_d$dtype, "\n")
  cat("Distribution       :", object$distribution, "\n")
  cat("Comparison         :", paste(comparisons, collapse = "; "), "\n")
  cat("Estimand           :", estimand, "\n")
  cat("Hypotheses         : H0: estimand <= L or >= U; H1: L < estimand < U\n")
  cat("\nEstimated Achieved Power:\n")
  print(out, row.names = FALSE)
  invisible(out)
}

#' Extract the Monte Carlo confidence interval from fixed-sample-size power
#' results
#'
#' @param object An object returned by [simPower()].
#' @param parm Unused; included for compatibility with [stats::confint()].
#' @param level Confidence level for the Monte Carlo interval.
#' @param ... Unused additional arguments.
#' @return A two-column matrix containing the lower and upper interval limits.
#' @export
#' @method confint simpower
confint.simpower <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number between 0 and 1.")
  if (inherits(object, "simpower_curve")) {
    return(data.frame(
      n = object$power_curve$n,
      n_total = object$power_curve$n_total,
      power = object$power_curve$power,
      Lower = object$power_curve$power_LCI,
      Upper = object$power_curve$power_UCI
    ))
  }
  if (!is.null(object$successes) && !is.null(object$nsim)) {
    interval <- stats::binom.test(object$successes, object$nsim,
                                   conf.level = level)$conf.int
  } else {
    if (!isTRUE(all.equal(level, 0.95)))
      stop("Only the stored 95% interval is available for this result.")
    interval <- c(object$power_LCI, object$power_UCI)
  }
  data.frame(
    n = if (!is.null(object$n)) object$n[[1L]] else NA_real_,
    n_total = if (!is.null(object$n_total)) object$n_total else NA_real_,
    power = object$power,
    Lower = interval[[1L]],
    Upper = interval[[2L]]
  )
}

#' Plot fixed-sample-size power results
#'
#' @param x An object returned by [simPower()].
#' @param target_power Target power shown as a horizontal reference line. The
#'   default is 0.80, unless the object stores a planning target, in which
#'   case that target is used when this argument is omitted. An explicit value
#'   always overrides the stored target.
#' @param display Character vector of comparator panels to display. Use
#'   `"all"` (the default) to display every comparator, or provide comparator
#'   names such as `"R_vs_T"` to display selected panels.
#' @param all Logical. If `TRUE` (the default), display only the aggregate
#'   `All comparators` result without a comparator facet. If `FALSE`, display
#'   the current comparator panels selected by `display`.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object showing estimated power and its Monte Carlo
#' confidence interval as a clearly visible vertical line.
#' @export
#' @method plot simpower_curve
plot.simpower_curve <- function(x, target_power = 0.80, display = "all", all = TRUE, ...) {
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1)
    stop("'target_power' must be a single number between 0 and 1.")

  curve <- x$power_curve
  if (is.null(curve) || nrow(curve) == 0L)
    stop("The simpower curve does not contain power results.")
  if (!is.character(display) || length(display) == 0L || anyNA(display))
    stop("'display' must be 'all' or a non-empty character vector of comparator names.")
  if (!is.logical(all) || length(all) != 1L || is.na(all))
    stop("'all' must be a single logical value.")
  display <- unique(gsub(" vs ", "_vs_", display, fixed = TRUE))

  first_result <- x$curve_results[[1L]]
  n_comparators <- if (!is.null(first_result$result) &&
                       !is.null(first_result$result$param$list_comparator)) {
    length(first_result$result$param$list_comparator)
  } else if (!is.null(first_result$comparisons)) {
    length(first_result$comparisons)
  } else {
    1L
  }
  detail <- list()
  if (n_comparators > 1L) {
    detail[[1L]] <- data.frame(
      n = curve$n_total, power = curve$power,
      power_LCI = curve$power_LCI, power_UCI = curve$power_UCI,
      Comparator = "All comparators", Endpoint = "Total",
      stringsAsFactors = FALSE
    )
  }

  power_interval <- function(successes, nsim) {
    interval <- stats::prop.test(
      x = successes, n = nsim, correct = TRUE
    )$conf
    c(interval[[1L]], interval[[2L]])
  }

  # Continuous scalar results retain the simulated trial table, allowing the
  # curve to show endpoint- and comparator-specific power as well as the
  # simultaneous all-comparisons power.
  for (i in seq_along(x$curve_results)) {
    scalar <- x$curve_results[[i]]
    if (is.null(scalar$result) || is.null(scalar$result$table.test)) next
    test_table <- scalar$result$table.test

    comp_cols <- grep("^totalyComp:", names(test_table), value = TRUE)
    endpoint_cols <- grep("Comp:", names(test_table), value = TRUE)
    endpoint_cols <- endpoint_cols[
      !grepl("^(totaly|mu_|sd_|eql_|equ_)", endpoint_cols)
    ]

    for (comp_col in comp_cols) {
      comparator <- sub("^totalyComp:", "", comp_col)
      comparator <- gsub(" vs ", "_vs_", comparator, fixed = TRUE)
      comparator_power <- mean(test_table[[comp_col]], na.rm = TRUE)
      comparator_ci <- power_interval(
        sum(test_table[[comp_col]], na.rm = TRUE), nrow(test_table)
      )
      detail[[length(detail) + 1L]] <- data.frame(
        n = curve$n_total[[i]], power = comparator_power,
        power_LCI = comparator_ci[[1L]], power_UCI = comparator_ci[[2L]],
        Comparator = comparator, Endpoint = "Total",
        stringsAsFactors = FALSE
      )

      for (endpoint_col in endpoint_cols) {
        endpoint_comparator <- sub(".*Comp:", "", endpoint_col)
        if (endpoint_comparator != sub("^totalyComp:", "", comp_col))
          next
        endpoint_power <- mean(test_table[[endpoint_col]], na.rm = TRUE)
        endpoint_ci <- power_interval(
          sum(test_table[[endpoint_col]], na.rm = TRUE), nrow(test_table)
        )
        detail[[length(detail) + 1L]] <- data.frame(
          n = curve$n_total[[i]], power = endpoint_power,
          power_LCI = endpoint_ci[[1L]], power_UCI = endpoint_ci[[2L]],
          Comparator = comparator,
          Endpoint = sub("Comp:.*", "", endpoint_col),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(detail) == 0L) {
    detail[[1L]] <- data.frame(
      n = curve$n_total, power = curve$power,
      power_LCI = curve$power_LCI, power_UCI = curve$power_UCI,
      Comparator = "Comparison", Endpoint = "Total",
      stringsAsFactors = FALSE
    )
  }

  plotdata <- do.call(rbind, detail)
  plotdata$Comparator <- gsub(" vs ", "_vs_", plotdata$Comparator, fixed = TRUE)
  if (all) {
    if (any(plotdata$Comparator == "All comparators")) {
      plotdata <- plotdata[plotdata$Comparator == "All comparators", , drop = FALSE]
    } else {
      plotdata$Comparator <- "All comparators"
    }
  }
  available <- unique(plotdata$Comparator)
  if (!any(display == "all")) {
    unknown <- setdiff(display, available)
    if (length(unknown))
      stop("Unknown comparator(s) in 'display': ", paste(unknown, collapse = ", "))
    plotdata <- plotdata[plotdata$Comparator %in% display, , drop = FALSE]
    if (nrow(plotdata) == 0L)
      stop("'display' did not select any comparator panels.")
  }
  plotdata$Endpoint <- factor(
    plotdata$Endpoint,
    levels = c(setdiff(unique(plotdata$Endpoint), "Total"), "Total")
  )
  colors <- stats::setNames(
    grDevices::hcl.colors(length(setdiff(levels(plotdata$Endpoint), "Total")),
                          palette = "Dark 3"),
    setdiff(levels(plotdata$Endpoint), "Total")
  )
  colors <- c(colors, Total = "black")

  p <- ggplot2::ggplot(
    plotdata,
    ggplot2::aes(x = n, y = power, color = Endpoint, group = Endpoint)
  ) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                        linewidth = 0.6, color = "#343A40") +
    ggplot2::geom_line(linewidth = 1.0, alpha = 0.6, na.rm = TRUE) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = power_LCI, ymax = power_UCI),
      linewidth = 0.8, color = "black", alpha = 1, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 1.25, alpha = 0.75, na.rm = TRUE
    ) +
    (if (all) ggplot2::facet_null() else ggplot2::facet_grid(. ~ Comparator, scales = "free_x")) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      name = "Estimated power (%)", limits = c(0, 1.12),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::label_percent(accuracy = 1)
    ) +
    ggplot2::scale_x_continuous(
      name = "Total sample size",
      labels = scales::label_number(accuracy = 1, big.mark = ",")
    ) +
    ggplot2::labs(
      color = "Endpoint"
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#E5E7EB",
                                                 linewidth = 0.35),
      panel.grid.major.y = ggplot2::element_line(color = "#D9DDE3",
                                                 linewidth = 0.45),
      panel.spacing.x = grid::unit(1.35, "lines"),
      strip.background = ggplot2::element_rect(fill = "#EAF1F8",
                                                color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 11,
                                         color = "#20364D"),
      axis.title = ggplot2::element_text(face = "bold", color = "#20364D"),
      axis.text = ggplot2::element_text(color = "#343A40"),
      legend.position = "none",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width = grid::unit(1.7, "lines"),
      legend.text = ggplot2::element_text(size = 10.5),
      plot.margin = ggplot2::margin(8, 10, 8, 10)
    )
  p
}

#' Plot fixed-sample-size power results
#'
#' @param x An object returned by [simPower()].
#' @param target_power Target power shown as a horizontal reference line.
#' @param display Character vector of comparator panels to display. This
#'   argument is retained for compatibility with curve results.
#' @param all Logical retained for compatibility with curve results. If
#'   `TRUE`, the aggregate result is displayed when available.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object showing estimated power and its Monte Carlo
#'   confidence interval.
#' @export
#' @method plot simpower
plot.simpower <- function(x, target_power = 0.80, display = "all", all = TRUE, ...) {
  if (inherits(x, "simpower_curve"))
    return(plot.simpower_curve(x, target_power = target_power,
                               display = display, all = all, ...))
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1)
    stop("'target_power' must be a single number between 0 and 1.")
  if (!is.logical(all) || length(all) != 1L || is.na(all))
    stop("'all' must be a single logical value.")
  data <- data.frame(
    n = x$n,
    power = x$power,
    power_LCI = x$power_LCI,
    power_UCI = x$power_UCI
  )
  p <- ggplot2::ggplot(data, ggplot2::aes(x = n, y = power)) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                        color = "#4B5563") +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = power_LCI, ymax = power_UCI),
      linewidth = 1.1, color = "#0072B2", alpha = 0.9
    ) +
    ggplot2::geom_point(
      size = 1.25, color = "#0072B2", alpha = 0.75
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1.12),
                                name = "Estimated power (%)",
                                labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(x = "Sample size per arm") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", color = "#20364D"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (length(unique(data$n)) == 1L) {
    span <- max(1, 0.30 * abs(data$n[[1]]))
    p <- p + ggplot2::coord_cartesian(
      xlim = c(max(0, data$n[[1]] - span), data$n[[1]] + span)
    ) + ggplot2::scale_x_continuous(
      breaks = data$n[[1]],
      labels = scales::label_number(accuracy = 1, big.mark = ",")
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      labels = scales::label_number(accuracy = 1, big.mark = ",")
    )
  }
  p
}

#' Plot count-outcome power results
#'
#' @param x An object returned by [simPower()].
#' @param target_power Target power shown as a horizontal reference line. The
#'   default is 0.80, unless the object stores a planning target, in which
#'   case that target is used when this argument is omitted. An explicit value
#'   always overrides the stored target.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object showing estimated power and its Monte Carlo
#' confidence interval. The interval is shown as a clearly visible vertical
#' line through each point.
#' @export
#' @method plot countpower
plot.countpower <- function(x, target_power = 0.80, ...) {
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1)
    stop("'target_power' must be a single number between 0 and 1.")
  data <- data.frame(
    n = x$n_total,
    power = x$power,
    power_LCI = x$power_LCI,
    power_UCI = x$power_UCI
  )
  ggplot2::ggplot(data, ggplot2::aes(x = n, y = power)) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                        color = "#4B5563") +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = power_LCI, ymax = power_UCI),
      linewidth = 0.65, fatten = 1.55, color = "#0072B2"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1.12),
                                name = "Estimated power (%)",
                                labels = scales::label_percent(accuracy = 1)) +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 1,
                                                              big.mark = ",")) +
    ggplot2::labs(x = "Total sample size", title = "Count-outcome power") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot count-outcome sample-size results
#'
#' @param x An object returned by [sampleSize()]. Count-specific dispatch is
#' retained through the secondary compatibility class.
#' @param target_power Target power shown as a horizontal reference line. The
#'   default is 0.80, unless the object stores a planning target, in which
#'   case that target is used when this argument is omitted. An explicit value
#'   always overrides the stored target.
#' @param display Character vector of comparator panels to display. Use
#'   `"all"` (the default) to display every comparator and, for joint count
#'   analyses, the all-comparisons panel.
#' @param all Logical. If `TRUE` (the default), display only the aggregate
#'   `All comparators` result without a comparator facet. If `FALSE`, display
#'   the comparator panels selected by `display`.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object showing the simulated power curve over the
#' evaluated candidate sample sizes, with the selected sample size highlighted.
#' Confidence intervals are shown as clearly visible vertical lines through
#' each point.
#' If an older `countss` object has no search history, the plot falls back to
#' the selected-sample-size point.
#' @export
#' @method plot countss
plot.countss <- function(x, target_power = 0.80, display = "all", all = TRUE, ...) {
  if (missing(target_power))
    target_power <- .stored_plot_target_power(x, fallback = target_power)
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1)
    stop("'target_power' must be a single number between 0 and 1.")
  if (!is.logical(all) || length(all) != 1L || is.na(all))
    stop("'all' must be a single logical value.")
  history <- x$table.test
  has_history <- is.data.frame(history) && nrow(history) > 0L
  decision_cols <- if (has_history)
    grep("Comp:", names(history), value = TRUE) else character()

  if (length(decision_cols)) {
    comparator_cols <- grep("^totalyComp:", decision_cols, value = TRUE)
    n_comparators <- length(comparator_cols)
    if (n_comparators > 1L && "totaly" %in% names(history))
      decision_cols <- c(decision_cols, "totaly")

    groups <- split(seq_len(nrow(history)), history$n_total)
    power_interval <- function(successes, trials) {
      interval <- suppressWarnings(
        stats::prop.test(successes, trials, correct = TRUE)$conf
      )
      c(successes / trials, interval[[1L]], interval[[2L]])
    }
    plotdata <- do.call(rbind, lapply(decision_cols, function(column) {
      do.call(rbind, lapply(groups, function(index) {
        interval <- power_interval(
          sum(history[index, column], na.rm = TRUE), length(index)
        )
        if (identical(column, "totaly")) {
          endpoint <- "Total"
          comparator <- "All comparators"
        } else {
          endpoint <- sub("Comp:.*", "", column)
          comparator <- sub(".*Comp:", "", column)
          if (identical(endpoint, "totaly")) endpoint <- "Total"
        }
        data.frame(
          n = as.numeric(history$n_total[index[[1L]]]),
          power = interval[[1L]], power_LCI = interval[[2L]],
          power_UCI = interval[[3L]], Endpoint = endpoint,
          Comparator = comparator, stringsAsFactors = FALSE
        )
      }))
    }))
    rownames(plotdata) <- NULL
    plotdata <- plotdata[order(plotdata$Comparator, plotdata$Endpoint,
                               plotdata$n), , drop = FALSE]

    selected_n <- .select_power_plot_n(
      plotdata, target_power = target_power, n_col = "n",
      fallback = x$n_total
    )

    if (!is.character(display) || length(display) == 0L || anyNA(display))
      stop("'display' must be 'all' or a non-empty character vector of comparator names.")
    display <- unique(gsub(" vs ", "_vs_", display, fixed = TRUE))
    plotdata$Comparator <- gsub(" vs ", "_vs_", plotdata$Comparator,
                                fixed = TRUE)
    if (all) {
      if (any(plotdata$Comparator == "All comparators")) {
        plotdata <- plotdata[plotdata$Comparator == "All comparators", , drop = FALSE]
      } else {
        plotdata$Comparator <- "All comparators"
      }
    }
    available <- unique(plotdata$Comparator)
    if (!any(display == "all")) {
      unknown <- setdiff(display, available)
      if (length(unknown))
        stop("Unknown comparator(s) in 'display': ", paste(unknown, collapse = ", "))
      plotdata <- plotdata[plotdata$Comparator %in% display, , drop = FALSE]
    }
    plotdata$Endpoint <- factor(
      plotdata$Endpoint,
      levels = c(setdiff(unique(plotdata$Endpoint), "Total"), "Total")
    )
    comparator_levels <- unique(plotdata$Comparator)
    plotdata$Comparator <- factor(plotdata$Comparator,
                                  levels = comparator_levels)
    endpoint_levels <- setdiff(levels(plotdata$Endpoint), "Total")
    endpoint_colors <- stats::setNames(
      grDevices::hcl.colors(length(endpoint_levels), palette = "Dark 3"),
      endpoint_levels
    )
    endpoint_colors <- c(endpoint_colors, Total = "black")
    selected <- plotdata[plotdata$n == selected_n, , drop = FALSE]
    errorbar_width <- if (length(unique(plotdata$n)) > 1L) {
      max(1, 0.01 * diff(range(plotdata$n, na.rm = TRUE)))
    } else {
      0
    }
    return(
      ggplot2::ggplot(
        plotdata,
        ggplot2::aes(x = n, y = power, color = Endpoint, group = Endpoint)
      ) +
        ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                            linewidth = 0.6, color = "#343A40") +
        ggplot2::geom_line(linewidth = 1.0, alpha = 0.6, na.rm = TRUE) +
        ggplot2::geom_point(size = 1.25, alpha = 0.75, na.rm = TRUE) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = power_LCI, ymax = power_UCI),
          width = max(0.4, errorbar_width * 0.35), linewidth = 0.8, color = "black", alpha = 1,
          na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = selected,
          ggplot2::aes(x = n, y = power),
          inherit.aes = FALSE, size = 2.2, shape = 21,
          fill = "#009E73", color = "#1F2937", stroke = 0.6
        ) +
        (if (all) ggplot2::facet_null() else ggplot2::facet_grid(. ~ Comparator, scales = "free_x")) +
        ggplot2::scale_color_manual(values = endpoint_colors, drop = FALSE) +
        ggplot2::scale_y_continuous(
          name = "Estimated power (%)", limits = c(0, 1.12),
          breaks = seq(0, 1, by = 0.2),
          labels = scales::label_percent(accuracy = 1)
        ) +
        ggplot2::scale_x_continuous(
          name = "Total sample size",
          labels = scales::label_number(accuracy = 1, big.mark = ",")
        ) +
        ggplot2::labs(
          color = "Endpoint", title = "Count-outcome sample size",
          subtitle = "Vertical bars show 95% Monte Carlo confidence intervals"
        ) +
        ggplot2::theme_minimal(base_size = 11.5) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(fill = "#EAF1F8",
                                                    color = NA),
          strip.text = ggplot2::element_text(face = "bold", color = "#20364D"),
          axis.title = ggplot2::element_text(face = "bold", color = "#20364D"),
          legend.position = "none"
        )
    )
  }

  # Compatibility fallback for countss objects created before component
  # decisions were retained.
  data <- data.frame(
    n = x$n_total, power = x$power,
    power_LCI = x$power_LCI, power_UCI = x$power_UCI
  )
  ggplot2::ggplot(data, ggplot2::aes(x = n, y = power)) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                        color = "#4B5563") +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = power_LCI, ymax = power_UCI),
      linewidth = 1.1, size = 1.55, color = "#009E73"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1.12),
                                name = "Achieved power (%)",
                                labels = scales::label_percent(accuracy = 1)) +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 1,
                                                              big.mark = ",")) +
    ggplot2::labs(x = "Total sample size", title = "Count-outcome sample size") +
    ggplot2::theme_minimal(base_size = 12)
}
