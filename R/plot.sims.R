#' @title Plot Power vs Sample Size for Simulation Results
#' @description Generates a detailed plot showing the relationship between power and total sample size for each comparator and the overall combined comparators. The combined-comparator panel is omitted when only one comparator is present.
#' The plot also includes confidence intervals for power estimates and highlights the target power with a dashed line for easy visual comparison.
#'
#' @param x An object of class `simss` containing simulation results.
#' @param target_power Target power shown as a horizontal reference line. The
#'   default is 0.80, unless the object stores a planning target, in which
#'   case that target is used when this argument is omitted. An explicit value
#'   always overrides the stored target.
#' @param display Character vector of comparator panels to display. Use
#'   `"all"` (the default) to display every comparator, or provide comparator
#'   names such as `"R_vs_T"` to display selected panels.
#' @param all Logical. If `TRUE` (the default), display only the aggregate
#'   `All comparators` result without a comparator facet. If `FALSE`, display
#'   the comparator panels selected by `display`.
#' @param \ldots Additional arguments to be passed to the `plot.simss` function for customization.
#'
#' @return A `ggplot` object illustrating:
#'   - Power (y-axis) vs. Total Sample Size (x-axis) for individual endpoints and comparators.
#'   - Clearly visible vertical lines representing the 95% confidence interval
#'     of the power estimates.
#'   - A dashed horizontal line indicating the target power for comparison.
#'   - A larger outlined point at the selected sample size.
#'   - Faceted panels for each comparator, making it easy to compare results across different groups.
#'
#' @details
#' The plot dynamically adjusts to exclude unnecessary components, such as redundant endpoints or comparators with insufficient data, ensuring clarity and simplicity.
#' The `ggplot2` framework is used for visualizations, allowing further customization if needed.
#'
#' @author
#' Johanna Muñoz \email{johanna.munoz@fromdatatowisdom.com}
#'
#' @importFrom data.table .SD
#' @export plot.simss
#'
#' @export
plot.simss <- function(x, target_power = 0.80, display = "all", all = TRUE, ...){
  if (missing(target_power))
    target_power <- .stored_plot_target_power(x, fallback = target_power)
  if (!is.numeric(target_power) || length(target_power) != 1L ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1) {
    stop("'target_power' must be a single number between 0 and 1.")
  }
  if (inherits(x, "countss")) {
    return(plot.countss(x, target_power = target_power,
                        display = display, all = all, ...))
  }

  qnam = n_iter = n_total = t_true = power = Endpoint = power_LCI = power_UCI = NULL # due to NSE notes in R CMD check

  table.test <- x$table.test
  nsim <- as.numeric(x$param.d["nsim"])
  # Calculate totaly test across all comparators= power
  qnam <- colnames(table.test)[grep("^totaly",colnames(table.test))]


  # Get a summary across all the n_iter with confidence interval power
  namexc <- c(colnames(table.test)[grep("^[^(mu_|sd_|eql_|equ_|n_)]",colnames(table.test))],"n_total")
  summary <- table.test[
    , c(lapply(.SD, FUN = function(x) sum(x, na.rm = TRUE)),
        list(n_trials = .N)),
    by = n_total
  ][, c(namexc, "n_trials"), with = FALSE]
  plotdata <- data.table::melt(summary,
                               id.vars = c("t_true", "n_total", "n_trials"))
  powerfun <- function(x, n_trials) {
    bin_test <- stats::prop.test(x = x, n = n_trials, correct = TRUE)
    c(bin_test$estimate[[1]],bin_test$conf[1],bin_test$conf[2])
  }

  powerv <- do.call(rbind, lapply(seq_along(plotdata$value), function(i) {
    powerfun(plotdata$value[[i]], plotdata$n_trials[[i]])
  }))
  colnames(powerv) <- c("power","power_LCI","power_UCI")
  plotdata <- cbind(plotdata,powerv)
  plotdata$Endpoint <- sub("(.*)Comp:.*", "\\1", plotdata$variable)
  plotdata$Comparator <- sub(".*Comp:", "", plotdata$variable)

  plotdata$Endpoint<-ifelse(plotdata$Endpoint=="totaly","Total",plotdata$Endpoint)
  plotdata$Comparator<-ifelse(plotdata$Comparator=="totaly","All comparators",plotdata$Comparator)

  # The optimizer does not necessarily evaluate candidates in increasing
  # order. Sort before drawing lines so that the plotted curve represents
  # power as a function of sample size rather than search order.
  plotdata <- plotdata[order(plotdata$Comparator, plotdata$Endpoint,
                             plotdata$n_total), , drop = FALSE]

  fallback_n <- if (!is.null(x$response) &&
                    !is.null(x$response[["n_total"]])) {
    as.numeric(x$response[["n_total"]][[1L]])
  } else if (!is.null(x$n_total)) {
    as.numeric(x$n_total[[1L]])
  } else {
    NA_real_
  }
  selected_n <- .select_power_plot_n(
    plotdata, target_power = target_power, n_col = "n_total",
    fallback = fallback_n
  )

  # Comparator names are stored with spaces internally but are displayed with
  # underscores. Accept either spelling in `display`.
  if (!is.character(display) || length(display) == 0L || anyNA(display))
    stop("'display' must be 'all' or a non-empty character vector of comparator names.")
  if (!is.logical(all) || length(all) != 1L || is.na(all))
    stop("'all' must be a single logical value.")
  display <- unique(gsub(" vs ", "_vs_", display, fixed = TRUE))

  if (all) {
    if (any(plotdata$Comparator == "All comparators")) {
      plotdata <- plotdata[plotdata$Comparator == "All comparators", , drop = FALSE]
    } else if (length(x$param$list_comparator) > 1L) {
      stop("'all = TRUE' requires an 'All comparators' result; use a simulation with multiple comparisons.")
    } else {
      plotdata$Comparator <- "All comparators"
    }
  }

  # The overall product of comparator results is redundant for a single
  # comparator and is omitted from the plot in that case.
  if (length(x$param$list_comparator) == 1L) {
    if (!all) plotdata <- plotdata[Comparator != "All comparators"]
  }
  plotdata$Comparator <- gsub(" vs ", "_vs_", plotdata$Comparator, fixed = TRUE)
  available <- unique(plotdata$Comparator)
  if (!any(display == "all")) {
    unknown <- setdiff(display, available)
    if (length(unknown))
      stop("Unknown comparator(s) in 'display': ", paste(unknown, collapse = ", "))
    plotdata <- plotdata[Comparator %in% display]
  }

  endpoint_levels <- setdiff(unique(plotdata$Endpoint), "Total")
  endpoint_colors <- stats::setNames(
    grDevices::hcl.colors(length(endpoint_levels), palette = "Dark 3"),
    endpoint_levels
  )
  endpoint_colors <- c(endpoint_colors, Total = "black")

  plotdata$Endpoint <- factor(
    plotdata$Endpoint,
    levels = c(endpoint_levels, "Total")
  )

  comparator_levels <- unique(plotdata$Comparator)
  plotdata$Comparator <- factor(plotdata$Comparator, levels = comparator_levels)
  selected <- plotdata[!is.na(selected_n) &
                         plotdata$n_total == selected_n, , drop = FALSE]
  errorbar_width <- if (length(unique(plotdata$n_total)) > 1L) {
    max(1.5, 0.01 * diff(range(plotdata$n_total, na.rm = TRUE)))
  } else {
    0
  }

  plot <- ggplot2::ggplot(plotdata,
                          ggplot2::aes(x = n_total, y = power,
                                       color = Endpoint, group = Endpoint))+
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                        linewidth = 0.6, color = "#343A40")+
    ggplot2::geom_line(linewidth = 1.05, lineend = "round", alpha = 0.6) +
    ggplot2::geom_point(size = 1.25, alpha = 0.75, na.rm = TRUE) +
    ggplot2::geom_point(
      data = selected,
      ggplot2::aes(x = n_total, y = power),
      inherit.aes = FALSE, size = 2.2, shape = 21,
      fill = "#009E73", color = "#1F2937", stroke = 0.6
    ) +
    # Repeat the interval on top of the curve/points so short intervals remain
    # visible in the aggregate continuous-outcome plot.
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = power_LCI, ymax = power_UCI),
      linewidth = 2.0, color = "black", alpha = 1, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = selected,
      ggplot2::aes(x = n_total, y = power),
      inherit.aes = FALSE, size = 1.9, shape = 21,
      fill = "#009E73", color = "#1F2937", stroke = 0.6
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = power_LCI, ymax = power_UCI),
      width = max(0.4, errorbar_width * 0.35), linewidth = 2.0, color = "black", alpha = 1,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 1.25, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(
      data = selected,
      ggplot2::aes(x = n_total, y = power),
      inherit.aes = FALSE, size = 2.2, shape = 21,
      fill = "#009E73", color = "#1F2937", stroke = 0.6
    ) +
    (if (all) ggplot2::facet_null() else ggplot2::facet_grid(. ~ Comparator, scales = "free_x")) +
    ggplot2::scale_color_manual(values = endpoint_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(name = "Estimated power (%)",
                                breaks = seq(0, 1, by = 0.2),
                                labels = scales::label_percent(accuracy = 1),
                                limits = c(0, 1.12),
                                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::scale_x_continuous(name = "Total sample size",
                                labels = scales::label_number(accuracy = 1,
                                                              big.mark = ","),
                                expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    ggplot2::labs(color = "Endpoint",
                  title = "Continuous-outcome sample size",
                  subtitle = "Vertical bars show 95% Monte Carlo confidence intervals") +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.35),
      panel.grid.major.y = ggplot2::element_line(color = "#D9DDE3", linewidth = 0.45),
      panel.spacing.x = grid::unit(1.35, "lines"),
      strip.background = ggplot2::element_rect(fill = "#EAF1F8", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 11, color = "#20364D"),
      axis.title = ggplot2::element_text(face = "bold", color = "#20364D"),
      axis.text = ggplot2::element_text(color = "#343A40"),
      legend.position = "none",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width = grid::unit(1.7, "lines"),
      legend.text = ggplot2::element_text(size = 10.5),
      plot.margin = ggplot2::margin(8, 10, 8, 10)
    )

  return(plot)
}
