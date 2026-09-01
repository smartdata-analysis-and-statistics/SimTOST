.type1_boundary_colors <- c(lower = "#F8766D", upper = "#00BFC4")

.type1_error_plot_data <- function(x) {
  objects <- if (inherits(x, "type1error_set")) x else list(x)
  boundaries <- names(objects)
  if (is.null(boundaries)) boundaries <- rep("boundary", length(objects))
  rows <- lapply(seq_along(objects), function(i) {
    object <- objects[[i]]
    if (inherits(object, "simpower_curve")) {
      data.frame(
        boundary = boundaries[[i]], n_total = object$power_curve$n_total,
        estimate = object$power_curve$power,
        Lower = object$power_curve$power_LCI,
        Upper = object$power_curve$power_UCI,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        boundary = boundaries[[i]], n_total = NA_real_,
        estimate = object$type1_error,
        Lower = object$power_LCI, Upper = object$power_UCI,
        stringsAsFactors = FALSE
      )
    }
  })
  do.call(rbind, rows)
}

.type1_error_alpha <- function(x) {
  object <- if (inherits(x, "type1error_set")) x[[1L]] else x
  result <- if (!is.null(object$result)) object$result else object
  alpha <- if (!is.null(result$param.d$alpha)) result$param.d$alpha else
    if (!is.null(object$param.d$alpha)) object$param.d$alpha else 0.05
  as.numeric(alpha[[1L]])
}

.type1_mark_worst <- function(data) {
  data$Worst <- FALSE
  if ("n_total" %in% names(data) && any(!is.na(data$n_total))) {
    groups <- split(seq_len(nrow(data)), data$n_total)
  } else {
    groups <- list(seq_len(nrow(data)))
  }
  for (group in groups) {
    data$Worst[group[[which.max(data$estimate[group])]]] <- TRUE
  }
  data$Worst <- factor(data$Worst, levels = c(FALSE, TRUE))
  data
}

#' Plot empirical Type I error
#'
#' The plot displays the complete-trial Type I error at all evaluated
#' boundaries. The worst-case scenario is highlighted with a thicker interval
#' and a larger point, and its estimated value is shown next to the point.
#' Error bars are 95% Monte Carlo confidence intervals.
#' @param x An object returned by [type1Error()].
#' @param ... Unused additional arguments.
#' @return A `ggplot` object.
#' @export
#' @method plot type1error
plot.type1error <- function(x, ...) {
  objects <- setNames(list(x), x$null_boundary)
  conf.level <- attr(x, "conf.level")
  if (is.null(conf.level)) conf.level <- 0.95
  global_upper <- x$simultaneous_upper
  if (is.null(global_upper)) global_upper <- NA_real_
  objects <- structure(objects, class = "type1error_set",
                       conf.level = conf.level,
                       global_upper = global_upper)
  plot.type1error_set(objects, ...)
}

#' @export
#' @method plot type1error_set
plot.type1error_set <- function(x, ...) {
  if (!inherits(x, "type1error_set"))
    stop("Input must be a type1error_set object.")
  data <- .type1_error_plot_data(x)
  data <- .type1_mark_worst(data)
  data$estimate_label <- sprintf("%.1f%%", 100 * data$estimate)
  worst <- data[data$Worst == TRUE, , drop = FALSE]
  alpha <- .type1_error_alpha(x)
  global_upper <- attr(x, "global_upper")
  if (is.null(global_upper) && inherits(x, "type1error_set")) {
    global_upper <- max(vapply(x, function(object) {
      if (!is.null(object$simultaneous_upper)) object$simultaneous_upper else NA_real_
    }, numeric(1)), na.rm = TRUE)
    if (!is.finite(global_upper)) global_upper <- NA_real_
  }
  is_curve <- any(!is.na(data$n_total))
  if (is_curve) {
    p <- ggplot2::ggplot(
      data, ggplot2::aes(x = n_total, y = estimate,
                         color = boundary, group = boundary)
    ) +
      ggplot2::geom_hline(yintercept = alpha, linetype = "dashed",
                          color = "#343A40", linewidth = 0.6) +
      ggplot2::geom_line(ggplot2::aes(linewidth = Worst), alpha = 0.7) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = Lower, ymax = Upper, linewidth = Worst),
        alpha = 0.7
      ) +
      ggplot2::geom_point(ggplot2::aes(size = Worst), alpha = 0.7) +
      ggplot2::geom_point(data = worst, shape = 21, fill = NA,
                          color = "black", stroke = 1.1, size = 3.4) +
      ggplot2::geom_text(
        data = worst,
        ggplot2::aes(label = estimate_label),
        vjust = -0.8, show.legend = FALSE
      ) +
      ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.7, `TRUE` = 1.6), guide = "none") +
      ggplot2::scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.0), guide = "none") +
      ggplot2::scale_x_continuous(name = "Total sample size")
  } else {
    p <- ggplot2::ggplot(
      data, ggplot2::aes(x = boundary, y = estimate)
    ) +
      ggplot2::geom_hline(yintercept = alpha, linetype = "dashed",
                          color = "#343A40", linewidth = 0.6) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = Lower, ymax = Upper, color = boundary,
                     linewidth = Worst), alpha = 0.7
      ) +
      ggplot2::geom_point(ggplot2::aes(color = boundary, size = Worst), alpha = 0.7) +
      ggplot2::geom_point(data = worst, shape = 21, fill = NA,
                          color = "black", stroke = 1.1, size = 3.6) +
      ggplot2::geom_text(
        data = worst,
        ggplot2::aes(label = estimate_label),
        vjust = -0.8, show.legend = FALSE
      ) +
      ggplot2::scale_linewidth_manual(values = c(`FALSE` = 1.0, `TRUE` = 2.0), guide = "none") +
      ggplot2::scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.2), guide = "none") +
      ggplot2::scale_color_manual(values = .type1_boundary_colors, drop = FALSE) +
      ggplot2::scale_x_discrete(name = "Least-favorable null boundary")
  }
  p <- if (is_curve) {
    p + ggplot2::scale_y_continuous(
      name = "Empirical Type I error", limits = c(0, 1),
      labels = scales::label_percent(accuracy = 1)
    )
  } else {
    y_max <- max(c(alpha, data$Upper), na.rm = TRUE)
    y_max <- max(y_max * 1.15, y_max + 0.005)
    p + ggplot2::scale_y_continuous(
      name = "Empirical Type I error", limits = c(0, y_max),
      labels = scales::label_percent(accuracy = 1)
    )
  }
  p +
    ggplot2::labs(subtitle = paste0(
      "Larger point and thicker interval indicate the worst-case scenario")) +
    .diagnostic_theme()
}
