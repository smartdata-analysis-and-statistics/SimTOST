#' @export
#' @method print type1error_joint
print.type1error_joint <- function(x, ...) {
  conf.level <- if (is.null(x$conf.level)) 0.95 else x$conf.level
  global_upper <- if (is.null(x$global_upper)) NA_real_ else x$global_upper
  cat("Joint empirical Type I error analysis\n")
  cat("Boundary configurations:", nrow(x$table), "\n")
  print(x$table[, c("Boundary", "Comparator", "Endpoint", "NullCount",
                    "Type1_Error", "Lower", "Upper",
                    "Simultaneous_Upper")], row.names = FALSE)
  global <- if (!is.null(x$global)) x$global else
    x$table[which.max(x$table$Type1_Error), , drop = FALSE]
  cat(sprintf(
    "Global worst-case Type I error: %.4f [%.4f, %.4f]\n",
    global$Type1_Error[[1L]], global$Lower[[1L]], global$Upper[[1L]]
  ))
  cat(sprintf(
    "Simultaneous %.1f%% upper bound for global Type I error: %.4f\n",
    100 * conf.level, global_upper
  ))
  cat(sprintf("Global scenario: %s, %s, %s boundary\n",
              global$Comparator[[1L]], global$Endpoint[[1L]],
              global$Boundary[[1L]]))
  invisible(x)
}

#' Plot joint Type I error scenarios
#'
#' Displays the complete-trial Type I error for the minimal required boundary
#' scenarios by default. The worst-case displayed scenario is highlighted with
#' a thicker interval and a larger point, and its estimated value is shown
#' directly above the point. The full scenario table, including larger null
#' counts, remains available in the returned object and can be displayed with
#' `null_count = "all"`.
#' The dashed line is the nominal alpha level. The dotted line is the
#' simultaneous one-sided Monte Carlo upper bound for the maximum.
#' @param x An object returned by `type1Error(..., joint = TRUE)`.
#' @param null_count Character string. The default, `"minimal"`, plots only
#'   the minimal required number of boundary endpoints (`m-k+1`) for each
#'   comparator. Use `"all"` to plot every boundary count evaluated by
#'   `type1Error()`.
#' @param ... Unused additional arguments.
#' @return A `ggplot` object.
#' @export
#' @method plot type1error_joint
plot.type1error_joint <- function(x, null_count = c("minimal", "all"), ...) {
  if (!inherits(x, "type1error_joint"))
    stop("Input must be a type1error_joint object.")
  null_count <- match.arg(null_count)
  conf.level <- if (is.null(x$conf.level)) 0.95 else x$conf.level
  data <- x$table
  if (null_count == "minimal" && nrow(data)) {
    minimum_by_comparator <- ave(
      data$NullCount, data$Comparator, FUN = min
    )
    data <- data[data$NullCount == minimum_by_comparator, , drop = FALSE]
  }
  data$Component <- factor(
    paste(data$Comparator, data$Endpoint, sep = "\n"),
    levels = unique(paste(data$Comparator, data$Endpoint, sep = "\n"))
  )
  data$Measure <- "Joint decision"
  data$Estimate <- data$Type1_Error
  data$Lower_CI <- data$Lower
  data$Upper_CI <- data$Upper
  data$Estimate_Label <- sprintf("%.1f%%", 100 * data$Estimate)
  data$Worst <- factor(
    seq_len(nrow(data)) == which.max(data$Type1_Error),
    levels = c(FALSE, TRUE)
  )
  boundary_levels <- unique(as.character(data$Boundary))
  boundary_colors <- stats::setNames(
    grDevices::hcl.colors(length(boundary_levels), palette = "Dark 3"),
    boundary_levels
  )
  if ("lower" %in% boundary_levels)
    boundary_colors[["lower"]] <- .type1_boundary_colors[["lower"]]
  if ("upper" %in% boundary_levels)
    boundary_colors[["upper"]] <- .type1_boundary_colors[["upper"]]
  dodge <- ggplot2::position_dodge(width = 0.55)
  global_upper <- if (is.null(x$global_upper)) NA_real_ else x$global_upper
  y_max <- max(c(x$alpha, data$Upper_CI, global_upper), na.rm = TRUE)
  y_max <- max(y_max * 1.15, y_max + 0.005)
  worst <- data[data$Worst == TRUE, , drop = FALSE]
  label_offset <- max(0.002, y_max * 0.015)
  # Supply the observed boundary groups for each component, including blank
  # labels for the non-highlighted groups.  Position dodge then uses the same
  # group count as the point layer, so the non-blank label sits exactly above
  # its highlighted point.
  label_data <- do.call(rbind, lapply(seq_len(nrow(worst)), function(i) {
    component <- as.character(worst$Component[[i]])
    component_boundaries <- unique(as.character(
      data$Boundary[data$Component == component]
    ))
    boundary <- as.character(worst$Boundary[[i]])
    data.frame(
      Component = factor(rep(component, length(component_boundaries)),
                         levels = levels(data$Component)),
      Boundary = factor(component_boundaries, levels = boundary_levels),
      label = ifelse(component_boundaries == boundary,
                     worst$Estimate_Label[[i]], ""),
      label_y = worst$Estimate[[i]] + label_offset,
      stringsAsFactors = FALSE
    )
  }))
  ggplot2::ggplot(
    data, ggplot2::aes(x = Component, y = Estimate, color = Boundary)
  ) +
    ggplot2::geom_hline(yintercept = x$alpha, linetype = "dashed",
                        color = "#343A40", linewidth = 0.6) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = Lower_CI, ymax = Upper_CI, linewidth = Worst,
                   group = Boundary),
      position = dodge, alpha = 0.75
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = Worst, group = Boundary),
      position = dodge, alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = data, ggplot2::aes(shape = Worst, group = Boundary),
      position = dodge, fill = NA, color = "black",
      stroke = 1.1, size = 3.6, show.legend = FALSE
    ) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = NA, `TRUE` = 21),
                                 guide = "none") +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = Component, y = label_y, label = label,
                   group = Boundary),
      position = dodge, vjust = 0, color = "#343A40",
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.9, `TRUE` = 2.0), guide = "none") +
    ggplot2::scale_size_manual(values = c(`FALSE` = 1.7, `TRUE` = 3.2), guide = "none") +
    ggplot2::scale_y_continuous(
      name = "Empirical Type I error", limits = c(0, y_max),
      labels = scales::label_percent(accuracy = 1)
    ) +
    ggplot2::scale_color_manual(values = boundary_colors, drop = FALSE) +
    ggplot2::labs(
      x = "Comparator and boundary endpoint(s)",
      color = "Boundary configuration",
      subtitle = paste0(
        "Larger point and thicker interval indicate the worst-case scenario; plot shows ",
        if (null_count == "minimal") "minimal null counts" else
          "all null counts"
      )
    ) +
    ggplot2::facet_wrap(~ Measure) +
    .diagnostic_theme() +
    ggplot2::theme(strip.text = ggplot2::element_blank())
}
