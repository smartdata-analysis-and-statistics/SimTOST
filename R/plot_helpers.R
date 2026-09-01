# Internal helpers shared by power and sample-size plots.

#' Select the smallest plotted sample size that reaches the target power.
#'
#' The sample-size search may evaluate candidates in a non-monotone order.
#' Plot selection must therefore be based on the plotted global total-power
#' rows, not on the order in which the optimizer evaluated candidates.
#'
#' @keywords internal
.select_power_plot_n <- function(data, target_power, n_col,
                                 fallback = NA_real_) {
  if (!is.data.frame(data) || !nrow(data) ||
      !all(c(n_col, "power", "Endpoint") %in% names(data))) {
    return(as.numeric(fallback)[1L])
  }

  total <- data[as.character(data$Endpoint) == "Total", , drop = FALSE]
  global_labels <- c("All comparators", "All comparisons")
  if ("Comparator" %in% names(total) &&
      any(total$Comparator %in% global_labels, na.rm = TRUE)) {
    total <- total[total$Comparator %in% global_labels, , drop = FALSE]
  }
  total_n <- as.numeric(total[[n_col]])
  reached <- total[is.finite(total_n) & is.finite(total$power) &
                     total$power >= target_power, , drop = FALSE]
  if (nrow(reached))
    return(min(as.numeric(reached[[n_col]]), na.rm = TRUE))

  as.numeric(fallback)[1L]
}

.stored_plot_target_power <- function(x, fallback = 0.80) {
  candidates <- list(
    x$target_power,
    if (!is.null(x$param.d)) x$param.d$power else NULL,
    if (!is.null(x$parameters)) x$parameters$power else NULL
  )
  for (value in candidates) {
    value <- as.numeric(value)[1L]
    if (is.finite(value) && value > 0 && value < 1)
      return(value)
  }
  fallback
}
