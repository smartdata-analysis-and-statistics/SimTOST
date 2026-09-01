#' Summary for Simulation Results
#'
#' @description Generates a summary of the simulation results, including per-arm and total sample sizes. The printed summary also states the equivalence null and alternative hypotheses for the selected outcome and design.
#' @param object An object of class `"simss"` returned by a sampleSize function.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly, a data frame containing the estimated sample size for
#' each arm (or sequence), plus the total sample size. The printed report also
#' includes the design, estimand, equivalence margins, target power, achieved
#' power, and Monte Carlo interval. Count-outcome results use this same report
#' and retain their count-specific fields in the returned object.
#'
#' @author
#' Johanna Muñoz \email{johanna.munoz@fromdatatowisdom.com}
#'
#' @export summary.simss
#' @method summary simss
#' @examples
#' \dontrun{
#' res <- sampleSize(mu_list = list(T = c(y1 = 1), R = c(y1 = 1)),
#'                   sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
#'                   list_comparator = list(c("T", "R")),
#'                   list_lequi.tol = list(c(y1 = .8)),
#'                   list_uequi.tol = list(c(y1 = 1.25)),
#'                   ctype = "ROM", distribution = "lnorm", nsim = 10,
#'                   lower = 2, upper = 4)
#' summary(res)
#' }
summary.simss <- function(object, ...) {

  if (inherits(object, "countss")) {
    return(summary.countss(object, ...))
  }

  if (!inherits(object, "simss")) {
    stop("Input must be of class 'simss'")
  }

  # Build the margin table from the comparator definitions so that summary()
  # uses the same treatment-reference direction and labels as plot.simss().
  comparator_list <- object[["param"]][["list_comparator"]]
  lower_list <- object[["param.d"]][["list_lequi.tol"]]
  upper_list <- object[["param.d"]][["list_uequi.tol"]]
  endpoint_list <- object[["param"]][["list_y_comparator"]]
  margins <- data.table::rbindlist(lapply(seq_along(comparator_list), function(i) {
    lower <- as.numeric(lower_list[[i]])
    upper <- as.numeric(upper_list[[i]])
    endpoint_names <- names(lower_list[[i]])
    if (is.null(endpoint_names)) endpoint_names <- names(upper_list[[i]])
    if (is.null(endpoint_names)) endpoint_names <- endpoint_list[[i]]
    if (is.null(endpoint_names) || length(endpoint_names) != length(lower)) {
      endpoint_names <- paste0("y", seq_along(lower))
    }
    data.frame(
      Comparison = paste(comparator_list[[i]], collapse = "_vs_"),
      Endpoint = endpoint_names,
      Lower = lower,
      Upper = upper,
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)

  # Header
  cat("Sample Size Summary\n")
  cat(strrep("-", 22), "\n")

  # Sample size results
  ss <- as.data.frame(object[["response"]][, !c("power","power_LCI", "power_UCI","n_iter", "n_drop"), with = FALSE])
  ss_names <- sub("^n_", "", colnames(ss))
  colnames(ss) <- ifelse(colnames(ss) == "total", "Total", colnames(ss))

  # Display design and summary
  cat("Design type        :", object[["param.d"]][["dtype"]], "\n")
  cat("Distribution       :", object[["param.d"]][["distribution"]], "\n")
  cat("Comparison type    :", object[["param.d"]][["ctype"]], "\n")
  if (identical(object[["param.d"]][["distribution"]], "norm")) {
    estimand <- if (object[["param.d"]][["ctype"]] == "DOM") "mean difference" else "mean ratio"
    cat("Estimand           :", estimand, "\n")
    cat("Hypotheses         : H0: estimand <= L or >= U; H1: L < estimand < U\n")
  } else if (identical(object[["param.d"]][["distribution"]], "lnorm")) {
    cat("Estimand           : arithmetic mean ratio\n")
    cat("Hypotheses         : H0: ratio <= L or >= U; H1: L < ratio < U (tested on the log scale)\n")
  }
  cat("Alpha              :", object[["param.d"]][["alpha"]], "\n")
  cat("Target power       :", sprintf("%.4f", object[["param.d"]][["power"]]), "\n")
  cat("Achieved power     :", sprintf("%.4f", object[["response"]][["power"]]), "\n")
  if (!is.null(object$method)) {
    cat("Method             :", object$method, "\n")
  }

  cat("\nEquivalence Margins:\n")
  print(as.data.frame(margins[, c("Comparison", "Endpoint", "Lower", "Upper"), with = FALSE]), row.names = FALSE)

  cat("\nEstimated Sample Size:\n")
  print(ss, row.names = FALSE)

  invisible(ss)
}
