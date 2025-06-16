#' Summary for Simulation Results
#'
#' @description Generates a summary of the simulation results, including per-arm and total sample sizes.
#' @param object An object of class `"simss"` returned by a sampleSize function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector with the sample size per arm and the total (Total) sample size.
#'
#' @author
#' Johanna Muñoz \email{johanna.munoz@fromdatatowisdom.com}
#'
#' @export summary.simss
#' @examples
#' # Assume `res` is a result from `sampleSize()`
#' # summary(res)
summary.simss <- function(object, ...) {

  if (!inherits(object, "simss")) {
    stop("Input must be of class 'simss'")
  }

  # Equivalent margins
  margins <- data.table::data.table(
    names = names(unlist(object[["param.d"]][["list_lequi.tol"]])),
    Lower = unlist(object[["param.d"]][["list_lequi.tol"]]),
    Upper = unlist(object[["param.d"]][["list_uequi.tol"]])
    )
  margins[, c("Comparison", "Endpoint") := data.table::tstrsplit(names, "\\.")]

  # Header
  cat("Sample Size Summary\n")
  cat(strrep("-", 22), "\n")

  # Sample size results
  ss <- as.data.frame(object[["response"]][, !c("power","power_LCI", "power_UCI","n_iter", "n_drop"), with = FALSE])
  ss_names <- sub("^n_", "", colnames(ss))
  colnames(ss) <- ifelse(colnames(ss) == "total", "Total", colnames(ss))

  # Display design and summary
  cat("Design type        :", object[["param.d"]][["dtype"]], "\n")
  cat("Comparison type    :", object[["param.d"]][["ctype"]], "\n")
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

