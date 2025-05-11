#' Summary for Simulation Results
#' @description Generates a summary of the simulation results, specifying the sample size for each comparator-endpoint
#' @param object An object of class `"simss"` returned by a sampleSize function
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector with the sample size of each arm and also the total (Total) sample size.
#' @export
#' @examples
#' # Assume `res` is a result from `sampleSize()`
#' # summary(res)
summary.simss <- function(object, ...) {
  # Equivalent margins
  margins <- data.table(names = names(unlist(object[["param.d"]][["list_lequi.tol"]])),
                        Lower= unlist(object[["param.d"]][["list_lequi.tol"]]),
                        Upper = unlist(object[["param.d"]][["list_uequi.tol"]]))
  margins[, c("Comparison", "Endpoint") := tstrsplit(names, "\\.")]

  cat("Sample Size Summary\n")
  cat("--------------------\n")

  # Sample size table
  ss <- as.data.frame(N_ss[["response"]][, !c("power","power_LCI", "power_UCI","n_iter", "n_drop"), with = FALSE])
  ss_names <- sub("^n_", "", colnames(ss))
  ss_names <- ifelse(ss_names == "total", "Total", ss_names)
  colnames(ss) <- ss_names

  if (!inherits(object, "simss")) {
    stop("Object must be of class 'simss'")
  }

  cat("Design:", object[["param.d"]][["dtype"]], "\n")
  cat("Comparison type:", object[["param.d"]][["ctype"]])
  cat("Equivalence Margins:\n")
  print(as.data.frame(margins[, c("Comparison", "Endpoint", "Lower", "Upper")]), row.names = FALSE)
  cat("Alpha:", object[["param.d"]][["alpha"]], "\n")
  cat("Target Power:", sprintf("%.4f",object[["param.d"]][["power"]]), "\n")
  cat("Achieved Power:", sprintf("%.4f",object[["response"]][["power"]]), "\n")
  cat("Estimated Sample Size:\n")
  print(ss, row.names = FALSE)
  if (!is.null(object$method)) {
    cat("Method:", object$method, "\n")
  }
  invisible(ss)
}

