#' Confidence Interval for Achieved Power from simss object
#'
#' @param object An object of class `"simss"` returned by a sampleSize function
#' @param parm Unused; included for compatibility with [stats::confint()].
#' @param level Confidence level for the Monte Carlo interval.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector with two elements:
#'   \describe{
#'     \item{Achieved Power}{Achieved power.}
#'     \item{Lower}{Lower bound of the confidence interval.}
#'     \item{Upper}{Upper bound of the confidence interval.}
#'   }
#' @export
#' @method confint simss
#' @examples
#' \dontrun{
#' # confint(res), where res is returned by sampleSize()
#' }
#'
confint.simss <- function(object, parm = NULL, level = 0.95, ...) {
  if (!is.null(parm) && missing(level) && is.numeric(parm) &&
      length(parm) == 1L && is.finite(parm) && parm > 0 && parm < 1) {
    level <- parm
    parm <- NULL
  }
  if (inherits(object, "countss")) {
    return(confint.countss(object, level = level, ...))
  }

  # Check if object is of class simss
  if (!inherits(object, "simss")) {
    stop("Object must be of class 'simss'")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number between 0 and 1.")
  if (!isTRUE(all.equal(level, 0.95)))
    stop("Only the stored 95% interval is available for this result.")

  # Look for CI directly (e.g., object$power.CI or object$power_ci)
  if (!is.null(object[["response"]][["power"]])) {
    ci <- c(object[["response"]][["power"]],object[["response"]][["power_LCI"]],object[["response"]][["power_UCI"]])
  }  else {
    stop("Confidence interval of achieved power not found in object.")
  }

  # Format as a named vector
  ci_out <- setNames(ci, c("Achieved Power", "Lower", "Upper"))

  cat("Confidence Interval for Achieved Power (95%):\n")
  cat(sprintf(" %0.4f [%0.4f, %0.4f]\n", ci_out[1], ci_out[2], ci_out[3]))

  invisible(ci_out)
}
