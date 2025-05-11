#' Confidence Interval for Achieved Power from simss object
#'
#' @param object An object of class `"simss"` returned by a sampleSize function
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector with two elements:
#'   \describe{
#'     \item{Achieved Power}{Achieved power.}
#'     \item{Lower}{Lower bound of the confidence interval.}
#'     \item{Upper}{Upper bound of the confidence interval.}
#'   }
#' @export
#' @examples
#' # Assume `res` is a result from `sampleSize()`
#' # confint(res)
#'
confint.simss <- function(object, ...) {
  # Check if object is of class simss
  if (!inherits(object, "simss")) {
    stop("Object must be of class 'simss'")
  }
  # Look for CI directly (e.g., object$power.CI or object$power_ci)
  if (!is.null(object[["response"]][["power"]])) {
    ci <- c(object[["response"]][["power"]],object[["response"]][["power_LCI"]],object[["response"]][["power_UCI"]])
  }  else {
    stop("Confidence interval of achieved power not found in object.")
  }

  # Format as a named vector
  ci_out <- setNames(ci, c("Achieved Power", "Lower", "Upper"))

  cat(sprintf("Confidence Interval for Achieved Power (%.0f%%):\n", object[["param.d"]][["alpha"]] * 100))
  cat(sprintf(" %0.4f [%0.4f, %0.4f]\n", ci_out[1], ci_out[2], ci_out[3] ))

  invisible(ci_out)
}
