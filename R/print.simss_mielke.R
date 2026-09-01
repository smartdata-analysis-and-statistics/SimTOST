#' Print a Mielke sample-size result
#' @param x A `simss_mielke` object.
#' @param ... Unused additional arguments.
#' @export
#' @method print simss_mielke
print.simss_mielke <- function(x, ...) {
  cat("Mielke sample-size calculation\n")
  cat("Required sample size:", x$SS, "\n")
  cat("Achieved power:", format(x$power.a, digits = 4), "\n")
  invisible(x)
}

#' @export
#' @method summary simss_mielke
summary.simss_mielke <- function(object, ...) {
  if (!inherits(object, "simss_mielke"))
    stop("Input must be of class 'simss_mielke'.")
  data.frame(
    design = object$parameters$design,
    n_per_sequence = object$n_per_sequence,
    n_total = object$n_total,
    target_power = object$target_power,
    achieved_power = object$power,
    m = object$parameters$m,
    k = object$parameters$k,
    rho = object$parameters$rho,
    nsim = object$parameters$nsim,
    row.names = NULL
  )
}

#' Coerce a Mielke result to a numeric sample size
#' @param x A `simss_mielke` object.
#' @param ... Unused additional arguments.
#' @export
as.numeric.simss_mielke <- function(x, ...) as.numeric(x$SS)

#' Preserve named-vector indexing for legacy Mielke examples
#' @param x A `simss_mielke` object.
#' @param i Index supplied to `[`. 
#' @param ... Unused additional arguments.
#' @export
`[.simss_mielke` <- function(x, i, ...) {
  if (missing(i)) return(x)
  if (is.character(i) && all(i %in% c("power.a", "SS"))) {
    return(stats::setNames(vapply(i, function(name) as.numeric(x[[name]]), numeric(1)), i))
  }
  NextMethod("[")
}
