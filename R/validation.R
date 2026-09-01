#' @title Check Sample Size Limits
#'
#' @author
#' Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @description Validates that the upper and lower limits are numeric and that the upper limit is greater than the lower limit.
#'
#' @param lower Numeric. The initial lower limit for the search range.
#' @param upper Numeric. The initial upper limit for the search range.
#'
#' @return NULL. If the checks pass, the function returns nothing. If the checks fail, it stops execution with an error message.
#' @keywords internal
validate_sample_size_limits <- function(lower, upper) {

  # Check if both lower and upper are numeric
  if (!is.numeric(upper) || !is.numeric(lower)) {
    stop("The upper and lower limits for the sample size must be numeric.")
  }

  # Check if both lower and upper are integers
  if (lower != as.integer(lower) || upper != as.integer(upper)) {
    stop("The upper and lower limits for the sample size must be integers.")
  }

  # Check lower is greater than 0
  if (lower <= 0) {
    stop("The lower limit for the sample size must be greater than 0.")
  }

  # Check if upper is greater than or equal to lower
  if (upper < lower) {
    stop("The upper limit for the sample size must be greater than or equal to the lower limit.")
  }
}


#' @title Validate Positive Semi-Definite Matrices
#'
#' @author
#' Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @description Validates that all matrices in a list are symmetric and positive semi-definite.
#'
#' @param varcov_list List of matrices. Each matrix is checked to ensure it is symmetric and positive semi-definite.
#'
#' @return NULL. If all matrices pass, the function returns nothing. If any matrix fails, it stops with an error message.
#' @keywords internal
validate_positive_definite <- function(varcov_list) {
  # Function to check if a matrix is positive semi-definite
  is_positive <- function(x) {
    resp <- tryCatch({
      matrixcalc::is.positive.semi.definite(x)
    }, error = function(e) {
      FALSE
    })
    return(resp)
  }

  # Apply the check to each matrix in the list
  positive_definite_list <- unlist(lapply(varcov_list, is_positive))

  # Stop if any matrix is not positive semi-definite
  if (!all(positive_definite_list)) {
    stop("All matrices in 'varcov_list' must be symmetric and positive semi-definite.")
  }
}

.validate_common_controls <- function(power = NULL, alpha, nsim, ncores,
                                      dropout = NULL, rho = NULL) {
  if (!is.null(power) && (!is.numeric(power) || length(power) != 1L ||
                          !is.finite(power) || power < 0 || power >= 1))
    stop("'power' must be a finite number in [0, 1).")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a finite number strictly between 0 and 1.")
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != as.integer(nsim))
    stop("'nsim' must be a positive integer.")
  if (!is.numeric(ncores) || length(ncores) != 1L ||
      (!is.na(ncores) && (!is.finite(ncores) || ncores < 1 ||
                          ncores != as.integer(ncores))))
    stop("'ncores' must be a positive integer or NA.")
  if (!is.null(dropout) && !all(is.na(dropout))) {
    if (!is.numeric(dropout) || any(!is.finite(dropout)) ||
        any(dropout < 0 | dropout >= 1))
      stop("'dropout' must contain values in [0, 1).")
  }
  if (!is.null(rho) && (!is.numeric(rho) || any(!is.finite(rho)) ||
                        any(rho < -1 | rho > 1)))
    stop("'rho' must contain finite correlations in [-1, 1].")
  invisible(TRUE)
}


validate_sample_size_inputs <- function(mu_list, sigma_list, varcov_list,
                                        arm_names, list_comparator,
                                        list_y_comparator) {
  if (!is.list(mu_list) || length(mu_list) < 2 ||
      any(!vapply(mu_list, is.numeric, logical(1)))) {
    stop("'mu_list' must be a list of numeric vectors with at least two arms.")
  }
  if (length(arm_names) != length(mu_list) || anyDuplicated(arm_names) > 0) {
    stop("'arm_names' must contain one unique name for each arm.")
  }
  if (!is.list(list_comparator) || length(list_comparator) == 0L ||
      any(vapply(list_comparator, length, integer(1)) != 2L)) {
    stop("'list_comparator' must be a non-empty list of two-arm comparisons.")
  }
  if (any(!vapply(list_comparator, function(x) all(x %in% arm_names), logical(1)))) {
    stop("Every arm in 'list_comparator' must be present in 'arm_names'.")
  }
  if (!is.list(list_y_comparator) || length(list_y_comparator) != length(list_comparator)) {
    stop("'list_y_comparator' must have one endpoint vector per comparator.")
  }
  for (i in seq_along(list_comparator)) {
    common <- intersect(names(mu_list[[list_comparator[[i]][1]]]),
                        names(mu_list[[list_comparator[[i]][2]]]))
    if (is.null(names(mu_list[[list_comparator[[i]][1]]])) ||
        is.null(names(mu_list[[list_comparator[[i]][2]]]))) {
      next
    }
    if (length(list_y_comparator[[i]]) == 0L ||
        any(!list_y_comparator[[i]] %in% common)) {
      stop("Each comparator endpoint must be present in both comparator arms.")
    }
  }
  if (!(length(sigma_list) == 1L && is.na(sigma_list))) {
    if (!is.list(sigma_list) || length(sigma_list) != length(mu_list) ||
        any(vapply(sigma_list, length, integer(1)) != vapply(mu_list, length, integer(1)))) {
      stop("'sigma_list' must contain one standard-deviation vector per arm and endpoint.")
    }
  }
  if (!(length(varcov_list) == 1L && is.na(varcov_list)) &&
      (!is.list(varcov_list) || length(varcov_list) != length(mu_list))) {
    stop("'varcov_list' must contain one covariance matrix per arm.")
  }
  invisible(TRUE)
}

# Normalize user-facing distribution labels to R's distribution-function names.
# The longer labels retain compatibility with older examples while the stored
# value uses the canonical levels norm, lnorm, pois, and nbinom.
.normalize_distribution <- function(distribution) {
  if (!is.character(distribution) || length(distribution) != 1L ||
      is.na(distribution) || !nzchar(trimws(distribution))) {
    stop("'distribution' must be one of: norm, lnorm, pois, nbinom.")
  }
  value <- tolower(trimws(distribution))
  aliases <- c(
    norm = "norm", normal = "norm",
    lnorm = "lnorm", lognormal = "lnorm",
    `log normal` = "lnorm", `log-normal` = "lnorm",
    pois = "pois", poisson = "pois",
    nbinom = "nbinom", `negative binomial` = "nbinom",
    `negative-binomial` = "nbinom"
  )
  if (!value %in% names(aliases)) {
    stop("'distribution' must be one of: norm, lnorm, pois, nbinom.")
  }
  unname(aliases[[value]])
}
