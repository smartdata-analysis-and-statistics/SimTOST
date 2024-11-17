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
validate_positive_definite <- function(varcov_list) {
  # Function to check if a matrix is positive semi-definite
  is_positive <- function(x) {
    resp <- tryCatch({
      matrixcalc::is.positive.semi.definite(round(x, 3))
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

#' Validate Treatment Allocation Rate (TAR)
#'
#' This function validates whether the length of the treatment allocation rate (`TAR`) matches
#' the number of treatment arms (`arm_names`). If the lengths do not match, it throws an error.
#'
#' @param TAR Numeric vector. Treatment allocation rates for each arm.
#' @param n_arms Numeric. The number of treatment arms that `TAR` must correspond to.
#'
#' @author
#' Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @return NULL (used for validation only).
validate_tar <- function(TAR = NULL, n_arms) {
  if (missing(TAR)) {
    return(NULL)
  }
  if (is.null(TAR)) {
    return(NULL)
  }

  # Validate TAR length
  if (length(TAR) != n_arms) {
    stop("Validation Error: The length of TAR must match the number of arms specified by arm_names.")
  }

  if (any(TAR <= 0, na.rm = TRUE)) {
    stop("Validation Error: TAR must contain only positive values. Negative or zero values are not allowed.")
  }


}

