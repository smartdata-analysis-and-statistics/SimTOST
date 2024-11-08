#' @title Parameter Configuration for Endpoints and Comparators
#' @description Constructs and returns a list of key parameters (mean vectors, variance-covariance matrices, and allocation rates) required for input into the `estSampleSize` function. This function ensures that the parameters for each endpoint and comparator are consistent, properly named, and formatted.
#'
#' @param mu_list A list of mean (\eqn{\mu}) vectors. Each element in the list represents a comparator, with the corresponding \eqn{\mu} vector having a length equal to the number of endpoints.
#' @param varcov_list A list of variance-covariance matrices. Each element corresponds to a comparator, with a matrix of size \eqn{(n \times n)}, where \eqn{n} is the number of endpoints.
#' @param TAR_list A list of treatment allocation rates (TARs) for each comparator. Each element contains a numeric value (can be fractional or integer) representing the allocation rate for the respective comparator.
#' @param type_y A numeric vector specifying the type of each endpoint. Use `1` for primary endpoints and `2` for secondary or other endpoints.
#' @param arm_names (Optional) A character vector containing names of the arms. If not provided, default names (e.g., T1, T2, ...) will be generated.
#' @param y_names (Optional) A character vector containing names of the endpoints. If not provided, default names (e.g., y1, y2, ...) will be generated.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\code{mu}}{A list of mean vectors, named according to \code{arm_names}.}
#'   \item{\code{varcov}}{A list of variance-covariance matrices, named according to \code{arm_names}.}
#'   \item{\code{tar}}{A list of treatment allocation rates (TARs), named according to \code{arm_names}.}
#'   \item{\code{type_y}}{A vector specifying the type of each endpoint.}
#'   \item{\code{weight_seq}}{A weight sequence calculated from \code{type_y}, used for endpoint weighting.}
#'   \item{\code{y_names}}{A vector of names for the endpoints, named as per \code{y_names}.}
#' }
#'
#' #' @details
#' This function ensures that all input parameters (\code{mu_list}, \code{varcov_list}, and \code{TAR_list}) are consistent across comparators and endpoints. It performs checks for positive semi-definiteness of variance-covariance matrices and automatically assigns default names for arms and endpoints if not provided.
#'
#' @examples
#' mu_list <- list(c(0.1, 0.2), c(0.15, 0.25))
#' varcov_list <- list(matrix(c(1, 0.5, 0.5, 1), ncol = 2), matrix(c(1, 0.3, 0.3, 1), ncol = 2))
#' TAR_list <- list(0.5, 0.5)
#' get_par(mu_list, varcov_list, TAR_list, type_y = c(1, 2), arm_names = c("Arm1", "Arm2"))
#'
#' @export
get_par <- function(mu_list, varcov_list, TAR_list, type_y=NA, arm_names=NA, y_names=NA){

  len_list <- c(length(mu_list),length(varcov_list),length(TAR_list))
  len_mu <- lapply(mu_list,length)
  len_cvar <- lapply(varcov_list,ncol)
  len_y <- unlist(c(len_mu ,len_cvar)) # number of endpoints

  positive <- function(x){
    if (length(x) == 1) {
      return(TRUE)
    }
    return(tryCatch({matrixcalc::is.positive.semi.definite(round(x,3))
    }, error = function(e) {
        FALSE
    }))
  }

  lis_pdef <- unlist(lapply(varcov_list,positive))


  if (max(len_list) != min(len_list)) {
    stop("Inconsistent lengths: 'mu_list', 'varcov_list', and 'TAR_list' must be defined for all arms.")
  }

  # Check consistency in the number of endpoints
  if (length(unique(len_y)) != 1) {
    stop("Inconsistent endpoint lengths: 'mu_list' and 'varcov_list' must have the same number of endpoints.")
  }


  if(max(len_y) != min(len_y)){
    stop("mu,varcov should be defined for all the endpoints")
  }

  if(!all(lis_pdef)){
    stop("all varcov should be symmetric and positive definite")
  }

  if(any(is.na(arm_names))){
    arm_names <- paste0("T",rep(1:length(mu_list)))
  }

  if (len_mu[[1]] == 1){
    mu_list <- lapply(mu_list,FUN = function(x){array(unlist(x))})
    varcov_list<-lapply(varcov_list,FUN = function(x){matrix(unlist(x))})}

  names(mu_list) <- arm_names
  names(varcov_list) <- arm_names
  names(TAR_list) <- arm_names

  if(any(is.na(type_y))){
    type_y <- rep(1,len_mu[[1]])}

  weight <- 1/table(type_y)
  weight_seq <-type_y

  for(x in unique(type_y)){
    weight_seq  [weight_seq == x] <- weight[x]
  }

  if(any(is.na(y_names))){
    y_names <- paste0("y",1:len_mu[[1]])}

  return(list(mu = mu_list, varcov = varcov_list, tar = TAR_list, type_y = type_y, weight_seq = weight_seq, y_names=y_names))
}
