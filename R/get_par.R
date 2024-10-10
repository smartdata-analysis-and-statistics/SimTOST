#' @title get_par
#' @description  Function to pass the parameters of endpoints and comparators in a list form required as input in the calopt function
#'
#' @param mu_list list of mu vectors, each element corresponds to a comparator with a mu vector of size number of endpoints.
#' @param varcov_list list of var-cov matrices, each element corresponds to a comparator with a varcov matrix of size number of endpoints X number of endpoints.
#' @param TAR_list list of allocation rates, each element corresponds to a comparator with a allocation number(integer or not).
#' @param type_y vector with the type of enpoints: primarly endpoint(1), otherwise (2).
#' @param arm_names vector with the names of the arms.
#' @param y_names vector with the names of the endpoints.
#'
#' @return named list with arm_names of mu, varcov and TAR to been passed to our function
#'
#' @export
get_par <- function( mu_list, varcov_list, TAR_list, type_y=NA, arm_names=NA, y_names=NA){

  len_list <- c(length(mu_list),length(varcov_list),length(TAR_list))
  len_mu <- lapply(mu_list,length)
  len_cvar <- lapply(varcov_list,ncol)
  len_y <- unlist(c(len_mu ,len_cvar)) # number of endpoints
  positive <- function(x){
    if(length(x)==1){
      resp <- TRUE
    }else{
      resp <- base::tryCatch({matrixcalc::is.positive.semi.definite(round(x,3))},
                             error=function(e) {FALSE})}
    return(resp)
  }
  lis_pdef <- unlist(lapply(varcov_list,positive))


  if(max(len_list) != min(len_list)){
    stop("mu,varcov and TAR should be defined for all the arms")
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
