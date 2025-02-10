#' @title Generate Simulated Endpoint Data for Parallel Group Design
#'
#' @description Generate simulated endpoint data for a parallel design, with options for normal and lognormal distributions.
#'
#' @author
#' Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @param n Integer. The sample size for the generated data.
#' @param mu.arithmetic Numeric vector. The arithmetic mean of the endpoints on the original scale.
#' @param mu.geometric Numeric vector. The geometric mean of the endpoints on the original scale. Only used if `dist = "lognormal"`.
#' @param Sigma Matrix. Variance-covariance matrix of the raw data on the original scale. If `dist = "lognormal"`, this matrix is transformed to the log scale.
#' @param CV Numeric vector. Coefficient of variation (CV) of the raw data. Only used when `dist = "lognormal"`, where it is transformed to the log scale.
#' @param seed Integer. Seed for random number generation, ensuring reproducibility.
#' @param dist Character. Assumed distribution of the endpoints: either `"normal"` or `"lognormal"`.
#'
#' @return A matrix of simulated endpoint values for a parallel design, with dimensions `n` by the number of variables in `mu.arithmetic` or `mu.geometric`.
#'
#' @export
#'
simParallelEndpoints <- function(n,
                                 mu.arithmetic,
                                 mu.geometric = NULL,
                                 Sigma,
                                 CV = NULL, # Vector of size (mu.arithmetic)
                                 seed,
                                 dist = "normal") {

  if (dist == "normal") {
    dmu <- mu.arithmetic
    dsigma <- Sigma
  } else if (dist == "lognormal" & !is.null(mu.arithmetic)) {
    dsigma <- log(Sigma/(mu.arithmetic %*% t(mu.arithmetic)) + 1)
    dmu   <- log(mu.arithmetic) - 1/2*diag(dsigma)
  } else if (dist == "lognormal" & !is.null(mu.geometric)) {
    dsigma <- log(1 + CV**2)
    dmu   <- log(mu.geometric)
  } else {
    stop("Invalid distribution")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  return(MASS::mvrnorm(n = n, mu = dmu, Sigma = dsigma))
}

#' @title Calculate the power across all comparators
#' @description  Internal function to calculate the power across all comparators
#'
#' @param n sample size
#' @param nsim number of simulated studies
#' @param param list of parameters (mean,sd,tar)
#' @param seed main seed
#' @param ncores number of cores
#' @param param.d design parameters
#'
#' @return power calculated from a global list of comparators
#' @keywords internal
#'
power_cal <- function(n,nsim,param,param.d,seed,ncores){

  if (param.d$dtype == "parallel") {
    TAR_used <- unlist(param$TAR_list)[unique(unlist(param$list_comparator))]
    size <- ceiling(n*TAR_used)
    size[size < 2] <- 2
    size_ndrop <- ceiling((1 - param.d$dropout[names(size)])*size)
    size_ndrop[size_ndrop < 2] <- 2
    n_drop <- sum(size)-sum(size_ndrop)

  } else if (param.d$dtype == "2x2") {
    # expected
    size <- NULL
    size_drop <- NULL
    for (j in seq(length(param$list_comparator))) {
      comp <- param$list_comparator[[j]]
      ns0i <- ceiling(n/2) # n/2 per sequence
      ns1i <- n - ns0i # n/2 per sequence
      ns0i <- ifelse(ns0i < 2, 2, ns0i)
      ns1i <- ifelse(ns1i < 2, 2, ns1i)
      # no drop out
      ns0 <- ceiling((1 - param.d$dropout[1])*ns0i)
      ns0 <- ifelse(ns0 < 2, 2, ns0)
      ns1 <- ceiling((1 - param.d$dropout[2])*ns1i)
      ns1 <- ifelse(ns1 < 2 ,2, ns1)

      # Expected per sequence
      sizej <- c(ns0i, ns1i)
      # Drop out per sequence
      size_dropj <- c(ns0i - ns0, ns1i - ns1)
      names(sizej) <- names(size_dropj) <- paste0(c("seq0_","seq1_"), paste0(comp, collapse = "vs"))
      size <- c(size,sizej)
      size_drop <- c(size_drop,size_dropj)
    }
    n_drop <- sum(size_drop)
  } else {
    stop("Invalid design type")
  }

  arm_names <- param$arm_names
  if (is.na(seed)) {
    seed <- sample(1:2^15,1)
  }
  set.seed(seed)

  # Draw a unique random seed for each arm in each simulation
  if (param.d$dtype == "parallel") {
    arm_seed <- matrix(sample(x = seq((length(arm_names)*nsim*100)),
                              size = length(arm_names)*nsim,
                              replace = FALSE),
                       ncol = length(arm_names))
    colnames(arm_seed) <- arm_names
  } else if (param.d$dtype == "2x2") {
    # Not same seed because the indiviudals change on each 2x2 study
    arm_seed <- matrix(sample(x = seq((length(param$list_comparator)*nsim*100)),
                              size = length(arm_names)*nsim,
                              replace = FALSE),
                       ncol = length(param$list_comparator))
  } else {
    stop("Invalid design type")
  }

  test_listcomp <- do.call("rbind",lapply(1:length(param$list_comparator),test_studies,nsim=nsim,n=n,param=param,param.d=param.d,arm_seed=arm_seed,ncores=ncores))
  tbiocom_listcomp <- test_listcomp[grep("^totaly",rownames(test_listcomp)),]

  if (is.null(nrow(tbiocom_listcomp))){ # only one comparator
    t_true <- sum(tbiocom_listcomp)
  }else{
    t_true <- sum(apply(tbiocom_listcomp, 2, prod))

  }


  # Filter only the TAR of arms used

  size <- c(size,total = sum(size))
  names(size) <- paste0("n_",names(size))
  output.test <- as.data.frame(t(rbind(test_listcomp)))
  output.test$n_iter <- n
  output.test$t_true <- t_true
  output.test$n_drop <- n_drop
  for(i in 1:length(size)){
    output.test[,names(size[i])] <- size[i]
  }
  return(list(power = t_true/nsim,
              output.test = output.test))
}

#' @title test_studies
#' @description  Internal function to estimate the bioequivalence test for nsim simulated studies given a sample size n
#' @param nsim number of simulated studies
#' @param n sample size
#' @param comp index comparator
#' @param param list of parameters (mean,sd,tar)
#' @param arm_seed seed for each endpoint to get consistent in simulations across all comparators
#' @param ncores number of cores used for the calculation
#' @param param.d design parameters
#'
#' @return a logical matrix of size  (nsim) X (number of endpoints + 1) function only replicates test_bioq nsim times.
#'
#' @keywords internal

test_studies <- function(nsim, n, comp, param, param.d, arm_seed, ncores){
  if(is.na(ncores)){
    ncores <- parallel::detectCores() -1}
  treat1 <- param$list_comparator[[comp]][[1]]
  treat2 <- param$list_comparator[[comp]][[2]]
  endp <- param$list_y_comparator[[comp]]

  m <- length(endp) # number of endpoints

  muT <- param$mu[[treat1]][,endp] # treatment given by user
  muR <- param$mu[[treat2]][,endp] # reference given by user

  if (param.d$dtype == "parallel"){
    SigmaT<- param$varcov[[treat1]][endp,endp]
    SigmaR<- param$varcov[[treat2]][endp,endp]
  }else{
    SigmaW <- param$varcov[[treat2]][endp,endp] # Within subjects variance in previous experiment
    #To be added on main list of parameters
    sigmaB <- param$sigmaB # Between subjects variance
    sigmaB <- ifelse(is.na(sigmaB), if (length(SigmaW)==1) 2* sqrt(SigmaW) else 2 * sqrt(max(diag(SigmaW))), sigmaB) # Assumes to be at least the double of the max within variance
  }

  # Set equivalence tolerance
  lequi.tol <- param.d$list_lequi.tol[[comp]][endp]
  uequi.tol <- param.d$list_uequi.tol[[comp]][endp]
  dropout <- param.d$dropout
  alphau <- param.d$alpha # alpha unique @Thomas, here we can divide by the number of comparators in case you want a bonberroni across comparators.Think about this!!
  adjust <- param.d$adjust
  k <- param.d$k[[comp]]
  # alpha vector
  if (adjust=="no") {alpha <- rep(alphau,m)}
  if (adjust=="bon") {alpha <- rep(alphau/(m),m)}
  if (adjust=="sid") {alpha <- rep(1-(1-alphau)^{1/m},m)}
  if (adjust=="k") { alpha <- rep(k*alphau/(m),m)}
  if (adjust=="seq"){
    alpha <- alphau*param$weight_seq[endp]
    }

  if(param.d$ctype=="ROM"&param.d$lognorm == TRUE){
    if (param.d$dtype == "parallel"){
      SigmaT <-  as.matrix(log(SigmaT/(muT%*%t(muT))+1))
      SigmaR <-  as.matrix(log(SigmaR/(muR%*%t(muR))+1))
      muT <- log(muT)-1/2*diag(SigmaT)
      muR <- log(muR)-1/2*diag(SigmaR)
    }else{
      SigmaW <-  as.matrix(log(SigmaW/(muR%*%t(muR))+1))
      sigmaB <-  sigmaB #log(sigmaB/(muR%*%t(muR))+1)
      muR <- log(muR)-1/2*diag(SigmaW)
      muT <- log(muT)-1/2*diag(SigmaW)
    }
    lequi.tol <- log(lequi.tol)
    uequi.tol <- log(uequi.tol)
  }

  # Get typey the positions of the primarly in C++
  if (any(param$type_y[endp] == 1) ) {
    typey <- which(param$type_y[endp] == 1) - 1
  } else { # in case no primary endpoint is specified
    typey = -1
  }


  result <- mcsapply(1:nsim, function(i){
    arm_seedx <- arm_seed[i,]
    if(param.d$dtype == "parallel" ) {
      if(param.d$ctype == "DOM"|(param.d$ctype=="ROM"&param.d$lognorm == TRUE) ) {
        outtest <- as.vector(test_par_dom(n = n, muT=muT, muR = muR,
                                          SigmaT = as.matrix(SigmaT),
                                          SigmaR = as.matrix(SigmaR),
                                          lequi_tol = lequi.tol ,uequi_tol = uequi.tol,
                                          alpha = alpha, k = k,
                                          dropout = as.numeric(c(dropout[treat1], dropout[treat2])),
                                          typey = typey,
                                          adseq = param.d$adjust == "seq",
                                          arm_seedT = arm_seedx[treat1],
                                          arm_seedR = arm_seedx[treat2],
                                          TART = param$TAR_list[[treat1]],
                                          TARR = param$TAR_list[[treat2]],
                                          vareq =param.d$vareq))


      }else{ #ROM & normal distribution
        outtest <- as.vector(test_par_rom(n=n, muT=muT, muR =muR,
                                          SigmaT = as.matrix(SigmaT),
                                          SigmaR = as.matrix(SigmaR),
                                          lequi_tol = lequi.tol,
                                          uequi_tol = uequi.tol,
                                          alpha = alpha, k = k,
                                          dropout = as.numeric(c(dropout[treat1],dropout[treat2])),
                                          typey = typey,
                                          adseq = param.d$adjust=="seq",
                                          arm_seedT = arm_seedx[treat1],
                                          arm_seedR = arm_seedx[treat2],
                                          TART = param$TAR_list[[treat1]],
                                          TARR = param$TAR_list[[treat2]],
                                          vareq =param.d$vareq))
      }

      names(outtest) <- paste0(c("totaly", endp,
                                 paste0("mu_",endp,"_",treat1),
                                 paste0("mu_",endp,"_",treat2),
                                 paste0("sd_",endp,"_",treat1),
                                 paste0("sd_",endp,"_",treat1)),"Comp:",treat1," vs ",treat2)

    }else{ # 2X2
      if(param.d$ctype == "DOM"|(param.d$ctype=="ROM"&param.d$lognorm == TRUE)){
        outtest <- as.vector(test_2x2_dom(n=n, muT=muT, muR=muR,
                                          SigmaW = as.matrix(SigmaW),
                                          lequi_tol = lequi.tol,
                                          uequi_tol = uequi.tol,
                                          sigmaB= sigmaB,
                                          dropout=dropout,
                                          typey = typey,
                                          adseq=param.d$adjust=="seq",
                                          Eper = param$Eper, Eco=param$Eco,
                                          arm_seed = arm_seedx[comp],
                                          alpha = alpha,
                                          k=k))
        names(outtest) <- paste0(c("totaly", endp,
                                   paste0("mu_",endp,"_",treat1),
                                   paste0("mu_",endp,"_",treat2),
                                   paste0("sdw_",endp,"_",treat1),
                                   paste0("sdb_",endp,"_",treat1)),"Comp:",treat1," vs ",treat2)

      }else{ #ROM & normal distribution
        outtest <- as.vector(test_2x2_rom(n=n, muT=muT, muR=muR,
                                          SigmaW=as.matrix(SigmaW),
                                          lequi_tol = lequi.tol,
                                          uequi_tol = uequi.tol,
                                          sigmaB= sigmaB,
                                          dropout=dropout,
                                          typey = typey,
                                          adseq=param.d$adjust=="seq",
                                          Eper=param$Eper, Eco=param$Eco, arm_seed=arm_seedx[comp],
                                          alpha=alpha,
                                          k=k))

      }
      names(outtest) <- paste0(c("totaly", endp,
                                 paste0("mu_",endp,"_",treat1),
                                 paste0("mu_",endp,"_",treat2),
                                 paste0("sdw_",endp,"_",treat1),
                                 paste0("sdb_",endp,"_",treat1)),"Comp:",treat1," vs ",treat2)

    }

    outtest
  }, mc.cores = ncores)
  return(result)
}



#' @title uniroot.integer.mod
#' @description  Optimizer Uniroot integer modified from the ssanv function package https://github.com/cran/ssanv/blob/master/R/uniroot.integer.R
#'
#' @param f  function for which a root is needed
#' @param power target power value
#' @param lower minimum allowable root
#' @param upper maximum allowable root
#' @param step.power initial step size is 2^step.power
#' @param step.up if TRUE steps up from 'lower', if FALSE steps down from 'upper'
#' @param pos.side if TRUE finds integer, i, closest to the root such that f(i) >0
#' @param maxiter maximum number of iterations
#' @param ... additional arguments to 'f'.
#'
#' @return a list with the following elements:
#' root= the integer on the correct side of the root
#' f.root	= value of f at root  (output= main output with final power, output.test= output table with test at endpoint level)
#' iter	= number of times f was evaluated
#' table.iter = data.frame with the estimated N and power at each iteration
#' table.test = data.frame with the test at endpoint level for a given N and each simulation draw
#'
#' @keywords internal

uniroot.integer.mod <-function (f, power, lower = lower, upper = upper, step.power=step.power, step.up=step.up, pos.side=pos.side, maxiter = maxiter,...) {
  iter <- 0
  table.test<-data.frame()
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper)
    stop("lower < upper  is not fulfilled")
  if (lower==-Inf && step.up==TRUE) stop("lower cannot be -Inf when step.up=TRUE")
  if (upper==Inf && step.up==FALSE) stop("upper cannot be Inf when step.up=FALSE")
  if (step.up){
    f.old<-f(lower,...)
    iter<-iter+1
    sign<-1
    xold<-lower }
  else{
    f.old<-f(upper,...)
    iter<-iter+1
    sign<- -1
    xold<-upper
  }

  ever.switched<-FALSE
  tried.extreme<-FALSE
  while (step.power>-1){

    if ((power-f.old$power)==0) break()
    if (iter>=maxiter) stop("reached maxiter without a solution")
    xnew<- xold + sign*(2^step.power)
    if ((step.up & xnew< upper) || (!step.up & xnew> lower) ){
      f.new<-f(xnew,...)
      iter<-iter+1
      if(!xold%in%c(table.test$n_iter)){
        table.test<-rbind(table.test,f.old$output.test)}
    }
    else{

      xnew<- xold
      f.new<-f.old
      step.power<-step.power-1
      if (tried.extreme==FALSE){
        if (step.up){ f.extreme <- f(upper,...); iter<-iter+1; x.extreme<-upper }
        else{ f.extreme <- f(lower,...); iter<-iter+1; x.extreme<-lower }
        tried.extreme <- TRUE
        xswitch <- x.extreme
        f.switch <- f.extreme
        if ((power-f.extreme$power)==0){
          xold<-x.extreme
          f.old<-f.extreme
          break()
        }

        if (((power-f.old$power)*(power-f.extreme$power))>=0){
          warning("f() at extremes not of opposite sign, try to set up upper level to a higher number")
          return(list(iter=iter,f.root=f(upper,...),root=upper,table.test=table.test))
        }
      }
    }

    if ( ((power-f.old$power)*(power-f.new$power))<0){
      sign<- sign*(-1)
      ever.switched<-TRUE
      xswitch<-xold
      f.switch<-f.old
    }
    if (ever.switched){
      step.power<-step.power-1
    }

    xold<- xnew
    f.old<-f.new
    if(step.power<0){
      if(!xold%in%c(table.test$n_iter)){
        table.test<-rbind(table.test,f.old$output.test)}
    }
  }

  if ((power-f.old$power)==0){
    root<-xold
    f.root<-f.old
  } else if ((power-f.new$power)==0){
    root<-xnew
    f.root<-f.new

  } else if ((power-f.switch$power)==0){
    root <- xswitch
    f.root <- f.switch
  } else if (pos.side){
    root <- if((power-f.new$power)>0) xnew else xswitch
    f.root<-if((power-f.new$power)>0) f.new else f.switch
  } else {
    root<-if((power-f.new$power)<0) xnew else xswitch
    f.root<-if((power-f.new$power)<0) f.new else f.switch
  }

  if(!root%in%c(table.test$n_iter)){
    table.test<-rbind(table.test,f.old$output.test)}

  power <- c(root,f.root$power)
  names(power) <- c("n_iter","power")

  return(list(power=power,
              table.test=table.test))

}

#' @title mcsapply
#' @description An mc-version of the sapply function. https://stackoverflow.com/questions/31050556/parallel-version-of-sapply
#'
#' @param X  vector of iterations
#' @param FUN function
#' @param ...  additional parameters to pass
#' @param simplify  simplify array
#' @param USE.NAMES  use names in array
#'
#' @return vector output
#'
#'@keywords internal
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- base::match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && base::is.character(X) && base::is.null(names(answer)))
    base::names(answer) <- X
  if (!base::isFALSE(simplify) && base::length(answer))
    base::simplify2array(answer, higher = (simplify == "array"))
  else answer
}

#' Helper function for conditional messages
#'
#' This function displays a message if the `verbose` parameter is set to `TRUE`.
#' It is useful for providing optional feedback to users during function execution.
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @param message A character string containing the message to display.
#' @param verbose Logical, if `TRUE`, the message is displayed; if `FALSE`, the message is suppressed.
#'
#' @return NULL (invisible). This function is used for side effects (displaying messages).
info_msg <- function(message, verbose) {
  if (verbose) message(message)
}

