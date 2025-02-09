#' @title Sample Size Calculation for Bioequivalence and Multi-Endpoint Studies
#' @description This function calculates the required sample size to achieve a target power in studies with multiple endpoints and treatment arms. The function leverages modified root-finding algorithms to estimate sample size while considering correlation structures, variance assumptions, and equivalence bounds across endpoints. It is especially useful for bioequivalence trials or multi-arm trials with complex endpoint structures.
#'
#' @param mu_list Named list of arithmetic means per treatment arm. Each element contains a vector (i.e., one per treatment arm) with the expected outcomes for all endpoints of interest.
#' @param varcov_list list of var-cov matrices, each element corresponds to a comparator with a varcov matrix of size number of endpoints X number of endpoints.
#' @param sigma_list  list of sigma vectors, each element corresponds to a comparator with a sigma vector of size number of endpoints.
#' @param cor_mat matrix specifying the correlation matrix between endpoints, used along with sigma_list  to calculate the varcov list in case it is not provided.
#' @param sigmaB number between subject variance only for 2x2 design.
#' @param Eper Optional. Vector of length 2 specifying the period effect in a `dtype = "2x2"` design, applied to c(Period 0, Period 1). Defaults to `c(0, 0)` if not provided. Ignored for `dtype = "parallel"`.
#' @param Eco Optional. Vector of length 2 specifying the carry-over effect for each arm in a `dtype = "2x2"` design, applied to c(Reference, Treatment). Defaults to `c(0, 0)` if not provided. Ignored for `dtype = "parallel"`.
#' @param rho Correlation parameter applied uniformly across all endpoint pairs, used with sigma_list to calculate varcov if cor_mat or varcov_list are not provided.
#' @param TAR Numeric vector. Treatment allocation rates for each arm, where the order of values corresponds to the order of `arm_names`. The length of `TAR` must match the number of arms. If not provided, a default equal allocation rate is assigned across all arms.
#' @param arm_names Optional vector with the treatment names. If not supplied, it will be derived from mu_list.
#' @param ynames_list Optional list of vectors with Endpoint names on each arm. When not all endpoint names are provided for each arm, arbitrary names (assigned by vector order) are used.
#' @param type_y vector with the type of endpoints: primary endpoint(1), otherwise (2).
#' @param list_comparator list of comparators, i.e each comparator is a vector of size 1 X 2 where are specified the name of treatments
#' @param list_y_comparator list of endpoints to be considered in each comparator. Each element of the list is a vector containing the names of the endpoints to compare. When it is not provided, all endpoints present in both compared arms are used.
#' @param power target power (default = 0.8)
#' @param alpha alpha level (default = 0.05)
#' @param lequi.tol lower equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints) of endpoint repeated on all endpoints and comparators
#' @param uequi.tol upper equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints) of endpoint repeated on all endpoints and comparators
#' @param list_lequi.tol list of lower equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints) of endpoint in comparator
#' @param list_uequi.tol list of upper equivalence bounds (e.g., -0.5) expressed in raw scale units (e.g., scalepoints) of endpoint in comparator
#' @param vareq Logical indicating whether variances are assumed equal across arms (default = FALSE).
#' @param dtype Character. Design type for the trial: `"parallel"` (default) for parallel group design or `"2x2"` for crossover design (applicable only for trials with 2 arms).
#' @param lognorm Is data log-normally distributed? (TRUE, FALSE)
#' @param k Vector with the number of endpoints that must be successful (integer) for global bioequivalence for each comparator. If no k vector is provided, it will be set to the total number of endpoints on each comparator.
#' @param adjust Character. Method for alpha adjustment: `"k"` (K-fold), `"bon"` (Bonferroni), `"sid"` (Sidak), `"no"` (no adjustment, default), or `"seq"` (sequential adjustment).
#' @param ctype Character. Specifies the type of hypothesis test for comparison: `"DOM"` for Difference of Means or `"ROM"` for Ratio of Means.
#' @param dropout vector with proportion of total population with dropout per arm
#' @param nsim number of simulated studies (default=5000)
#' @param seed main seed
#' @param ncores Integer. Number of processing cores to use for parallel computation. Defaults to one less than the total number of detected cores.
#' @param optimization_method Character. Method for determining the required sample size: "fast" (using modified root-finding algorithms) or "step-by-step". Defaults to "fast".
#' @param lower Integer. Initial value of `N` for the search range. Defaults to 2.
#' @param upper Integer. Maximum value of `N` for the search range. Defaults to 500.
#' @param step.power Numeric. The initial step size for the sample size search, defined as `2^step.power`. Relevant when `optimization_method` is `"fast"`.
#' @param step.up Logical. If `TRUE` (default), the sample size search increments upward from the `lower` limit; if `FALSE`, it decrements downward from the `upper` limit. Used only when `optimization_method` is `"fast"`.
#' @param pos.side Logical. If `TRUE`, finds the smallest integer, `i`, closest to the root such that `f(i) > 0`. Used only when `optimization_method` is `"fast"`.
#' @param maxiter Integer. Maximum number of iterations allowed for finding the sample size. Defaults to 1000. Used only when `optimization_method` is `"fast"`.
#' @param verbose Logical. If `TRUE`, the function displays progress and informational messages during execution. Defaults to `FALSE`.
#' @return An object simss that contains the following elements :
#' \describe{
#'  \item{"response"}{ An array summarizing the results of the simulation, including the estimated sample sizes for each arm, the approximated achieved power, and the corresponding confidence interval.}
#'  \item{"table.iter"}{A data frame detailing the estimated sample size for each arm and the calculated power at each iteration during the sample size searching process.}
#'  \item{"table.test"}{A data frame containing the test results for all simulated trials, including relevant metrics for each tested sample size}
#'  \item{"param.u"}{The original set of parameters provided by the user for the simulation.}
#'  \item{"param"}{ The final set of parameters used for the sample size calculation. These are based on \code{param.u} but adjusted to address any inconsistencies or missing information.}
#'  \item{"param.d"}{ The design parameters used in the simulation, including details relevant to the trial design.}
#'}
#'
#' @references
#' Mielke, J., Jones, B., Jilma, B., & König, F. (2018). Sample size for multiple hypothesis testing in biosimilar development. Statistics in Biopharmaceutical Research, 10(1), 39-49.
#'
#' Berger, R. L., & Hsu, J. C. (1996). Bioequivalence trials, intersection-union tests and equivalence confidence sets. Statistical Science, 283-302.
#'
#' @author
#' Johanna Muñoz \email{johanna.munoz@fromdatatowisdom.com}
#'
#' @examples
#'
#' mu_list <- list(SB2 = c(AUCinf = 38703, AUClast = 36862, Cmax = 127.0),
#'                 EUREF = c(AUCinf = 39360, AUClast = 37022, Cmax = 126.2),
#'                 USREF = c(AUCinf = 39270, AUClast = 37368, Cmax = 129.2))
#'
#' sigma_list <- list(SB2 = c(AUCinf = 11114, AUClast = 9133, Cmax = 16.9),
#'                    EUREF = c(AUCinf = 12332, AUClast = 9398, Cmax = 17.9),
#'                    USREF = c(AUCinf = 10064, AUClast = 8332, Cmax = 18.8))
#'
#' # Equivalent boundaries
#' lequi.tol <- c(AUCinf = 0.8, AUClast = 0.8, Cmax = 0.8)
#' uequi.tol <- c(AUCinf = 1.25, AUClast = 1.25, Cmax = 1.25)
#'
#' # Arms to be compared
#' list_comparator <- list(EMA = c("SB2", "EUREF"),
#'                         FDA = c("SB2", "USREF"))
#'
#' # Endpoints to be compared
#' list_y_comparator <- list(EMA = c("AUCinf", "Cmax"),
#'                           FDA = c("AUClast", "Cmax"))
#'
#' # Equivalence boundaries for each comparison
#' lequi_lower <- c(AUCinf = 0.80, AUClast = 0.80, Cmax = 0.80)
#' lequi_upper <- c(AUCinf = 1.25, AUClast = 1.25, Cmax = 1.25)
#'
#'# Run the simulation
#' sampleSize(power = 0.9, alpha = 0.05, mu_list = mu_list,
#'            sigma_list = sigma_list, list_comparator = list_comparator,
#'            list_y_comparator = list_y_comparator,
#'            list_lequi.tol = list("EMA" = lequi_lower, "FDA" = lequi_lower),
#'            list_uequi.tol = list("EMA" = lequi_upper, "FDA" = lequi_upper),
#'            adjust = "no", dtype = "parallel", ctype = "ROM", vareq = FALSE,
#'            lognorm = TRUE, ncores = 1, nsim = 50, seed = 1234)
#' @export
sampleSize <- function(mu_list, varcov_list = NA, sigma_list = NA, cor_mat = NA,
                       sigmaB = NA, Eper, Eco, rho = 0, TAR = NULL,
                       arm_names = NA,
                    ynames_list=NA,
                    type_y=NA,
                    list_comparator=NA,
                    list_y_comparator=NA,
                    power = 0.8,
                    alpha=0.05,
                    lequi.tol=NA,
                    uequi.tol=NA,
                    list_lequi.tol=NA,
                    list_uequi.tol=NA,
                    dtype="parallel",
                    ctype = "ROM",
                    vareq = TRUE,
                    lognorm = TRUE,
                    k=NA,
                    adjust="no",
                    dropout = NA,
                    nsim=5000,
                    seed=1234,
                    ncores=NA,
                    optimization_method = "fast",
                    lower=2,
                    upper=500,
                    step.power=6,
                    step.up=TRUE,
                    pos.side=FALSE,
                    maxiter = 1000, verbose = FALSE
){

  # Derive the Number of Arms
  n <- length(mu_list)

  # Assign default values for Eper and Eco
  if (missing(Eper)) {
    Eper <- c(0, 0)
    info_msg("Eper not provided. Defaulting to c(0, 0).", verbose)
  }
  if (missing(Eco)) {
    Eco <- c(0, 0)
    info_msg("Eco not provided. Defaulting to c(0, 0).", verbose)
  }

  # is mu provided?
  if (all(is.na(mu_list))) {
    stop("mu_list must be provided")
  }

  # Conduct validations
  validate_sample_size_limits(lower = lower, upper = upper)
  TAR <- derive_TAR(TAR = TAR, n_arms = n)


  # Derive the Arm Names
  arm_names <- derive_arm_names(arm_names = arm_names, mu_list = mu_list,
                                verbose = verbose)

  # Derive the Endpoint Names
  ynames_list <- derive_endpoint_names(ynames_list = ynames_list,
                                       mu_list = mu_list, verbose = verbose)

  # Derive the Treatment Allocation Rate
  TAR_list <- derive_allocation_rate(TAR = TAR, arm_names = arm_names,
                                     verbose = verbose)


  for (i in 1:n) {
    mu <- t(as.matrix(mu_list[[i]]))
    mu_list[[i]] <- mu
  }

  # Derive the list of covariance matrices
  varcov_list <- derive_varcov_list(mu_list = mu_list, sigma_list = sigma_list,
                                    ynames_list = ynames_list,
                                    varcov_list = varcov_list, cor_mat = cor_mat,
                                    rho = rho)




  weight_seq <- NA
  param.u <- list(mu = mu_list, varcov = varcov_list, TAR_list = TAR_list,
                  type_y = type_y, weight_seq = weight_seq, arm_names = arm_names,
                  ynames_list = ynames_list, list_comparator = list_comparator,
                  list_y_comparator = list_y_comparator)

  len_list <- c(length(mu_list), length(varcov_list), length(TAR_list)) # length in terms of arms

  if (max(len_list) != min(len_list)) {
    stop("'mu', 'varcov', and 'TAR' must be defined for all arms.")
  }

  # Remove endpoints with NA values in mu or varcov
  if(any(is.na(unlist(mu_list)))|any(is.na(unlist(varcov_list)))){
    for (i in 1:n){
      while(any(is.na(c(mu_list[[i]],as.vector(varcov_list[[i]]))))){
        mui <- mu_list[[i]]
        vari <- varcov_list[[i]]
        ynami <- ynames_list[[i]]
        ind_mu.na <- ind_cov.na <- NA
        mu.na <- which(is.na(mui))
        if(length(mu.na)!=0){
          ind_mu.na <- mu.na
        }
        cov.na <- which(is.na(vari), arr.ind = T)[,2]
        if(length(cov.na)!=0){
          ind_cov.na <- max(table(cov.na))
        }
        ind_na <- unique(c(ind_mu.na,ind_cov.na))
        ind_na <- ind_na[!is.na(ind_na)]
        ynames_list[[i]] <- ynami[-ind_na]
        mu_list[[i]] <- t(as.matrix(mui[,-ind_na]))
        vari <- vari[,-ind_na]
        vari <- vari[-ind_na,]
        varcov_list[[i]] <- vari
      }
    }
  }


  # name mu_list and varcov_list
  uynames <- NULL # unique ynames
  for (i in 1:n) {
    colnames(mu_list[[i]]) <- rownames(varcov_list[[i]]) <- colnames(varcov_list[[i]]) <- ynames_list[[i]]
    uynames <- c(uynames,ynames_list[[i]])
  }
  uynames <- unique(uynames)

  # length in terms of endpoints
  len_mu <- lapply(mu_list,length)
  len_cvar <- lapply(varcov_list,ncol)

  # test positive defined varcov
  validate_positive_definite(varcov_list)

  # Define weights according to type of endpoint,i.e. primary, secondary

  if (any(is.na(type_y)) | length(type_y) != length(uynames)) {
    type_y <- rep(1,length(uynames))
  }

  names(type_y) <- uynames

  weight <- 1/table(type_y)
  weight_seq <- type_y

  for (x in unique(type_y)) {
    weight_seq[weight_seq == x] <- weight[x]
  }
  names(weight_seq) <- uynames

  #if (len_mu[[1]] == 1){
  #  mu_list <- lapply(mu_list,FUN = function(x){array(unlist(x))})
  #  varcov_list <- lapply(varcov_list,FUN = function(x){matrix(unlist(x))})}

  names(mu_list) <-  names(varcov_list) <- names(TAR_list) <- names(ynames_list) <- arm_names

  # Get list of comparators if it is not provided
  if (any(is.na(list_comparator))) {
    comb_mat <- utils::combn(arm_names,2)
    list_comparator <- list()
    for (i in 1:ncol(comb_mat)) {
      list_comparator[[i]] <- c(comb_mat[1,i],comb_mat[2,i])
    }
  }

  if (any(!unique(unlist(list_comparator)) %in% arm_names)) {
    stop("All arm names specified in 'list_comparator' must be present in 'arm_names'.")
  }

  # Get the endpoints to be compared on each comparator
  if (any(is.na(list_y_comparator)) | length(list_comparator) != length(list_y_comparator)) {
    list_y_comparator <- list()
  }

  for (i in 1:length(list_comparator)) {
    treat1 <- list_comparator[[i]][[1]]
    treat2 <- list_comparator[[i]][[2]]
    y_comp <- intersect(colnames(mu_list[[treat1]]),colnames(mu_list[[treat2]]))

    if (i > length(list_y_comparator)) { # No element assigned
      # warning("As no list_y_comparator was provided, it will be compared all endpoints in commun on the compared arms")
      list_y_comparator[[i]] <- y_comp
    }

    if(any(!list_y_comparator[[i]]%in%y_comp)){
      warning(paste0("It will be compared only the endpoints included both arms of the comparator",paste(list_comparator[[i]],collapse="")))
      y_compi <- list_y_comparator[[i]]
      if(length(y_compi[y_compi%in%y_comp])==0){
        list_y_comparator[[i]] <- NULL
        list_comparator[[i]] <- NULL}
      else{
        list_y_comparator[[i]] <- y_compi[y_compi%in%y_comp]
      }
    }
  }


  if(any(unlist(len_mu) != unlist(len_cvar))){
    stop("mu,varcov should be defined for all the endpoints")
  }

  # list equivalence boundaries

  # check if list or equitol vector are not provided
  if ((all(is.na(lequi.tol))&all(is.na(uequi.tol)))&
      (all(is.na(list_lequi.tol))&all(is.na(list_uequi.tol)))
      ){
    warning("No boundaries were provided so standard values will be used")

  }

  # when only equitol vector is provided
  if(!any(is.na(c(lequi.tol,uequi.tol)))){
    if(any(lequi.tol>=uequi.tol)){
      warning("some lequitol>=uequi.tol, so reference values will be used")
      lequi.tol <- NA
      uequi.tol <- NA
    }
    if ((length(lequi.tol) != length(uynames)) | (length(uequi.tol) != length(uynames))) {
      warning("Insufficient number of tolerance values supplied (One needed for each endpoint), so reference values will be used for all comparators")
      lequi.tol <- NA
      uequi.tol <- NA
    }
  }

  if (all(is.na(list_lequi.tol),is.na(list_uequi.tol)) & any(!is.na(c(lequi.tol,uequi.tol)))){
    warning("Using the same tolerance boundaries (lequi.tol, uequi.tol) across all comparators.")
    list_lequi.tol <- list()
    list_uequi.tol <- list()
      if (any(lequi.tol >= uequi.tol)) {
        warning("Inconsistent tolerance bounds: some values in lequi.tol are greater than or equal to uequi.tol. Reference values will be used instead.")
        lequi.tol <- NA
        uequi.tol <- NA
      }
      if ((length(lequi.tol) != length(uynames)) | (length(uequi.tol) != length(uynames))) {
        warning("Insufficient number of tolerance values provided: one value per endpoint is required. Reference values will be used instead.")
        lequi.tol <- NA
        uequi.tol <- NA
      }

    for (i in 1:length(list_comparator)){
      muend <-  mu_list[[list_comparator[[i]][[2]]]]
      if ((length(lequi.tol) == length(muend)) & (length(uequi.tol) == length(muend)) & all(!is.na(lequi.tol) & !is.na(uequi.tol))){
        lequi.toli <- lequi.tol
        uequi.toli <- uequi.tol
      }else if(ctype == "DOM"){
        lequi.toli <- - 0.2*muend
        uequi.toli <- 0.2*muend
      }else{
        lequi.toli <- rep(0.80,length(muend))
        uequi.toli <- rep(1.25,length(muend))
      }

      if(is.null(names(lequi.toli))|is.null(names(uequi.toli))){
        lequi.toli<- as.vector(lequi.toli)
        uequi.toli<- as.vector(uequi.toli)
        names(lequi.toli)<-names(uequi.toli)<-colnames(muend)
      }
      list_lequi.tol[[i]] <-lequi.toli
      list_uequi.tol[[i]] <-uequi.toli

    }

  }


  for (i in 1:length(list_comparator)){
    muend <-  mu_list[[list_comparator[[i]][[2]]]]
    if (is.null(names(list_lequi.tol[[i]]))|is.null(names(list_uequi.tol[[i]]))){
      names(list_lequi.tol[[i]])<- names(list_uequi.tol[[i]])<-colnames(muend)
    }
  }

  # check list equitol

  for (i in 1:length(list_comparator)){
    muend <-  mu_list[[list_comparator[[i]][[2]]]]

    namerep <- unique(names(which((list_lequi.tol[[i]]>=list_uequi.tol[[i]])|is.na(list_lequi.tol[[i]])| is.na(list_uequi.tol[[i]])))) # name to be replaced by reference

    if(length(namerep)>0){
      list_lequi.tol[[i]][namerep]
      if(ctype == "DOM"){
        list_lequi.tol[[i]][namerep] <- - 0.2*muend[namerep]
        list_uequi.tol[[i]][namerep] <- 0.2*muend[namerep]
      }else{
        list_lequi.tol[[i]][namerep] <- 0.8
        list_uequi.tol[[i]][namerep] <- 1.25
      }
    }
  }

  if( any(is.na(dropout))){
    if(dtype=="parallel"){
      dropout <- rep(0,length(arm_names))
      names(dropout) <- arm_names
    }else{
      dropout <- rep(0,2)
    }
  }


  if (dtype == "parallel") {
    if (length(arm_names) != length(dropout)) {
      warning("The number of dropout values provided does not match the number of arms specified in 'arm_names'. A default dropout rate of 0 will be assigned to each arm.")
      dropout <- rep(0,length(arm_names))
    }
    if (is.null(names(dropout))) {
      names(dropout) <- arm_names
    }

    # Check if dtype is "parallel" and Eper or Eco are non-default
    if (any(Eper != c(0, 0)) || any(Eco != c(0, 0))) {
      warning("Eper and Eco are only applicable for dtype = '2x2'. Non-default values for Eper or Eco will be ignored in parallel design.")
    }
  } else {
    if(length(dropout)!=2){
      warning("Incorrect number of dropout supplied (One needed for each sequence),so it will be assigned a dropout=0")
      dropout <- rep(0,2)
    }
  }



  # k is na
  kmax <- sapply(list_y_comparator,length)

  if(any(is.na(k))){
    #warning("No k vector provided,it will be set to the total number of endpoints on each comparator")
    k <- kmax
  }

  if(length(k)!=length(kmax)){
    #warning("No k vector provided,it will be set to the total number of endpoints on each comparator")
    k <- rep(k[[1]],length(kmax))
  }

  if (any(k>kmax)){
    warning("The specified k on a comparator is larger than the number of endpoints,
              k is assigned to the total number of endpoints in the comparator")
    k <- pmin(k,kmax)
  }


  # Save endpoints related information on a parameter list

  param <- list(mu = mu_list, varcov = varcov_list, sigmaB = sigmaB,
                TAR_list = TAR_list, type_y = type_y, weight_seq = weight_seq,
                arm_names = arm_names,  ynames_list = ynames_list,
                list_comparator = list_comparator,
                list_y_comparator = list_y_comparator,
                list_lequi.tol = list_lequi.tol,
                list_uequi.tol = list_uequi.tol,
                sigmaB = sigmaB,
                Eper = Eper, Eco = Eco)

  # Parameters related to design ----

  if (lognorm == TRUE & ctype == "DOM"){
    stop("Testing is not supported for DOM when variables follow a log-normal distribution.")
  }

  if( lognorm == FALSE & ctype == "ROM" & dtype =="2x2"){
    #stop("Test not available here")
  }

  if(!(all(unlist(list_comparator)%in%arm_names)|all(arm_names%in%unlist(list_comparator)))){
    stop("Names included in the list_comparator should be coherent with the names in the arm_names vector")
  }

  param.d <- list(nsim=nsim,
                  power=power,
                  alpha=alpha,
                  dtype=dtype,
                  ctype=ctype,
                  lognorm=lognorm,
                  vareq=vareq,
                  k=k,
                  adjust=adjust,
                  dropout=dropout,list_lequi.tol=list_lequi.tol, list_uequi.tol=list_uequi.tol)

  if(is.na(ncores)) {
    ncores <- parallel::detectCores() - 1
  }


  if (optimization_method == "fast") {
    opt.response <- uniroot.integer.mod(function(x) (power_cal(n=x,nsim=nsim,param=param,param.d=param.d,seed=seed,ncores=ncores)),
                                        power = power,
                                        lower = lower,
                                        upper = upper,
                                        step.power = step.power,
                                        step.up = step.up,
                                        pos.side = pos.side,
                                        maxiter = maxiter)

    table.test <- data.table::as.data.table(opt.response$table.test)
  } else if (optimization_method == "step-by-step") {
    table.test <- NULL
    for (n in lower:upper) { # increase one by one until the desired power or the
      # maximal sample size is reached

      powercal <- power_cal(n=n,nsim=nsim,param=param,param.d=param.d,seed=seed,ncores=ncores)
      if( any(is.na(powercal))){
        powercal <- 0
      }
      if (powercal$power > power){
        table.test <- rbind(table.test,powercal$output.test)
        break
      }
      table.test <- rbind(table.test,powercal$output.test)
    }

    table.test <- data.table::as.data.table(table.test)
  } else {
    stop ("Invalid search way")
  }


  # Calculate totaly test across all comparators= power
  qnam <- colnames(table.test)[grep("^totaly",colnames(table.test))]
  table.test$totaly <- apply(table.test[,qnam,with=FALSE], 1, prod)


  # Get a summary across all the n_iter with confidence interval power
  n_iter <- NULL
  namexc <- colnames(table.test)[grep("^[^(mu_|sd_|eql_|equ_)]",colnames(table.test))]
  summary <- table.test[, lapply(.SD, FUN=function(x){sum(x, na.rm=TRUE)/nsim}), by= n_iter ][,c(namexc),with=FALSE]
  powerfun <- function(x) {
    bin_test <- stats::prop.test(x = x, n = nsim, correct = TRUE)
    c(bin_test$estimate[[1]],bin_test$conf[1],bin_test$conf[2])
  }

  powerv <- do.call(rbind,lapply(summary$t_true, powerfun))
  colnames(powerv) <- c("power","power_LCI","power_UCI")
  summary <- cbind(summary,powerv)
  sumcol <- colnames(summary)[grep("^(n_|power)",colnames(summary))]
  table.iter <- summary[,sumcol,with=FALSE]

  if (optimization_method == "fast") {
    if (is.null(opt.response$power)){
      response <- NA
    } else{
      response <- table.iter[n_iter == opt.response$power["n_iter"],]
    }

    out <- list( response = response,
                 table.iter = table.iter,
                 table.test = table.test,
                 param.u = param.u,
                 param = param,
                 param.d = param.d)
  } else if (optimization_method == "step-by-step") {
    out <- list(   response = table.iter[n_iter == n,],
                   table.iter = table.iter,
                   table.test = table.test,
                   param = param,
                   param.d = param.d)
  }

  class(out) <- "simss"
  return(out)

}

#' Derive or Assign Arm Names
#'
#' This function checks if `arm_names` is provided. If `arm_names` is missing, it attempts to derive names
#' from `mu_list`. If `mu_list` does not contain names, it assigns default names ("A1", "A2", etc.) to each arm.
#' Informational messages are displayed if `verbose` is set to `TRUE`.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @param arm_names Optional vector of arm names.
#' @param mu_list Named list of means per treatment arm, from which arm names may be derived.
#' @param verbose Logical, if `TRUE`, displays messages about the derivation process.
#'
#' @return A vector of arm names.
derive_arm_names <- function(arm_names, mu_list, verbose = FALSE) {

  # Check if arm_names is missing and attempt to derive from mu_list
  if (any(is.na(arm_names))) {
    if (!is.null(names(mu_list))) {
      arm_names <- names(mu_list)
      info_msg(paste("Arm names derived from mu_list: ", paste(arm_names, collapse = ", ")), verbose)
    } else {
      arm_names <- paste0("A",seq(mu_list))
      info_msg(paste("Arm names not provided and could not be derived from mu_list. Assigning default names: ", paste(arm_names, collapse = ", ")), verbose)
    }
  } else {
    info_msg(paste("Using user-provided arm names: ", paste(arm_names, collapse = ", ")), verbose)
  }

  return(arm_names)
}

#' Derive Endpoint Names
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' This function derives endpoint names (`ynames_list`) from `mu_list` if `ynames_list`
#' is missing. If `ynames_list` is already provided, it confirms the names to the user when
#' `verbose` is set to `TRUE`.
#'
#' @param ynames_list Optional list of vectors with endpoint names for each arm.
#' @param mu_list Named list of means per treatment arm, where names can be used as endpoint names.
#' @param verbose Logical, if `TRUE`, displays messages about the derivation process.
#'
#' @return A list of endpoint names for each arm.
derive_endpoint_names <- function(ynames_list, mu_list, verbose = FALSE) {

  # Check if ynames_list is missing and attempt to derive from mu_list
  if (any(is.na(ynames_list))) {

    # Try to derive the ynames from mu_list
    ynames_list <- lapply(mu_list, function(x) names(x))
    info_msg("Attempting to derive endpoint names (ynames_list) from mu_list.", verbose)

    # Check if ynames were successfully derived
    if (length(names(ynames_list)) == 0 || any(sapply(ynames_list, is.null))) {
      info_msg("Not all endpoint names were provided. Assigning arbitrary names (y1, y2, etc.) to endpoints for each arm.", verbose)
      ynames_list <- lapply(mu_list, function(x) paste0("y", 1:length(x)))
    } else {
      info_msg("Endpoint names derived from mu_list.", verbose)
    }
  } else {
    info_msg("Using user-provided endpoint names (ynames_list).", verbose)
  }

  return(ynames_list)
}

#' Derive Treatment Allocation Rate (TAR)
#'
#' This function checks if `TAR` (treatment allocation rate) is provided. If `TAR` is missing, it assigns a default
#' equal allocation rate across all arms and ensures the `TAR` values correspond to the specified `arm_names`. It then
#' converts `TAR` to a named list for further use. Informational messages are displayed if `verbose` is set to `TRUE`.
#'
#' @param TAR Optional numeric vector specifying the allocation rate for each treatment arm. If missing, a default
#' equal allocation rate is assigned.
#' @param arm_names Character vector specifying the names of the treatment arms. Used to name the elements of `TAR`.
#' @param verbose Logical, if `TRUE`, displays messages about the status of `TAR` derivation or assignment.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @return A list representing the treatment allocation rate for each arm.
derive_allocation_rate <- function(TAR = NULL, arm_names, verbose = FALSE) {

  n_arms <- length(arm_names)

  # Check if TAR is missing and assign default if necessary
  if (any(is.na(TAR)) | is.null(TAR)) {
    TAR <- rep(1, n_arms)  # Default equal allocation across all arms
    names(TAR) <- arm_names  # Assign names to correspond with arm names
    info_msg(paste("TAR not provided. Assigning equal allocation rate across all arms: ", paste(TAR, collapse = ":")), verbose)
  } else {
    names(TAR) <- arm_names  # Ensure TAR values are named to match arm names
    info_msg(paste("Using user-provided TAR: ", paste(TAR, collapse = ":")), verbose)
  }

  # Convert TAR to list format
  TAR_list <- as.list(TAR)

  return(TAR_list)
}

#' Derive Variance-Covariance Matrix List
#'
#' Constructs a list of variance-covariance matrices for multiple treatment arms based on provided standard deviations,
#' means, and correlation structures.
#'
#' @param mu_list A list of numeric vectors representing the means (\eqn{\mu}) for each treatment arm. Each element corresponds to one arm.
#' @param sigma_list A list of numeric vectors representing the standard deviations (\eqn{\sigma}) for each treatment arm. Each element corresponds to one arm.
#' @param ynames_list A list of character vectors specifying the names of the endpoints for each arm. Each element corresponds to one arm.
#' @param varcov_list (Optional) A pre-specified list of variance-covariance matrices for each arm. If provided, it will override the construction of variance-covariance matrices.
#' @param cor_mat (Optional) A correlation matrix to be used for constructing the variance-covariance matrices when there are multiple endpoints. If dimensions do not match the number of endpoints, a warning is issued.
#' @param rho (Optional) A numeric value specifying the constant correlation coefficient to be used between all pairs of endpoints if no correlation matrix is provided. Default is 0 (uncorrelated endpoints).
#'
#' @details
#' This function creates a list of variance-covariance matrices for multiple treatment arms. If the \code{varcov_list} is not provided,
#' the function uses the \code{sigma_list} to compute the matrices. For single endpoints, the variance is simply the square of the standard deviation.
#' For multiple endpoints, the function constructs the matrices using either a provided \code{cor_mat} or the constant correlation coefficient \code{rho}.
#'
#' The function ensures that the lengths of \code{mu_list}, \code{sigma_list}, and \code{ynames_list} match for each arm. If dimensions mismatch,
#' or if neither a variance-covariance matrix (\code{varcov_list}) nor a standard deviation list (\code{sigma_list}) is provided, an error is raised.
#'
#' @return A list of variance-covariance matrices, one for each treatment arm.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
derive_varcov_list <- function(mu_list, sigma_list, ynames_list, varcov_list = NULL, cor_mat = NULL, rho = 0) {
  # Check if variance-covariance matrix is missing
  if (any(is.na(varcov_list))) {
    if (any(is.na(sigma_list))) {
      stop("No variance-covariance matrix provided, and a standard deviation list is also missing. Either a variance-covariance matrix or a standard deviation list is required.")
    }

    # Validate lengths of inputs
    len_mu <- sapply(mu_list, length)
    len_sd <- sapply(sigma_list, length)
    len_y <- sapply(ynames_list, length) # Number of endpoints

    if (any((len_mu != len_sd) | (len_mu != len_y))) {
      stop("In each arm, 'mu', 'sigma', and 'y_name' must have the same length.")
    }

    # Initialize the varcov_list
    varcov_list <- vector("list", length(mu_list))

    # Loop through each arm to construct the variance-covariance matrices
    for (i in seq_along(sigma_list)) {
      m <- length(sigma_list[[i]]) # Number of endpoints

      if (m == 1) {
        # Single endpoint: variance is the square of the standard deviation
        varcov <- as.matrix(sigma_list[[i]]^2)
      } else {
        # Multiple endpoints: construct the correlation matrix
        R <- matrix(rho, m, m)
        diag(R) <- 1 # Set diagonal to 1 (self-correlation)

        # Use the provided correlation matrix if dimensions match
        if (any(!is.na(cor_mat))) {
          if ((nrow(cor_mat) == m) & (ncol(cor_mat) == m)) {
            R <- as.matrix(cor_mat)
          } else {
            warning("An uncorrelated matrix will be used as the provided matrix does not have the expected dimensions.")
          }
        }

        # Construct the variance-covariance matrix
        varcov <- diag(sigma_list[[i]]) %*% R %*% diag(sigma_list[[i]])
      }

      # Add the matrix to the list
      varcov_list[[i]] <- varcov
    }
  }

  return(varcov_list)
}

#' Derive Treatment Allocation Rate (TAR)
#'
#' This function validates and adjusts the treatment allocation rate (`TAR`) to ensure it matches
#' the number of treatment arms (`n_arms`). If `TAR` is missing or NULL, it is set to a default
#' vector of ones. If `TAR` is shorter than `n_arms`, missing values are replaced with ones.
#' If `TAR` contains NA values, they are replaced with ones. If `TAR` exceeds `n_arms` in length
#' or contains non-positive values, an error is raised.
#'
#' @param TAR Numeric vector. Treatment allocation rates for each arm.
#' @param n_arms Numeric. The number of treatment arms that `TAR` must correspond to.
#'
#' @details
#' - If `TAR` is NULL or missing, it is set to a vector of ones of length `n_arms`.
#' - If `TAR` is shorter than `n_arms`, missing values are replaced with ones and a warning is issued.
#' - If `TAR` is longer than `n_arms`, an error is raised.
#' - If `TAR` contains NA values, they are replaced with ones and a warning is issued.
#' - If `TAR` contains negative or zero values, an error is raised.
#'
#' @author
#' Thomas Debray \email{tdebray@fromdatatowisdom.com}
#'
#' @return A numeric vector of length `n_arms`, where any missing or NA values in `TAR`
#' have been replaced with ones.
derive_TAR <- function(TAR = NULL, n_arms) {

  if (missing(TAR) || is.null(TAR)) {
    warning("Warning: TAR is missing or NULL. Setting TAR to a default vector of ones.")
    TAR <- rep(1, n_arms)
  }

  # Ensure TAR has the correct length by filling missing values with 1
  if (length(TAR) < n_arms) {
    warning("Warning: TAR length is shorter than the number of arms. Missing values will be replaced with 1.")
    TAR <- c(TAR, rep(1, n_arms - length(TAR)))
  } else if (length(TAR) > n_arms) {
    stop("Validation Error: The length of TAR must not exceed the number of arms specified by n_arms.")
  }

  # Replace any NA values with 1
  if (any(is.na(TAR))) {
    warning("Warning: NA values detected in TAR. These will be replaced with 1.")
    TAR[is.na(TAR)] <- 1
  }

  # Validate that all values are positive
  if (any(TAR <= 0, na.rm = TRUE)) {
    stop("Validation Error: TAR must contain only positive values. Negative or zero values are not allowed.")
  }

  return(TAR)
}

