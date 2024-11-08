  #' @title Sample Size Calculation for Target Power step by step
  #' @description  Alternative way to estimate the sample size calculation given a target power (step by step)----
  #'
  #' @param mu_list Named list of arithmetic means per treatment arm. Each element contains a vector (i.e., one per treatment arm) with the expected outcomes for all endpoints of interest.
  #' @param varcov_list list of var-cov matrices, each element corresponds to a comparator with a varcov matrix of size number of endpoints X number of endpoints.
  #' @param sigma_list  list of sigma vectors, each element corresponds to a comparator with a sigma vector of size number of endpoints.
  #' @param cor_mat matrix specifying the correlation matrix between endpoints, used along with sigma_list  to calculate the varcov list in case it is not provided.
  #' @param sigmaB number between subject variance only for 2x2 design.
  #' @param Eper  vector of size 2, effect of period on c(0,1).
  #' @param Eco vector of size 2, carry over effect of arm c(Reference, Treatment).
  #' @param rho correlation parameter to fix on all the pair of endpoints, used along with sigma_list to calculate the varcov list in case neither cor_mat or var-cov list is not provided.
  #' @param TAR vector of allocation rates with allocation rates of the arm, default is equivalent rate.
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
  #' @param vareq variance equivalence assumption (FALSE, TRUE)
  #' @param dtype design type ("parallel","2x2")
  #' @param lognorm Is data log-normally distributed? (TRUE, FALSE)
  #' @param k Vector with the number of endpoints that must be successful (integer) for global bioequivalence for each comparator. If no k vector is provided, it will be set to the total number of endpoints on each comparator.
  #' @param adjust alpha adjustment ( "k", "bon","sid","no","seq")
  #' @param ctype comparison dtype ("DOM"(Difference of means), "ROM"(Ratio of means))
  #' @param dropout vector with proportion of total population with dropout per arm
  #' @param nsim number of simulated studies (default=5000)
  #' @param lower initial value of N to be searched (default=2)
  #' @param upper max value of N to be searched (default=500)
  #' @param seed main seed
  #' @param ncores Number of processing cores for parallel computation; defaults to the total detected cores minus one.
  #'
  #' @return An object simss that contains the following elements :
  #' \describe{
  #'  \item{"response"}{ array with the sample sizes for each arm and aproximated achieved power with confidence intervals}
  #'  \item{"table.iter"}{data frame with the estimated sample size for each arm and power calculated at each searching iteration}
  #'  \item{"table.test"}{data frame that collects the total information of the simulation at each iteration}
  #'  \item{"param.u"}{parameters provided by the user}
  #'  \item{"param"}{ parameters used for the sample size calculation; as param.u are checked and modified in case of any inconsistent or missing information provided}
  #'  \item{"param.d"}{ parameters of design}
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
  #' #' @examples
  #'
  #' mu_list <- list(SB2 = c(AUCinf = 38703, AUClast = 36862, Cmax = 127.0),
  #'                 EUREF = c(AUCinf = 39360, AUClast = 37022, Cmax = 126.2),
  #'                 USREF = c(AUCinf = 39270, AUClast = 37368, Cmax = 129.2))
  #'
  #' sigma_list <- list(SB2 = c(AUCinf = 11114, AUClast = 9133, Cmax = 16.9),
  #'                    EUREF = c(AUCinf = 12332, AUClast = 9398, Cmax = 17.9),
  #'                    USREF = c(AUCinf = 10064, AUClast = 8332, Cmax = 18.8))
  #'
  #'# Equivalent boundaries
  #'lequi.tol <- c(AUCinf = 0.8, AUClast = 0.8, Cmax = 0.8)
  #'uequi.tol <- c(AUCinf = 1.25, AUClast = 1.25, Cmax = 1.25)
  #'
  #'
  #' # arms to be compared
  #' list_comparator <- list(EMA = c("SB2", "EUREF"),
  #'                         FDA = c("SB2", "USREF"))
  #'
  #'# Endpoints to be compared
  #'list_y_comparator <- list(EMA = c("AUCinf", "Cmax"),
  #'                          FDA = c("AUClast", "Cmax"))
  #'
  #'# Run the simulation
  #'result <- calsbs(power = 0.9, # target power
  #'                 alpha = 0.05,
  #'                 mu_list = mu_list,
  #'                 sigma_list = sigma_list,
  #'                 lequi.tol = lequi.tol,
  #'                 uequi.tol = uequi.tol,
  #'                 list_comparator = list_comparator,
  #'                 list_y_comparator = list_y_comparator,
  #'                 adjust = "no",
  #'                 dtype = "parallel",
  #'                 ctype = "ROM",
  #'                 vareq = FALSE,
  #'                 lognorm = TRUE,
  #'                 ncores = 1,
  #'                 lower= 4,
  #'                 upper = 100,
  #'                 nsim = 50,
  #'                 seed = 1234)
  #'
  #' @export
  calsbs <- function( mu_list,
                      varcov_list=NA,
                      sigma_list=NA,
                      cor_mat=NA,
                      sigmaB =NA,
                      Eper = c(0,0),
                      Eco = c(0,0),
                      rho=0,
                      TAR=NA,
                      arm_names=NA,
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
                      vareq=T,
                      lognorm=T,
                      k=NA,
                      adjust="no",
                      dropout = NA,
                      nsim=5000,
                      lower=2,
                      upper=500,
                      seed=1234,
                      ncores=NA
  ){

    # is mu provided?
    if (all(is.na(mu_list))) {
      stop("mu_list must be provided")
    }

    # Parameters endpoints -----
    n <- length(mu_list) # number of arms

    # Derive the arm names
    if (any(is.na(arm_names))) {
      if (!is.null(names(mu_list))) {
        arm_names <- names(mu_list)
      }
    }

    # Derive the names of the endpoints in each arm
    #
    # Example output for a trial with 2 treatments and 3 outcomes:
    # $SB2
    # [1] "AUCinf"  "AUClast" "Cmax"
    # $EUREF
    # [1] "AUCinf"  "AUClast" "Cmax"
    if (any(is.na(ynames_list))) {

      # Try to derive the ynames from mu_list
      ynames_list <- lapply(mu_list, function(x) names(x))

      if (length(names(ynames_list)) == 0 | any(sapply(ynames_list, is.null))) {
        #warning("no all endpoints names provided for each arm, so arbitrary names are assigned")
        ynames_list <- lapply(mu_list, function(x) paste0("y", 1:length(x)))
      }
    }


    # Treatment allocation rate
    if (any(is.na(TAR))) {
      TAR <- rep(1,n)
    }

    TAR_list <- as.list(TAR)


    for (i in 1:n) {
      mu <- t(as.matrix(mu_list[[i]]))
      mu_list[[i]] <- mu
    }


    # Varcov specfication
    if (any(is.na(varcov_list))) {
      if (any(is.na(sigma_list))) {
        stop("No variance-covariance matrix provided, and a standard deviation list is also missing. Either a variance-covariance matrix or a standard deviation list is required.")
      }

      # length in terms of endpoints
      len_mu <- sapply(mu_list,length)
      len_sd <- sapply(sigma_list,length)
      len_y <- sapply(ynames_list,length) # number of endpoints

      if (any((len_mu != len_sd) | (len_mu != len_y))) {
        stop("In each arm, 'mu', 'sigma', and 'y_name' must have the same length.")
      }

      varcov_list <- NULL

      for (i in 1:n) {
        m <- length(sigma_list[[i]]) # number of endpoints

        if (m == 1) {
          varcov <- as.matrix(sigma_list[[i]]^2)

        } else {
          R <- matrix(rho, m, m)
          diag(R) <- 1 # correlation matrix

          if (any(!is.na(cor_mat))) {
            if ((nrow(cor_mat) == m) & (ncol(cor_mat) == m)) {
              R <- as.matrix(cor_mat)
            } else {
              warning("An uncorrelated matrix will be used as the provided matrix does not have the expected dimensions.")
            }
          }

          varcov <- diag(sigma_list[[i]]) %*% R %*% diag(sigma_list[[i]])
        }

        varcov_list[[i]] <-  varcov
      }
    }

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
    positive <- function(x) {
      resp <- base::tryCatch({matrixcalc::is.positive.semi.definite(round(x,3))},
                             error = function(e) {FALSE})
    }


    lis_pdef <- unlist(lapply(varcov_list, positive))

    if (!all(lis_pdef)) {
      stop("All 'varcov' matrices must be symmetric and positive definite.")
    }

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

    # Give to list the arm names
    if (any(is.na(arm_names))) {
      arm_names <- paste0("A",rep(1:n))
    }

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
    if (all(all(is.na(lequi.tol)),all(is.na(uequi.tol)),all(is.na(list_lequi.tol)),all(is.na(list_uequi.tol)))){
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
        warning("Insufficient number of tolerance values supplied (One needed for each endpoint), so reference values will be used")
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


    if(dtype=="parallel"){
      if(length(arm_names)!=length(dropout)){
        warning("Incorrect number of dropout supplied (One needed for each arm),so it will be assigned a dropout=0")
        dropout <- rep(0,length(arm_names))
      }
      if(is.null(names(dropout))){
        names(dropout) <- arm_names}
    }else{
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


    # Save endopoints related information on a parameter list

    param <- list(mu = mu_list, varcov = varcov_list, sigmaB=sigmaB, TAR_list = TAR_list, type_y = type_y, weight_seq = weight_seq, arm_names=arm_names,  ynames_list =ynames_list, list_comparator = list_comparator, list_y_comparator= list_y_comparator,sigmaB=sigmaB,Eper=Eper,Eco=Eco)

    # Parameters related to design ----

    if (lognorm == TRUE & ctype == "DOM"){
      stop("No test available for DOM when variables are log normal distributed")
    }

    if( lognorm == FALSE & ctype == "ROM" & dtype =="2x2"){
      #stop("Test not available here")
    }

    if(!(all(unlist(list_comparator)%in%arm_names)|all(arm_names%in%unlist(list_comparator)))){
      stop("Names included in the list_comparator should be coherent with the names in the arm_names vector")
    }

    param.d <- list(nsim=nsim,power=power,alpha=alpha,dtype=dtype,ctype=ctype,lognorm=lognorm,vareq=vareq,k=k,adjust=adjust,dropout=dropout,list_lequi.tol=list_lequi.tol, list_uequi.tol=list_uequi.tol)

    if (is.na(ncores)){
      ncores <- parallel::detectCores() - 1}



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

  # Calculate totaly test across all comparators= power
  qnam <- colnames(table.test)[grep("^totaly",colnames(table.test))]
  table.test$totaly <- apply(table.test[,qnam,with=FALSE], 1, prod)

  # Get a summary across all the n_iter with confidence interval power
  n_iter= NULL
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


  out <- list(   response = table.iter[n_iter==n,],
                 table.iter = table.iter,
                 table.test = table.test,
                 param = param,
                 param.d = param.d)

  class(out) <- "simss"
  return(out)
}
