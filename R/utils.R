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
#' @param param.d design parameters
#'
#' @return power calculated from a global list of comparators
#' @keywords internal
#'
.normalize_internal_distribution <- function(param.d) {
  if (is.null(param.d$distribution) || length(param.d$distribution) == 0L ||
      is.na(param.d$distribution[[1L]])) {
    param.d$distribution <- if (isTRUE(param.d$lognorm)) "lnorm" else "norm"
  }
  param.d
}

power_cal <- function(n,nsim,param,param.d,seed,ncores){

  param.d <- .normalize_internal_distribution(param.d)

  if (param.d$dtype == "parallel") {
    TAR_used <- unlist(param$TAR_list)[unique(unlist(param$list_comparator))]
    size <- ceiling(n*TAR_used)
    size[size < 2] <- 2
    size_ndrop <- ceiling((1 - param.d$dropout[names(size)])*size)
    size_ndrop[size_ndrop < 2] <- 2
    n_drop <- sum(size) - sum(size_ndrop)

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

# Normalize the user-facing multiplicity-adjustment labels.
.normalize_adjustment <- function(adjust, allow_sequential = TRUE) {
  if (!is.character(adjust) || length(adjust) != 1L || is.na(adjust) ||
      !nzchar(trimws(adjust)))
    stop("'adjust' must be a single character value.")
  value <- tolower(trimws(adjust))
  aliases <- c(
    no = "no", none = "no",
    bon = "bon", bonferroni = "bon",
    sid = "sid", sidak = "sid",
    k = "k",
    t = "t", pc = "t", partial = "t", partial_conjunction = "t",
    "partial-conjunction" = "t", "t-adjustment" = "t",
    seq = "seq", sequential = "seq"
  )
  if (!value %in% names(aliases))
    stop("'adjust' must be one of 'no', 'bon', 'sid', 'k', 't', or 'seq'.")
  result <- unname(aliases[[value]])
  if (!allow_sequential && result %in% c("k", "seq"))
    stop("This outcome type supports only 'no', 'bon', 'sid', and 't' adjustment.")
  result
}

# Derive the endpoint-wise significance level for one comparator. The endpoint
# count is deliberately supplied by the caller: it is the number of endpoints
# actually tested for that comparator, not the number stored in arm parameters.
.endpoint_alpha <- function(alpha, m, k, adjust) {
  if (length(m) != 1L || !is.finite(m) || m < 1L || m != as.integer(m))
    stop("'m' must be a positive integer.")
  if (length(k) != 1L || !is.finite(k) || k < 1L || k > m ||
      k != as.integer(k))
    stop("'k' must be an integer between 1 and 'm'.")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a number strictly between 0 and 1.")
  switch(adjust,
    no = rep(alpha, m),
    bon = rep(alpha / m, m),
    sid = rep(1 - (1 - alpha)^(1 / m), m),
    k = rep(k * alpha / m, m),
    t = rep(alpha / (m - k + 1), m),
    seq = rep(alpha, m),
    stop("Unknown multiplicity adjustment: ", adjust)
  )
}

# Validate and align the hierarchy vector after comparator endpoints have been
# resolved. Named vectors may cover only endpoints actually used in a
# comparison family; endpoints outside that family are assigned type 1 but are
# never tested by the sequential branch.
.prepare_type_y <- function(type_y, all_endpoints, selected_endpoints,
                            adjust) {
  supplied <- !(is.null(type_y) ||
                (length(type_y) == 1L && is.na(type_y)))
  # -1 is the package's internal sentinel for "no hierarchy" in objects
  # created by older versions; it should not trigger an ignored-input warning.
  if (supplied && !identical(adjust, "seq") && is.numeric(type_y) &&
      all(!is.na(type_y) & type_y < 0))
    supplied <- FALSE
  all_endpoints <- as.character(all_endpoints)
  selected_endpoints <- unique(as.character(selected_endpoints))
  default <- stats::setNames(rep(-1L, length(all_endpoints)), all_endpoints)

  if (!identical(adjust, "seq")) {
    if (supplied) {
      warning("'type_y' is used only with adjust = 'seq'; the supplied hierarchy is ignored.",
              call. = FALSE)
    }
    return(list(type_y = default, supplied = supplied, active = FALSE))
  }

  if (!supplied) {
    warning("adjust = 'seq' requires 'type_y'; no hierarchical alpha adjustment will be applied.",
            call. = FALSE)
    return(list(type_y = default, supplied = FALSE, active = FALSE))
  }
  if (!is.numeric(type_y) || any(!is.finite(type_y)) ||
      any(type_y != as.integer(type_y)) || any(!type_y %in% c(1, 2)))
    stop("When adjust = 'seq', 'type_y' must contain only integer values 1 (primary) or 2 (secondary).")
  if (!is.null(names(type_y))) {
    if (anyNA(names(type_y)) || any(!nzchar(names(type_y))) ||
        anyDuplicated(names(type_y)))
      stop("Names of 'type_y' must be non-empty and unique.")
    unknown <- setdiff(names(type_y), all_endpoints)
    if (length(unknown))
      stop("'type_y' contains unknown endpoint(s): ", paste(unknown, collapse = ", "))
    missing <- setdiff(selected_endpoints, names(type_y))
    if (length(missing))
      stop("'type_y' must classify every endpoint selected in 'list_y_comparator'; missing: ",
           paste(missing, collapse = ", "))
    aligned <- rep(1L, length(all_endpoints))
    names(aligned) <- all_endpoints
    aligned[names(type_y)] <- as.integer(type_y)
  } else {
    if (length(type_y) != length(all_endpoints))
      stop("Unnamed 'type_y' must have one value for every available endpoint; use names to classify a selected subset.")
    aligned <- stats::setNames(as.integer(type_y), all_endpoints)
  }
  list(type_y = aligned, supplied = TRUE, active = TRUE)
}

# Return the hierarchy values in the order used by one comparator.
.type_y_for_endpoints <- function(type_y, endpoints, adjust) {
  if (!identical(adjust, "seq") || is.null(type_y) || !length(type_y) ||
      anyNA(type_y) || any(type_y < 0))
    return(rep(-1L, length(endpoints)))
  if (is.null(names(type_y))) return(as.integer(type_y[seq_along(endpoints)]))
  values <- type_y[as.character(endpoints)]
  if (anyNA(values))
    stop("'type_y' does not classify every selected endpoint.")
  as.integer(values)
}

# Sequential testing uses the primary family at alpha and allocates the
# secondary-family alpha across its endpoints. The calculation is per
# comparator because list_y_comparator may define different endpoint families.
.sequential_endpoint_alpha <- function(alpha, type_y, k) {
  if (length(type_y) < 1L || any(!type_y %in% c(1L, 2L)))
    stop("Sequential testing requires endpoint types 1 or 2.")
  if (length(alpha) == 1L) alpha <- rep(alpha, length(type_y))
  if (length(alpha) != length(type_y) || any(!is.finite(alpha)) ||
      any(alpha <= 0) || any(alpha >= 1))
    stop("'alpha' must be a scalar or contain one value between 0 and 1 for each endpoint.")
  out <- as.numeric(alpha)
  secondary <- type_y == 2L
  if (any(secondary))
    out[secondary] <- alpha[secondary] * min(k, length(type_y)) / sum(secondary)
  out
}

# Count kernels use the same user-facing adjustment aliases as the continuous
# interface, but retain their historical long labels in returned objects.
.normalize_count_adjustment <- function(adjust) {
  normalized <- .normalize_adjustment(adjust)
  if (normalized == "k")
    stop("Count outcomes support 'no', 'bon', 'sid', 't', and 'seq' adjustment; the k-adjustment is not implemented for count outcomes.")
  switch(normalized,
         no = "none", bon = "bonferroni", sid = "sidak", t = "t",
         seq = "sequential")
}

.count_endpoint_alpha <- function(alpha, m, k, adjust, type_y = NULL,
                                  type_y_active = FALSE) {
  if (length(alpha) == 1L) alpha <- rep(alpha, m)
  if (length(alpha) != m || any(!is.finite(alpha)) ||
      any(alpha <= 0) || any(alpha >= 1))
    stop("'alpha' must be a scalar or contain one value between 0 and 1 for each endpoint.")
  if (identical(adjust, "sequential") && isTRUE(type_y_active))
    return(.sequential_endpoint_alpha(alpha, type_y, k))
  switch(adjust,
         none = alpha,
         bonferroni = alpha / m,
         sidak = 1 - (1 - alpha)^(1 / m),
         t = alpha / (m - k + 1),
         sequential = alpha,
         stop("Unknown count multiplicity adjustment: ", adjust))
}

# Warn about adjustment choices that are redundant or do not by themselves
# provide the intended error-rate protection.
.warn_adjustment_configuration <- function(k, m, adjust,
                                           type_y = NULL,
                                           type_y_supplied = FALSE,
                                           n_comparators = 1L,
                                           independent = NULL,
                                           context = "selected endpoints") {
  adjust <- .normalize_adjustment(adjust)
  messages <- character()
  all_required <- k == m
  if (any(all_required) && adjust != "no") {
    if (adjust == "k") {
      messages <- c(messages,
        paste0("The k-adjustment has no effect when 'k' equals the number of ",
               context, " (k = m); each endpoint remains at alpha."))
    } else if (adjust == "t") {
      messages <- c(messages,
        paste0("The t-adjustment has no effect when ",
               "'k' equals the number of ", context,
               " (k = m); each endpoint remains at alpha."))
    } else {
      reason <- if (isTRUE(independent))
        "All selected endpoints are independent. However, " else ""
      messages <- c(messages,
               paste0(reason, "formal endpoint-wise Type I-error adjustment is not necessary ",
               "when all selected endpoints are required (k = m). The requested ",
               "adjustment will nevertheless be applied. This conclusion does not ",
               "depend on endpoint independence."))
    }
  }
  if (adjust == "no" && any(k < m)) {
    messages <- c(messages,
      paste0("adjust = 'no' does not provide multiplicity control for a k-of-m ",
             "decision when k < m; interpret the resulting Type I error as ",
             "uncalibrated unless this is intentional."))
  }
  if (adjust == "seq" && type_y_supplied && !is.null(type_y) &&
      all(type_y == 1L)) {
    messages <- c(messages,
      "All selected endpoints are classified as primary; sequential adjustment has no secondary family to gate or reallocate.")
  }
  if (adjust == "seq" && type_y_supplied && !is.null(type_y) &&
      all(type_y == 2L)) {
    messages <- c(messages,
      "All selected endpoints are classified as secondary; sequential adjustment has no primary gate, so review the prespecified hierarchy and alpha allocation.")
  }
  if (n_comparators > 1L && adjust != "no") {
    messages <- c(messages,
      "The adjustment is applied within each comparator's selected endpoint family, not across comparators; the global decision requires all comparator decisions to pass.")
  }
  if (length(messages))
    warning(paste(unique(messages), collapse = " "), call. = FALSE)
  invisible(NULL)
}

# Backward-compatible wrapper used by the count helpers.
.warn_redundant_bonferroni <- function(k, m, adjust, context = "selected endpoints") {
  if (is.null(adjust)) adjust <- "no"
  .warn_adjustment_configuration(k = k, m = m,
                                 adjust = .normalize_adjustment(adjust),
                                 context = context)
}

# Warn when omitted endpoint selections make the tested family smaller than
# the endpoint set stored for one or more arms. This is intentional behavior:
# m is the number of endpoints actually available in both arms (or explicitly
# selected), not the largest number supplied anywhere in the input.
.warn_inferred_endpoint_reduction <- function(comparators, arm_endpoints,
                                              endpoint_sets,
                                              requested = NULL,
                                              context = "list_y_comparator") {
  omitted <- is.null(requested) ||
    (length(requested) == 1L && all(is.na(requested)))
  if (!omitted) return(invisible(NULL))
  reduced <- vapply(seq_along(comparators), function(i) {
    arms <- comparators[[i]]
    any(length(endpoint_sets[[i]]) <
          c(length(arm_endpoints[[arms[[1L]]]]),
            length(arm_endpoints[[arms[[2L]]]])))
  }, logical(1))
  if (any(reduced)) {
    m_text <- paste(lengths(endpoint_sets), collapse = ", ")
    warning("'", context, "' was omitted, so only endpoints common to both ",
            "arms are tested; the effective endpoint count m is ", m_text,
            " by comparator.", call. = FALSE)
  }
  invisible(NULL)
}

# Resolve the endpoint set for every comparator. A supplied list is retained
# per comparator; otherwise all endpoints common to its two arms are used.
.resolve_comparator_endpoints <- function(comparators, arm_endpoints,
                                          requested = NULL) {
  if (!is.list(comparators) || !length(comparators))
    stop("'comparators' must be a non-empty list.")
  comparison_names <- names(comparators)
  if (is.null(comparison_names))
    comparison_names <- paste0("comparison_", seq_along(comparators))
  missing_requested <- is.null(requested) ||
    (length(requested) == 1L && all(is.na(requested)))
  if (missing_requested) {
    requested <- lapply(comparators, function(arms)
      intersect(arm_endpoints[[arms[[1L]]]], arm_endpoints[[arms[[2L]]]]))
  } else {
    if (!is.list(requested) || length(requested) != length(comparators))
      stop("'list_y_comparator' must contain one endpoint vector per comparator.")
    if (!is.null(names(requested))) {
      if (!setequal(names(requested), comparison_names))
        stop("Names of 'list_y_comparator' must match 'list_comparator'.")
      requested <- requested[comparison_names]
    }
    requested <- lapply(seq_along(comparators), function(i) {
      value <- requested[[i]]
      if (!is.character(value) || !length(value) || anyNA(value) ||
          anyDuplicated(value))
        stop("Each element of 'list_y_comparator' must be a non-empty, unique character vector.")
      common <- intersect(arm_endpoints[[comparators[[i]][[1L]]]],
                          arm_endpoints[[comparators[[i]][[2L]]]])
      unknown <- setdiff(value, common)
      if (length(unknown)) {
        stop("Endpoint(s) not available in both arms of comparator ",
             comparison_names[[i]], ": ", paste(unknown, collapse = ", "))
      }
      value
    })
  }
  names(requested) <- comparison_names
  if (any(lengths(requested) < 1L))
    stop("Every comparator must have at least one common endpoint.")
  requested
}

.prepare_count_endpoint_subset <- function(rates, comparators,
                                           list_y_comparator = NULL,
                                           exposure = NULL, dispersion = NULL,
                                           cor_mat = NULL,
                                           list_lequi.tol = NULL,
                                           list_uequi.tol = NULL) {
  arm_endpoints <- lapply(rates, function(value) {
    value_names <- names(value)
    if (!is.null(value_names)) value_names else
      paste0("endpoint_", seq_along(value))
  })
  requested_list_y_comparator <- list_y_comparator
  endpoint_sets <- .resolve_comparator_endpoints(
    comparators = comparators, arm_endpoints = arm_endpoints,
    requested = list_y_comparator
  )
  .warn_inferred_endpoint_reduction(
    comparators = comparators, arm_endpoints = arm_endpoints,
    endpoint_sets = endpoint_sets, requested = requested_list_y_comparator,
    context = "list_y_comparator"
  )
  if (length(endpoint_sets) > 1L && any(vapply(
      endpoint_sets[-1L], function(x) !identical(x, endpoint_sets[[1L]]),
      logical(1)))) {
    stop("For joint count outcomes, all comparators must use the same endpoints in 'list_y_comparator'.")
  }
  endpoints <- endpoint_sets[[1L]]
  all_endpoints <- arm_endpoints[[1L]]
  indices <- match(endpoints, all_endpoints)
  subset_value <- function(value) {
    if (is.null(value)) return(value)
    if (is.list(value)) return(lapply(value, subset_value))
    if (length(value) == length(all_endpoints) && length(value) > 1L)
      return(if (!is.null(names(value))) value[endpoints] else value[indices])
    value
  }
  list(
    endpoints = endpoints,
    rates = lapply(rates, function(value) {
      if (!is.null(names(value))) value[endpoints] else value[indices]
    }),
    exposure = subset_value(exposure), dispersion = subset_value(dispersion),
    cor_mat = if (is.matrix(cor_mat) && all(dim(cor_mat) == length(all_endpoints)))
      cor_mat[indices, indices, drop = FALSE] else cor_mat,
    list_lequi.tol = subset_value(list_lequi.tol),
    list_uequi.tol = subset_value(list_uequi.tol)
  )
}

# Validate and cap the number of count endpoints required for equivalence.
# Returning NULL preserves the public default, which lets the downstream
# kernel set k to the number of endpoints after its own input validation.
.normalize_count_k <- function(k, m, allow_vector = FALSE) {
  if (length(m) != 1L || !is.numeric(m) || !is.finite(m) ||
      m < 1L || m != as.integer(m))
    stop("'m' must be a positive integer.")
  if (is.null(k) || (length(k) == 1L && is.na(k)))
    return(NULL)
  if (!is.numeric(k) || !length(k) || anyNA(k) || any(!is.finite(k)) ||
      any(k != as.integer(k)) || any(k < 1L) ||
      (!allow_vector && length(k) != 1L))
    stop("'k' must contain positive integers.")
  oversized <- which(k > m)
  if (length(oversized)) {
    warning("'k' is larger than the number of selected count endpoints; ",
            "setting it to the maximum possible value.", call. = FALSE)
    k[oversized] <- m
  }
  k
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
  param.d <- .normalize_internal_distribution(param.d)
  if (is.na(ncores)) {
    ncores <- parallel::detectCores() - 1
  }
  treat1 <- param$list_comparator[[comp]][[1]]
  treat2 <- param$list_comparator[[comp]][[2]]
  endp <- param$list_y_comparator[[comp]]

  m <- length(endp) # number of endpoints

  # Set equivalence tolerance
  lequi.tol <- param.d$list_lequi.tol[[comp]][endp]
  uequi.tol <- param.d$list_uequi.tol[[comp]][endp]

  dropout <- param.d$dropout
  alphau <- param.d$alpha # alpha unique @Thomas, here we can divide by the number of comparators in case you want a bonberroni across comparators.Think about this!!
  adjust <- param.d$adjust
  k <- param.d$k[[comp]]

  # alpha vector
  type_y_comp <- .type_y_for_endpoints(param$type_y, endp, adjust)
  alpha <- if (adjust == "seq" && all(type_y_comp %in% c(1L, 2L))) {
    .sequential_endpoint_alpha(alphau, type_y_comp, k)
  } else {
    .endpoint_alpha(alphau, m, k, adjust)
  }

  if (param.d$dtype == "parallel") {
    # Derive mu
    muT <- param$mu[[treat1]][,endp] # treatment given by user
    muR <- param$mu[[treat2]][,endp] # reference given by user

    # Derive Sigma
    SigmaT <- param$varcov[[treat1]][endp,endp]
    SigmaR <- param$varcov[[treat2]][endp,endp]

    # Fail before entering compiled code if any endpoint-wise input has been
    # subset inconsistently. This prevents opaque Armadillo errors.
    input_lengths <- c(length(muT), length(muR), length(endp),
                       length(lequi.tol), length(uequi.tol), length(alpha))
    covariance_dims <- c(nrow(SigmaT), ncol(SigmaT), nrow(SigmaR), ncol(SigmaR))
    if (any(input_lengths != length(endp)) ||
        any(covariance_dims != length(endp))) {
      stop(sprintf(
        "Inconsistent endpoint dimensions for %s vs %s: endpoints=%d, muT=%d, muR=%d, lower=%d, upper=%d, alpha=%d, covariance=%s.",
        treat1, treat2, length(endp), length(muT), length(muR),
        length(lequi.tol), length(uequi.tol), length(alpha),
        paste(covariance_dims, collapse = "x")
      ))
    }

    if (param.d$ctype == "DOM" & param.d$distribution == "norm") {
      result <- run_simulations_par_dom(nsim = nsim, n = n, muT = muT, muR = muR,
                                        SigmaT = as.matrix(SigmaT),
                                        SigmaR = as.matrix(SigmaR),
                                        lequi_tol = lequi.tol, uequi_tol = uequi.tol,
                                        alpha = alpha,
                                        dropout = as.numeric(c(dropout[treat1], dropout[treat2])),
                                        typey = type_y_comp,
                                        adseq = param.d$adjust == "seq", k = k,
                                        arm_seed_T = arm_seed[,treat1],
                                        arm_seed_R = arm_seed[,treat2],
                                        TART = param$TAR_list[[treat1]],
                                        TARR = param$TAR_list[[treat2]],
                                        vareq = param.d$vareq)
    } else if (param.d$ctype == "ROM" & param.d$distribution == "lnorm") {
      # Convert data to lognorm scale so we can perform a DOM test instead
      SigmaT <-  as.matrix(log(SigmaT/(muT %*% t(muT)) + 1))
      SigmaR <-  as.matrix(log(SigmaR/(muR %*% t(muR)) + 1))
      muT <- log(muT) - 1/2*diag(SigmaT)
      muR <- log(muR) - 1/2*diag(SigmaR)

      result <- run_simulations_par_dom(nsim = nsim, n = n, muT = muT, muR = muR,
                                        SigmaT = as.matrix(SigmaT),
                                        SigmaR = as.matrix(SigmaR),
                                        lequi_tol = log(lequi.tol),
                                        uequi_tol = log(uequi.tol),
                                         alpha = alpha,
                                         dropout = as.numeric(c(dropout[treat1], dropout[treat2])),
                                         typey = type_y_comp,
                                         adseq = param.d$adjust == "seq", k = k,
                                         arm_seed_T = arm_seed[,treat1],
                                         arm_seed_R = arm_seed[,treat2],
                                         TART = param$TAR_list[[treat1]],
                                         TARR = param$TAR_list[[treat2]],
                                         vareq = param.d$vareq)
    } else if (param.d$ctype == "ROM" & param.d$distribution == "norm") {
      result <- run_simulations_par_rom(nsim = nsim, n = n, muT = muT, muR = muR,
                                        SigmaT = as.matrix(SigmaT),
                                        SigmaR = as.matrix(SigmaR),
                                        lequi_tol = lequi.tol, uequi_tol = uequi.tol,
                                        alpha = alpha,
                                        dropout = as.numeric(c(dropout[treat1], dropout[treat2])),
                                        typey = type_y_comp,
                                        adseq = param.d$adjust == "seq", k = k,
                                        arm_seed_T = arm_seed[,treat1],
                                        arm_seed_R = arm_seed[,treat2],
                                        TART = param$TAR_list[[treat1]],
                                        TARR = param$TAR_list[[treat2]],
                                        vareq = param.d$vareq)
    } else {
      stop(paste("Error: Unsupported test type:", param.d$ctype,
                 "with distribution =", param.d$distribution))
    }
  } else if (param.d$dtype == "2x2") {
    # Derive mu
    muT <- param$mu[[treat1]][,endp] # treatment given by user
    muR <- param$mu[[treat2]][,endp] # reference given by user


    SigmaW <- param$varcov[[treat2]][endp,endp] # Within subjects variance in previous experiment
    #To be added on main list of parameters
    sigmaB <- param$sigmaB # Between subjects variance
    sigmaB <- ifelse(is.na(sigmaB), if (length(SigmaW)==1) 2* sqrt(SigmaW) else 2 * sqrt(max(diag(SigmaW))), sigmaB) # Assumes to be at least the double of the max within variance

    if (param.d$ctype == "DOM" & param.d$distribution == "norm") {
      result <- run_simulations_2x2_dom(nsim = nsim,
                                        n = n, muT = muT, muR = muR,
                                        SigmaW = as.matrix(SigmaW),
                                        lequi_tol = lequi.tol, uequi_tol = uequi.tol,
                                        alpha = alpha, sigmaB = sigmaB,
                                        dropout = dropout,
                                        Eper = param$Eper, Eco = param$Eco,
                                        typey = type_y_comp,
                                        adseq = param.d$adjust == "seq", k = k,
                                        arm_seed = arm_seed[,comp])
    } else if (param.d$ctype == "ROM" & param.d$distribution == "lnorm"){
      # Convert data to lognorm scale
      SigmaW <-  as.matrix(log(SigmaW/(muR%*%t(muR))+1))
      sigmaB <-  sigmaB #log(sigmaB/(muR%*%t(muR))+1)
      muR <- log(muR)-1/2*diag(SigmaW)
      muT <- log(muT)-1/2*diag(SigmaW)

      result <- run_simulations_2x2_dom(nsim = nsim,
                                        n = n, muT = muT, muR = muR,
                                        SigmaW = as.matrix(SigmaW),
                                        lequi_tol = log(lequi.tol),
                                        uequi_tol = log(uequi.tol),
                                        alpha = alpha, sigmaB = sigmaB,
                                        dropout = dropout,
                                        Eper = param$Eper, Eco = param$Eco,
                                        typey = type_y_comp,
                                        adseq = param.d$adjust == "seq", k = k,
                                        arm_seed = arm_seed[,comp])
    } else if (param.d$ctype == "ROM" & param.d$distribution == "norm") {
      result <- run_simulations_2x2_rom(nsim = nsim,
                                        n = n, muT = muT, muR = muR,
                                        SigmaW = as.matrix(SigmaW),
                                        lequi_tol = lequi.tol,
                                        uequi_tol = uequi.tol,
                                        alpha = alpha, sigmaB = sigmaB,
                                        dropout = dropout,
                                        Eper = param$Eper, Eco = param$Eco,
                                        typey = type_y_comp,
                                        adseq = param.d$adjust == "seq", k = k,
                                        arm_seed = arm_seed[,comp])
    } else {
    stop(paste("Error: Unsupported test type:", param.d$ctype,
               "with distribution =", param.d$distribution))
    }
  } else {
    stop(paste("Error: Unsupported design:", param.d$dtype))
  }

  rownames(result) <- paste0(c("totaly", endp,
                                 paste0("mu_",endp,"_",treat1),
                                 paste0("mu_",endp,"_",treat2),
                                 paste0("sd_",endp,"_",treat1),
                                 paste0("sd_",endp,"_",treat1)),"Comp:",treat1," vs ",treat2)
  return(result)
}



#' @title Optimizer for Uniroot Integer (Modified)
#'
#' @description A modified integer-based root-finding algorithm for determining the sample size required to achieve a target power.
#' This function extends the uniroot integer search method to handle cases with stepwise power searches while considering constraints on search limits.
#'
#' @param f Function for which a root is needed.
#' @param power Numeric. Target power value.
#' @param lower Integer. Minimum allowable root value.
#' @param upper Integer. Maximum allowable root value.
#' @param step.power Numeric. Initial step size defined as \code{2^step.power}.
#' @param step.up Logical. If \code{TRUE}, the search increments from \code{lower}; if \code{FALSE}, it decrements from \code{upper}.
#' @param pos.side Logical. If \code{TRUE}, finds the closest integer \code{i} such that \code{f(i) > 0}.
#' @param maxiter Integer. Maximum number of iterations allowed.
#' @param ... Additional arguments passed to \code{f}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{root}}{The integer value closest to the root on the correct side.}
#'   \item{\code{f.root}}{Value of \code{f} at the estimated root.}
#'   \item{\code{iter}}{Number of function evaluations performed.}
#'   \item{\code{table.iter}}{A data frame showing estimated sample size (\code{N}) and corresponding power at each iteration.}
#'   \item{\code{table.test}}{A data frame containing endpoint-level test results for each simulation and corresponding \code{N}.}
#' }
#'
#' @keywords internal
uniroot.integer.mod <-function (f, power, lower = lower, upper = upper, step.power=step.power, step.up=step.up, pos.side=pos.side, maxiter = maxiter,...) {
  # Function adapted from ssanv https://github.com/cran/ssanv/blob/master/R/uniroot.integer.R
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

#' Helper function for conditional messages
#'
#' This function displays a message if the `verbose` parameter is set to `TRUE`.
#' It is useful for providing optional feedback to users during function execution.
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @param message A character string containing the message to display.
#' @param verbose Logical, if `TRUE`, the message is displayed; if `FALSE`, the message is suppressed.
#'
#' @return NULL (invisible). This function is used for side effects (displaying messages).
#' @keywords internal
info_msg <- function(message, verbose) {
  if (verbose) message(message)
}
