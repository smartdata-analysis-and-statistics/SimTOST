#' @title Sample Size Calculation for Bioequivalence and Multi-Endpoint Studies
#'
#' @description Computes the required sample size to achieve a target power in studies with multiple endpoints and treatment arms.
#' The function employs modified root-finding algorithms to estimate sample size while accounting for correlation structures, variance assumptions,
#' and equivalence bounds across endpoints. It is particularly useful for bioequivalence trials and multi-arm studies with complex endpoint structures.
#'
#' @param distribution Outcome distribution. Choose the R distribution names
#' `"norm"`, `"lnorm"`, `"pois"`, or `"nbinom"`. Matching is case-insensitive.
#' The longer labels such as `"normal"`, `"lognormal"`, and `"poisson"` are
#' accepted for compatibility; results store the R-style labels above.
#' @param mu_list Named list of arithmetic means per treatment arm. Each element is a vector representing expected outcomes for all endpoints in that arm.
#' @param varcov_list List of variance-covariance matrices, where each element corresponds to a comparator. Each matrix has dimensions: number of endpoints × number of endpoints.
#' @param sigma_list List of standard deviation vectors, where each element corresponds to a comparator and contains one standard deviation per endpoint.
#' @param cor_mat Matrix specifying the correlation structure between
#' endpoints. For continuous outcomes it is used with `sigma_list` to
#' calculate `varcov_list`; for count outcomes it is used by the joint count
#' engine.
#' @param sigmaB Numeric. Between-subject standard deviation parameter for the
#' continuous 2×2 design; for count outcomes it is the log-rate standard
#' deviation in the 2×2 kernel.
#' @param rate_list Named list of equal-length endpoint-rate vectors, one per
#' count-outcome arm.
#' @param exposure Exposure per subject for count outcomes. Supply a scalar
#' or endpoint vector shared by arms, or a named list of arm-specific values.
#' @param dispersion Positive negative-binomial dispersion parameter. The
#' per-subject negative-binomial size is `1 / dispersion`; parallel-arm
#' totals use size `n / dispersion`. Supply a scalar or endpoint vector shared
#' by arms, or a named list of arm-specific values.
#' @param Eper Optional numeric vector of length 2 specifying period effects.
#' For count outcomes these are log-rate effects applied to periods 1 and 2.
#' @param Eco Optional numeric vector of length 2 specifying carry-over effects
#' in the order reference carry-over and treatment carry-over. For count
#' outcomes these are log-rate effects in period 2.
#' @param rho Numeric. Correlation parameter applied uniformly across all endpoint pairs. Used with \code{sigma_list} to compute \code{varcov_list} when \code{cor_mat} or \code{varcov_list} are not provided.
#' @param TAR Numeric vector specifying treatment allocation rates per arm. The order must match \code{arm_names}. Defaults to equal allocation across arms if not provided.
#' @param arm_names Optional character vector of treatment names. If not supplied, names are derived from \code{mu_list}.
#' @param ynames_list Optional list of vectors specifying endpoint names per arm. If names are missing, arbitrary names are assigned based on order.
#' @param type_y Integer vector indicating endpoint types: \code{1} for co-primary endpoints, \code{2} for secondary endpoints.
#' @param list_comparator List of comparators. Each element must be a vector
#'   of length 2 in the form `c(test, reference)`. The first arm is treated as
#'   the treatment arm and the second arm as the reference arm throughout the
#'   simulation and estimand calculations.
#' @param list_y_comparator List of endpoint sets per comparator. Each element is a vector containing endpoint names to compare. If not provided, all endpoints common to both comparator arms are used.
#'   For count outcomes, the selected endpoints define the count multiplicity
#'   and effective `k`; joint count analyses require the same endpoint set for
#'   every comparison.
#' @param power Numeric. Target power (default = 0.8).
#' @param alpha Numeric. Significance level (default = 0.05).
#' @param lequi.tol Numeric. Lower equivalence bounds (e.g., -0.5) applied uniformly across all endpoints and comparators.
#' @param uequi.tol Numeric. Upper equivalence bounds (e.g., 0.5) applied uniformly across all endpoints and comparators.
#' @param list_lequi.tol List of numeric vectors specifying lower equivalence bounds per comparator.
#' @param list_uequi.tol List of numeric vectors specifying upper equivalence bounds per comparator.
#' @param vareq Logical. Assumes equal variances across arms if \code{TRUE} (default = \code{FALSE}).
#' @param dtype Character. Trial design: \code{"parallel"} (default) for parallel-group or \code{"2x2"} for crossover (only for 2-arm studies).
#' @param k Integer vector. Minimum number of successful endpoints required for global bioequivalence per comparator. Defaults to all endpoints per comparator.
#' @param adjust Character. Alpha adjustment method: \code{"k"} (K-fold),
#' \code{"bon"} (Bonferroni across all selected endpoints), \code{"sid"}
#' (Sidak), \code{"t"} (Mielke's strong \\(k\\)-out-of-\\(m\\) adjustment using
#' \code{alpha / (m - k + 1)}; legacy \code{"pc"} aliases are accepted),
#' \code{"no"} (none), or
#' \code{"seq"} (sequential).
#' @param ctype Character. Continuous-outcome test type: \code{"DOM"} (Difference of Means) or \code{"ROM"} (Ratio of Means). For Poisson and negative-binomial outcomes, \code{"RR"} (event-rate ratio) is used; an unavailable value triggers a warning and the applicable default is used.
#' @param dropout Numeric vector specifying dropout proportion per arm.
#' @param nsim Integer. Number of simulated studies (default = 5000).
#' @param keep_sim_data Logical. If `TRUE`, retain model-scale observations
#'   for every simulated trial in `sim_data` for distribution diagnostics.
#'   Defaults to `FALSE` because retained data can be large.
#' @param seed Integer. Seed for reproducibility.
#' @param ncores Integer. Number of processing cores for parallel computation. Defaults to \code{1}. Set to \code{NA} for automatic detection (\code{ncores - 1}). For count outcomes, Monte Carlo trials are split into independent seeded chunks and evaluated by the count C++ kernel on each worker.
#' @param optimization_method Character. Sample size optimization method: \code{"fast"} (default, root-finding algorithm) or \code{"step-by-step"}.
#' @param lower Integer. Minimum sample size for search range (default = 2).
#' @param upper Integer. Maximum sample size for the search range (default =
#'   500). For count outcomes, this is the maximum number of subjects per arm;
#'   the plotted and returned total sample size is this value multiplied by the
#'   number of trial arms.
#' @param step.power Numeric. Initial step size for sample size search, defined as \code{2^step.power}. Used when \code{optimization_method = "fast"}.
#' @param step.up Logical. If \code{TRUE} (default), search increments upward from \code{lower}; if \code{FALSE}, decrements downward from \code{upper}. Used when \code{optimization_method = "fast"}.
#' @param pos.side Logical. If \code{TRUE}, finds the smallest integer \code{i} closest to the root such that \code{f(i) > 0}. Used when \code{optimization_method = "fast"}.
#' @param maxiter Integer. Maximum iterations allowed for sample size estimation (default = 1000). Used when \code{optimization_method = "fast"}.
#' @param verbose Logical. If \code{TRUE}, prints progress and messages during execution (default = \code{FALSE}).
#' @param .warn_redundant_bon Logical. If `TRUE`, warn when a requested
#'   multiplicity adjustment is redundant or uncalibrated for the selected
#'   endpoint decision.
#'
#' @details
#' The common planning arguments are \code{power}, \code{alpha},
#' \code{list_comparator}, \code{list_lequi.tol}, \code{list_uequi.tol},
#' \code{k}, \code{adjust}, \code{dtype}, \code{dropout}, \code{nsim},
#' \code{seed}, \code{lower}, and \code{upper}. Use the following
#' distribution-specific arguments in addition to those common arguments:
#' \describe{
#'   \item{Normal and Log Normal}{Supply \code{mu_list}, \code{sigma_list}
#'   or \code{varcov_list}; use \code{cor_mat} or \code{rho} for endpoint
#'   dependence. The \code{ctype} argument selects DOM or ROM testing.}
#'   \item{Poisson and Negative Binomial}{Supply \code{rate_list},
#'   \code{list_comparator}, and comparator-specific \code{list_lequi.tol}
#'   and \code{list_uequi.tol}. Use \code{exposure} and, for negative-binomial
#'   outcomes, \code{dispersion}; both may be scalar, endpoint-specific, or
#'   named arm-specific lists. Continuous-outcome arguments are ignored.}
#' }
#' For count outcomes, `optimization_method = "fast"` brackets the first
#' sample size whose simulated power reaches the target and refines the
#' bracket by integer bisection. The `"step-by-step"` option remains available
#' when a complete candidate-by-candidate power table is preferred. The fast
#' method assumes the usual approximately monotone power curve and uses the
#' same seed at each candidate to reduce simulation noise.
#' The effective endpoint count is comparator-specific: when
#' `list_y_comparator` is omitted, only endpoints present in both arms are
#' tested; when it is supplied, only the listed endpoints are tested. `k` is
#' validated against that comparator-specific count and oversized values are
#' capped with a warning. Formal endpoint-wise adjustment is unnecessary when
#' all selected endpoints are required (`k = m`), although requested
#' Bonferroni or Sidak adjustment remains available with a warning. For
#' `k < m`, `adjust = "no"` is explicitly reported as an uncalibrated choice.
#' For a `k`-of-`m` decision, `adjust = "t"` applies Mielke's strong
#' \\(k\\)-out-of-\\(m\\) calibration `alpha / (m - k + 1)`. The legacy
#' `adjust = "pc"` label is accepted as an alias.
#' For continuous and count outcomes, `type_y` is used with
#' `adjust = "seq"`; named endpoint vectors are aligned to the selected
#' comparator endpoints. Count analyses use the same primary-gate and
#' secondary-family decision rule as the continuous kernels.
#' The unified function returns primary class \code{simss} for all outcome
#' distributions. Count results retain \code{countss} as a secondary
#' compatibility class. Use \code{summary()} and \code{plot()} to inspect the
#' result.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{response}}{Array summarizing simulation results, including estimated sample sizes, achieved power, and confidence intervals.}
#'   \item{\code{table.iter}}{Data frame showing estimated sample sizes and calculated power at each iteration. For count outcomes, one row is retained for every evaluated candidate.}
#'   \item{\code{table.test}}{Data frame containing test results for all simulated trials. For count outcomes, this contains complete-trial, comparator, and endpoint decision indicators for each simulated trial and candidate; the count kernel returns aggregate decision counts rather than raw endpoint-level test statistics.}
#'   \item{\code{param.u}}{Original input parameters.}
#'   \item{\code{param}}{Final adjusted parameters used in sample size calculation.}
#'   \item{\code{param.d}}{Trial design parameters used in the simulation.}
#'   \item{\code{sim_data}}{Optional long-format simulated observations,
#'   returned when \code{keep_sim_data = TRUE}.}
#' }
#'
#' @references
#' Schuirmann, D. J. (1987). A comparison of the Two One-Sided Tests procedure and the
#' Power approach for assessing the equivalence of average bioavailability.
#' \emph{Journal of Pharmacokinetics and Biopharmaceutics, 15}(6), 657-680.
#' <doi:10.1007/BF01068419>
#'
#' Mielke, J., Jones, B., Jilma, B., & König, F. (2018). Sample size for multiple hypothesis
#' testing in biosimilar development. \emph{Statistics in Biopharmaceutical Research, 10}(1), 39-49.
#' <doi:10.1080/19466315.2017.1371071>
#'
#' Berger, R. L., & Hsu, J. C. (1996). Bioequivalence trials, intersection-union tests, and
#' equivalence confidence sets. \emph{Statistical Science}, 283-302.
#'
#' Sozu, T., Sugimoto, T., Hamasaki, T., & Evans, S. R. (2015). "Sample Size Determination in
#' Clinical Trials with Multiple Endpoints." \emph{SpringerBriefs in Statistics}.
#' <doi:10.1007/978-3-319-22005-5>
#'
#' @export
#'
#' @examples
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
#'            distribution = "lnorm", ncores = 1, nsim = 50, seed = 1234)
#'
#' # The same entry point for a two-arm Poisson count-rate calculation:
#' sampleSize(power = 0.80, distribution = "Poisson",
#'            rate_list = list(TEST = 0.21, REF = 0.20),
#'            list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
#'            list_lequi.tol = list(TEST_vs_REF = 0.80),
#'            list_uequi.tol = list(TEST_vs_REF = 1.25),
#'            exposure = 10, lower = 20, upper = 500,
#'            nsim = 100, seed = 1234)
sampleSize <- function(distribution = c("norm", "lnorm", "pois", "nbinom"),
                       mu_list = NULL, varcov_list = NA, sigma_list = NA, cor_mat = NA,
                       sigmaB = NA, rate_list = NULL, exposure = 1,
                       dispersion = 0.1, Eper = c(0, 0), Eco = c(0, 0), rho = 0,
                       TAR = rep(1, length(mu_list)), arm_names = NA,
                       ynames_list = NA,
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
                    k=NA,
                    adjust="no",
                    dropout = NA,
                    nsim=5000,
                    seed=1234,
                    ncores = 1,
                    optimization_method = "fast",
                    lower=2,
                    upper=500,
                    step.power=6,
                    step.up=TRUE,
                    pos.side=FALSE,
                    maxiter = 1000, verbose = FALSE, keep_sim_data = FALSE,
                    .warn_redundant_bon = TRUE
){

  keep_sim_data <- .check_keep_sim_data(keep_sim_data)
  .validate_common_controls(power, alpha, nsim, ncores, dropout, rho)

  if (missing(distribution) || is.null(distribution)) {
    # Preserve the historical sampleSize() default while making distribution
    # the only public outcome-family selector.
    distribution <- "lnorm"
  } else {
    distribution <- .normalize_distribution(distribution)
  }
  if (distribution %in% c("pois", "nbinom")) {
    outcome <- "count"
    model <- if (distribution == "pois") "poisson" else "negative-binomial"
  } else {
    outcome <- "continuous"
  }
  adjust <- .normalize_adjustment(adjust, allow_sequential = TRUE)
  type_y_supplied <- !(is.null(type_y) ||
                       (length(type_y) == 1L && is.na(type_y)))
  if (type_y_supplied && is.numeric(type_y) &&
      all(!is.na(type_y) & type_y < 0))
    type_y_supplied <- FALSE

  if (outcome == "count") {
    if (missing(ctype)) {
      ctype <- "RR"
    } else if (length(ctype) != 1L || !is.character(ctype) ||
               is.na(ctype) || !identical(toupper(ctype), "RR")) {
      warning("For count outcomes, ctype = 'RR' is used. The supplied ctype is not available for count outcomes.")
      ctype <- "RR"
    } else {
      ctype <- "RR"
    }
  } else {
        if (length(ctype) != 1L || !is.character(ctype) ||
            is.na(ctype) || !toupper(ctype) %in% c("DOM", "ROM")) {
      warning("For continuous outcomes, ctype must be 'DOM' or 'ROM'. Default ctype = 'ROM' is used.")
      ctype <- "ROM"
    } else {
      ctype <- toupper(ctype)
    }
  }

  # The count extension uses the same top-level planning arguments as the
  # continuous interface (power, alpha, k, adjust, nsim, seed, lower, upper).
  # Scalar/vector inputs use sampleSize_count(); named multi-arm inputs use
  # sampleSize_count_joint(). Keeping this dispatch here preserves all
  # existing sampleSize() calls while providing one public entry point.
  if (outcome == "count") {
    count_ncores <- if (length(ncores) == 1L && is.na(ncores)) {
      max(1L, parallel::detectCores() - 1L)
    } else ncores
    count_design <- match.arg(dtype, c("parallel", "2x2"))
    count_endpoint_corr <- if (length(cor_mat) == 1L && is.na(cor_mat)) {
      NULL
    } else {
      cor_mat
    }
    count_sigmaB <- if (missing(sigmaB)) 0 else sigmaB
    count_Eper <- if (missing(Eper)) c(0, 0) else Eper
    count_Eco <- if (missing(Eco)) c(0, 0) else Eco
    count_dropout <- if (all(is.na(dropout))) c(0, 0) else dropout
    count_k <- if (length(k) == 1L && is.na(k)) NULL else k
    count_adjust <- .normalize_count_adjustment(adjust)

    if (is.null(rate_list) || !is.list(rate_list))
      stop("For count outcomes, supply a named 'rate_list' with one rate vector per arm.")
    if (is.null(names(rate_list)) || any(!nzchar(names(rate_list))))
      stop("'rate_list' must be named so arm-specific parameters can be matched.")
    check_arm_parameter <- function(x, name) {
      if (is.list(x) && (is.null(names(x)) ||
                         !setequal(names(x), names(rate_list))))
        stop(sprintf("Names of arm-specific '%s' must match 'rate_list'.", name))
    }
    check_arm_parameter(exposure, "exposure")
    check_arm_parameter(dispersion, "dispersion")
    if (length(list_comparator) == 1L && is.na(list_comparator))
      stop("For count outcomes, supply 'list_comparator'.")
    if (!is.list(list_lequi.tol) || !is.list(list_uequi.tol))
      stop("For count outcomes, supply comparator-specific 'list_lequi.tol' and 'list_uequi.tol'.")
    comparisons <- list_comparator
    comparison_names <- names(comparisons)
    if (is.null(comparison_names)) {
      comparison_names <- paste0("comparison_", seq_along(comparisons))
      names(comparisons) <- comparison_names
    }
    if (is.null(names(list_lequi.tol)) || is.null(names(list_uequi.tol)) ||
        !setequal(names(list_lequi.tol), comparison_names) ||
        !setequal(names(list_uequi.tol), comparison_names))
      stop("Names of 'list_lequi.tol' and 'list_uequi.tol' must match 'list_comparator'.")

    requested_count_endpoints <- list_y_comparator
    count_all_endpoints <- names(rate_list[[1L]])
    if (is.null(count_all_endpoints))
      count_all_endpoints <- paste0("endpoint_", seq_along(rate_list[[1L]]))
    count_subset <- .prepare_count_endpoint_subset(
      rates = rate_list, comparators = comparisons,
      list_y_comparator = list_y_comparator, exposure = exposure,
      dispersion = dispersion, cor_mat = count_endpoint_corr,
      list_lequi.tol = list_lequi.tol, list_uequi.tol = list_uequi.tol
    )
    rate_list <- count_subset$rates
    exposure <- count_subset$exposure
    dispersion <- count_subset$dispersion
    count_endpoint_corr <- count_subset$cor_mat
    list_lequi.tol <- count_subset$list_lequi.tol
    list_uequi.tol <- count_subset$list_uequi.tol
    count_m <- length(count_subset$endpoints)
    count_type_info <- .prepare_type_y(
      type_y = type_y, all_endpoints = count_all_endpoints,
      selected_endpoints = count_subset$endpoints,
      adjust = if (adjust == "seq") "seq" else "no"
    )
    count_type_y <- if (count_type_info$active)
      count_type_info$type_y[count_subset$endpoints] else NULL
    count_k <- .normalize_count_k(count_k, count_m, allow_vector = TRUE)
    if (!is.null(count_k) && length(count_k) > 1L) {
      if (length(count_k) != length(comparisons) ||
          length(unique(count_k)) != 1L)
        stop("Count joint outcomes require one common 'k' across comparisons.")
      count_k <- count_k[[1L]]
    }

    # A single two-arm comparison uses the count kernel that also supports
    # the 2x2 design. Multiple arms/comparisons use the joint count kernel.
    if (length(comparisons) == 1L && length(rate_list) == 2L) {
      arms <- comparisons[[1L]]
      result <- sampleSize_count(
        power = power, rate_test = rate_list[[arms[1L]]],
        rate_reference = rate_list[[arms[2L]]], exposure = exposure,
        margin_lower = list_lequi.tol[[comparison_names[1L]]],
        margin_upper = list_uequi.tol[[comparison_names[1L]]],
        model = model, dispersion = dispersion, alpha = alpha, nsim = nsim,
        seed = seed, lower = lower, upper = upper, design = count_design,
        k = count_k, endpoint_corr = count_endpoint_corr,
        type_y = count_type_y, adjust = count_adjust, sigmaB = count_sigmaB,
        Eper = count_Eper, Eco = count_Eco, dropout = count_dropout,
        optimization_method = optimization_method, step.power = step.power,
        step.up = step.up, pos.side = pos.side, maxiter = maxiter,
        ncores = count_ncores,
        .warn_redundant_bon = .warn_redundant_bon
      )
      if (is.data.frame(result$table.test) &&
          length(comparison_names) == 1L) {
        names(result$table.test) <- gsub(
          "Comp:comparison_1",
          paste0("Comp:", comparison_names[[1L]]),
          names(result$table.test), fixed = TRUE
        )
      }
      if (keep_sim_data) {
        result$sim_data <- .retained_count_data(
          rate_list = rate_list, exposure = exposure, dispersion = dispersion,
          distribution = distribution, design = count_design,
          n = result$n_per_arm, nsim = nsim, seed = seed,
          endpoint_corr = count_endpoint_corr
        )
      }
      result$endpoint_corr <- count_endpoint_corr
      class(result) <- c("simss", class(result))
      result <- .decorate_count_sample_size(
        result = result, distribution = distribution, rate_list = rate_list,
        exposure = exposure, dispersion = dispersion,
        comparisons = comparisons,
        selected_endpoints = setNames(
          rep(list(count_subset$endpoints), length(comparisons)),
          comparison_names
        ),
        requested_endpoints = requested_count_endpoints,
        lower_list = list_lequi.tol, upper_list = list_uequi.tol,
        type_y = count_type_y, endpoint_corr = count_endpoint_corr,
        design = count_design, power = power, alpha = alpha, k = result$k,
        adjust = count_adjust, dropout = count_dropout,
        sigmaB = count_sigmaB, Eper = count_Eper, Eco = count_Eco,
        nsim = nsim, seed = seed, ncores = count_ncores,
        optimization_method = optimization_method, lower = lower, upper = upper,
        step.power = step.power, step.up = step.up, pos.side = pos.side,
        maxiter = maxiter, verbose = verbose, keep_sim_data = keep_sim_data
      )
      return(result)
    }

    result <- sampleSize_count_joint(
      power = power, rates = rate_list, comparisons = comparisons,
      exposure = exposure, margin_lower = NULL, margin_upper = NULL,
      model = model, dispersion = dispersion,
      alpha = alpha, endpoint_corr = count_endpoint_corr, k = count_k,
      type_y = count_type_y, adjust = count_adjust,
      nsim = nsim, seed = seed, lower = lower, upper = upper,
      design = count_design,
      list_margin_lower = if (is.list(list_lequi.tol)) list_lequi.tol else NULL,
      list_margin_upper = if (is.list(list_uequi.tol)) list_uequi.tol else NULL,
      optimization_method = optimization_method, step.power = step.power,
      step.up = step.up, pos.side = pos.side, maxiter = maxiter,
      ncores = count_ncores,
      .warn_redundant_bon = .warn_redundant_bon
    )
    if (keep_sim_data && !is.null(result$n_per_arm)) {
      result$sim_data <- .retained_count_data(
        rate_list = rate_list, exposure = exposure, dispersion = dispersion,
        distribution = distribution, design = count_design,
        n = result$n_per_arm, nsim = nsim, seed = seed,
        endpoint_corr = count_endpoint_corr
      )
    }
    result$endpoint_corr <- count_endpoint_corr
    class(result) <- c("simss", class(result))
    result <- .decorate_count_sample_size(
      result = result, distribution = distribution, rate_list = rate_list,
      exposure = exposure, dispersion = dispersion,
      comparisons = comparisons,
      selected_endpoints = setNames(
        rep(list(count_subset$endpoints), length(comparisons)),
        comparison_names
      ),
      requested_endpoints = requested_count_endpoints,
      lower_list = list_lequi.tol, upper_list = list_uequi.tol,
      type_y = count_type_y, endpoint_corr = count_endpoint_corr,
      design = count_design, power = power, alpha = alpha, k = result$k,
      adjust = count_adjust, dropout = count_dropout,
      sigmaB = count_sigmaB, Eper = count_Eper, Eco = count_Eco,
      nsim = nsim, seed = seed, ncores = count_ncores,
      optimization_method = optimization_method, lower = lower, upper = upper,
      step.power = step.power, step.up = step.up, pos.side = pos.side,
      maxiter = maxiter, verbose = verbose, keep_sim_data = keep_sim_data
    )
    return(result)
  }

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

  # Resolve hierarchy after the comparator-specific endpoint families are
  # known. This makes named type_y robust when comparisons select different
  # endpoints from the arm-level input.
  selected_endpoints <- if (is.list(list_y_comparator)) {
    unique(unlist(list_y_comparator, use.names = FALSE))
  } else {
    character()
  }
  type_info <- .prepare_type_y(
    type_y = type_y, all_endpoints = uynames,
    selected_endpoints = selected_endpoints, adjust = adjust
  )
  type_y <- type_info$type_y
  weight_seq <- stats::setNames(rep(1, length(uynames)), uynames)

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

  # Resolve the endpoints actually tested for each comparator. This is the
  # source of m_i for multiplicity adjustment and k defaults.
  arm_endpoints <- lapply(mu_list, colnames)
  requested_list_y_comparator <- list_y_comparator
  list_y_comparator <- .resolve_comparator_endpoints(
    comparators = list_comparator, arm_endpoints = arm_endpoints,
    requested = list_y_comparator
  )
  .warn_inferred_endpoint_reduction(
    comparators = list_comparator, arm_endpoints = arm_endpoints,
    endpoint_sets = list_y_comparator, requested = requested_list_y_comparator,
    context = "list_y_comparator"
  )

  validate_sample_size_inputs(mu_list = mu_list, sigma_list = sigma_list,
                              varcov_list = varcov_list, arm_names = arm_names,
                              list_comparator = list_comparator,
                              list_y_comparator = list_y_comparator)


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
    selected_endpoints_i <- list_y_comparator[[i]]
    available_endpoints_i <- colnames(muend)
    align_bounds <- function(value) {
      if (is.null(value)) return(value)
      if (!is.null(names(value)) && all(selected_endpoints_i %in% names(value)))
        return(value[selected_endpoints_i])
      if (is.null(names(value)) && length(value) == length(available_endpoints_i) &&
          length(value) != length(selected_endpoints_i))
        return(value[match(selected_endpoints_i, available_endpoints_i)])
      value
    }
    list_lequi.tol[[i]] <- align_bounds(list_lequi.tol[[i]])
    list_uequi.tol[[i]] <- align_bounds(list_uequi.tol[[i]])
    if (is.null(names(list_lequi.tol[[i]])) ||
        is.null(names(list_uequi.tol[[i]]))){
      names(list_lequi.tol[[i]]) <- names(list_uequi.tol[[i]]) <- selected_endpoints_i
    }
    if (length(list_lequi.tol[[i]]) != length(list_y_comparator[[i]]) ||
        length(list_uequi.tol[[i]]) != length(list_y_comparator[[i]])) {
      stop(sprintf(
        "Equivalence bounds for comparator %s must have one lower and upper value for each of its %d endpoints.",
        paste(list_comparator[[i]], collapse = " vs "),
        length(list_y_comparator[[i]])
      ))
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



  # k is defined per comparator, after list_y_comparator has been resolved.
  kmax <- lengths(list_y_comparator)
  k_scalar_input <- !is.list(k) && length(k) == 1L && !is.na(k)
  if (is.list(k)) k <- unlist(k, use.names = FALSE)
  missing_k <- is.null(k) || !length(k) || all(is.na(k))
  if (missing_k) {
    k <- kmax
  } else {
    if (!is.numeric(k) || anyNA(k))
      stop("'k' must be a positive integer, one value per comparator or one value for all comparators.")
    if (length(k) == 1L) k <- rep(k, length(kmax))
    if (length(k) != length(kmax))
      stop("'k' must have one value per comparator or be a single value.")
    if (any(!is.finite(k)) || any(k != as.integer(k)) || any(k < 1L))
      stop("'k' must contain positive integers.")
  }
  oversized <- which(k > kmax)
  if (length(oversized)) {
    warning("'k' is larger than the number of selected endpoints for comparator(s) ",
            paste(oversized, collapse = ", "), "; setting it to the maximum possible value.",
            call. = FALSE)
    k[oversized] <- kmax[oversized]
  }
  if (k_scalar_input && length(unique(kmax)) > 1L) {
    warning("The selected endpoint counts differ by comparator; the supplied scalar 'k' is interpreted separately for each comparator and capped at that comparator's m.",
            call. = FALSE)
  }
  if (adjust == "seq" && type_info$active) {
    selected_types <- type_y[selected_endpoints]
    secondary <- selected_types == 2L
    if (any(secondary))
      weight_seq[selected_endpoints[secondary]] <- min(k) / sum(secondary)
  }
  independent <- all(vapply(seq_along(list_comparator), function(i) {
    arms <- list_comparator[[i]]
    endpoints <- list_y_comparator[[i]]
    is_diagonal <- function(x) {
      x <- as.matrix(x)
      if (nrow(x) <= 1L) return(TRUE)
      off_diagonal <- x
      diag(off_diagonal) <- 0
      max(abs(off_diagonal)) <= sqrt(.Machine$double.eps)
    }
    is_diagonal(varcov_list[[arms[[1L]]]][endpoints, endpoints, drop = FALSE]) &&
      is_diagonal(varcov_list[[arms[[2L]]]][endpoints, endpoints, drop = FALSE])
  }, logical(1)))
  if (isTRUE(.warn_redundant_bon))
    .warn_adjustment_configuration(
      k = k, m = kmax, adjust = adjust,
      type_y = if (type_info$active) type_y[selected_endpoints] else NULL,
      type_y_supplied = type_info$supplied,
      n_comparators = length(list_comparator), independent = independent,
      context = "selected endpoints"
    )


  # Save endpoints related information on a parameter list

  param <- list(mu = mu_list, varcov = varcov_list, sigmaB = sigmaB,
                TAR_list = TAR_list, type_y = type_y, weight_seq = weight_seq,
                arm_names = arm_names,  ynames_list = ynames_list,
                list_comparator = list_comparator,
                list_y_comparator = list_y_comparator,
                list_lequi.tol = list_lequi.tol,
                list_uequi.tol = list_uequi.tol,
                Eper = Eper, Eco = Eco)

  if (distribution == "lnorm" && ctype == "DOM"){
    stop("DOM is not supported for lognormal outcomes: the current implementation defines lognormal testing as a ratio-of-means test on the log scale. Use distribution = 'normal' with ctype = 'DOM', or distribution = 'lognormal' with ctype = 'ROM'.")
  }

  param.d <- list(nsim=nsim,
                  power=power,
                  alpha=alpha,
                  dtype=dtype,
                  ctype=ctype,
                  distribution=distribution,
                  vareq=vareq,
                  k=k,
                  adjust=adjust,
                  dropout=dropout,list_lequi.tol=list_lequi.tol, list_uequi.tol=list_uequi.tol,
                  seed=seed, ncores=ncores,
                  optimization_method=optimization_method, lower=lower, upper=upper,
                  step.power=step.power, step.up=step.up, pos.side=pos.side,
                  maxiter=maxiter, verbose=verbose)

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

  if (keep_sim_data) {
    selected_response <- out$response
    selected_n <- if (!is.null(selected_response) &&
                       !is.logical(selected_response) &&
                       !is.na(selected_response$n_iter[[1L]])) {
      as.integer(selected_response$n_iter[[1L]])
    } else {
      as.integer(lower)
    }
    out$sim_data <- .retained_continuous_data(
      param = param, param.d = param.d, n = selected_n,
      nsim = nsim, seed = seed
    )
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
#' @keywords internal
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
#' This function derives endpoint names (\code{ynames_list}) from \code{mu_list} if \code{ynames_list}
#' is missing. If \code{ynames_list} is already provided, it confirms the names to the user when
#' \code{verbose} is set to \code{TRUE}.
#'
#' @param ynames_list Optional list of vectors with endpoint names for each arm.
#' @param mu_list Named list of means per treatment arm, where names can be used as endpoint names.
#' @param verbose Logical, if \code{TRUE}, displays messages about the derivation process.
#'
#' @return A list of endpoint names for each arm.
#' @keywords internal
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

#' Derive and Validate Treatment Allocation Rate (TAR)
#'
#' This function validates and adjusts the treatment allocation rate (\code{TAR}) to ensure it is correctly specified
#' for the given number of treatment arms (\code{n_arms}). If \code{TAR} is missing or NULL, it is assigned a default
#' vector of ones, ensuring equal allocation across all arms. The function also handles cases where \code{TAR}
#' is shorter than \code{n_arms}, contains NA values, or has invalid values.
#'
#' @param TAR Optional numeric vector specifying the allocation rate for each treatment arm. If missing, a default
#' equal allocation rate is assigned.
#' @param arm_names Character vector specifying the names of the treatment arms. Used to name the elements of \code{TAR}.
#' @param verbose Logical, if \code{TRUE}, displays messages about the status of \code{TAR} derivation or assignment.
#'
#' @return A named list representing the treatment allocation rate for each arm.
#'
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @keywords internal
derive_allocation_rate <- function(TAR = NULL, arm_names, verbose = FALSE) {

  n_arms <- length(arm_names)

  # Handle missing or NULL TAR
  if (missing(TAR) || is.null(TAR)) {
    info_msg("Warning: TAR is missing or NULL. Setting TAR to a default vector of ones.", verbose = verbose)
    TAR <- rep(1, n_arms)
  }

  # Check for incorrect length
  if (length(TAR) > n_arms) {
    stop("Validation Error: TAR cannot exceed the number of treatment arms.")
  } else if (length(TAR) < n_arms) {
    warning("TAR length is shorter than the number of arms. Missing values will be replaced with 1.")
    TAR <- c(TAR, rep(1, n_arms - length(TAR)))
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

  # Assign names and return as a named list
  return(stats::setNames(as.list(TAR), arm_names))
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
#' @keywords internal
derive_varcov_list <- function(mu_list, sigma_list, ynames_list = NULL, varcov_list = NULL, cor_mat = NULL, rho = 0) {
  # Check if variance-covariance matrix is missing
  if (is.null(varcov_list) | any(is.na(varcov_list))) {
    if (any(is.na(sigma_list))) {
      stop("No variance-covariance matrix provided, and a standard deviation list is also missing. Either a variance-covariance matrix or a standard deviation list is required.")
    }

    # Validate lengths of inputs
    len_mu <- sapply(mu_list, length)
    len_sd <- sapply(sigma_list, length)

    # Check ynames_list, allowing it to be NULL
    if (!is.null(ynames_list)) {
      len_y <- sapply(ynames_list, length)

      if (any((len_mu != len_sd) | (len_mu != len_y))) {
        stop("In each arm, 'mu', 'sigma', and 'y_name' must have the same length.")
      }
    } else {
      if (any(len_mu != len_sd)) {
        stop("In each arm, 'mu' and 'sigma' must have the same length.")
      }
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

      # Matrix multiplication can introduce tiny, machine-precision
      # differences between mirrored off-diagonal entries. Restore the
      # mathematical symmetry before the positive-semidefinite check.
      varcov <- (varcov + t(varcov)) / 2

      # Add the matrix to the list
      varcov_list[[i]] <- varcov
    }
  }

  return(varcov_list)
}
