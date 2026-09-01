# Optional raw-data retention for distribution diagnostics. This is deliberately
# separate from the compiled testing path so the default memory footprint is
# unchanged.
.check_keep_sim_data <- function(keep_sim_data) {
  if (!is.logical(keep_sim_data) || length(keep_sim_data) != 1L ||
      is.na(keep_sim_data))
    stop("'keep_sim_data' must be TRUE or FALSE.")
  keep_sim_data
}

.restore_rng <- function(seed, fun) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed) && is.finite(seed)) set.seed(as.integer(seed))
  force(fun())
}

.retained_continuous_data <- function(param, param.d, n, nsim, seed) {
  arms <- param$arm_names
  endpoints <- param$ynames_list[[arms[[1L]]]]
  if (is.null(endpoints)) endpoints <- paste0("y", seq_along(param$mu[[1L]]))
  m <- length(endpoints)
  rows <- vector("list", nsim * length(arms))
  index <- 0L
  .restore_rng(seed, function() {
    for (trial in seq_len(nsim)) {
      for (arm in arms) {
        index <<- index + 1L
        n_arm <- if (param.d$dtype == "parallel") {
          ceiling(n * param$TAR_list[[arm]])
        } else {
          ceiling(n / 2)
        }
        n_arm <- max(2L, as.integer(n_arm))
        mu <- as.numeric(param$mu[[arm]])
        Sigma <- as.matrix(param$varcov[[arm]])
        if (identical(param.d$distribution, "lnorm")) {
          Sigma <- log(Sigma / (mu %*% t(mu)) + 1)
          mu <- log(mu) - 0.5 * diag(Sigma)
        }
        values <- MASS::mvrnorm(n = n_arm, mu = mu, Sigma = Sigma)
        if (identical(param.d$distribution, "lnorm"))
          values <- exp(values)
        rows[[index]] <<- data.frame(
          trial = trial, arm = arm,
          subject = rep(seq_len(n_arm), times = m),
          endpoint = rep(endpoints, each = n_arm),
          value = as.vector(values),
          distribution = param.d$distribution, design = param.d$dtype,
          stringsAsFactors = FALSE
        )
      }
    }
  })
  do.call(rbind, rows)
}

.count_parameter <- function(x, arm, endpoint_index, arms) {
  if (is.list(x)) {
    value <- x[[arm]]
  } else {
    value <- x
  }
  if (length(value) == 1L) return(as.numeric(value))
  if (!is.null(names(value))) return(as.numeric(value[[endpoint_index]]))
  as.numeric(value[[endpoint_index]])
}

.retained_count_data <- function(rate_list, exposure, dispersion, distribution,
                                 design, n, nsim, seed, endpoint_corr = NULL) {
  arms <- names(rate_list)
  endpoints <- names(rate_list[[1L]])
  if (is.null(endpoints)) endpoints <- paste0("y", seq_along(rate_list[[1L]]))
  m <- length(rate_list[[1L]])
  if (is.null(endpoint_corr)) endpoint_corr <- diag(m)
  endpoint_corr <- as.matrix(endpoint_corr)
  if (!all(dim(endpoint_corr) == c(m, m)) ||
      max(abs(endpoint_corr - t(endpoint_corr))) > 1e-10 ||
      max(abs(diag(endpoint_corr) - 1)) > 1e-10 ||
      min(eigen(endpoint_corr, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("'endpoint_corr' must be a positive-definite correlation matrix.")
  rows <- vector("list", nsim * length(arms))
  index <- 0L
  .restore_rng(seed, function() {
    for (trial in seq_len(nsim)) {
      for (arm in arms) {
        index <<- index + 1L
        n_arm <- if (design == "parallel") ceiling(n) else ceiling(n / 2)
        n_arm <- max(2L, as.integer(n_arm))
        latent <- MASS::mvrnorm(n = n_arm, mu = rep(0, m),
                                Sigma = endpoint_corr)
        uniforms <- stats::pnorm(latent)
        values <- matrix(numeric(0), nrow = n_arm, ncol = m)
        true_rates <- true_exposure <- true_dispersion <- numeric(m)
        for (j in seq_len(m)) {
          rate <- .count_parameter(rate_list, arm, j, arms)
          exposure_arm <- .count_parameter(exposure, arm, j, arms)
          dispersion_arm <- .count_parameter(dispersion, arm, j, arms)
          true_rates[[j]] <- rate
          true_exposure[[j]] <- exposure_arm
          true_dispersion[[j]] <- dispersion_arm
          mean_count <- rate * exposure_arm
          values[, j] <- if (distribution == "pois") {
            stats::qpois(uniforms[, j], mean_count)
          } else {
            stats::qnbinom(uniforms[, j], size = 1 / dispersion_arm,
                           mu = mean_count)
          }
        }
        rows[[index]] <<- data.frame(
          trial = trial, arm = arm,
          subject = rep(seq_len(n_arm), times = m),
          endpoint = rep(endpoints, each = n_arm),
          value = as.vector(values),
          true_rate = rep(true_rates, each = n_arm),
          exposure = rep(true_exposure, each = n_arm),
          true_dispersion = rep(true_dispersion, each = n_arm),
          distribution = distribution, design = design,
          stringsAsFactors = FALSE
        )
      }
    }
  })
  do.call(rbind, rows)
}
