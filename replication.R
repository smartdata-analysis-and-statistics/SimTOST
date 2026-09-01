# Reproducibility script for SimTOST examples and Figure 2
#
# Run from the package root with:
#   Rscript replication.R
#
# The script intentionally uses modest nsim values for a quick smoke run.
# Increase `nsim` below for publication-quality Monte Carlo precision.

if (!requireNamespace("SimTOST", quietly = TRUE))
  stop("Install SimTOST before running this script.")

set.seed(20240818)
out_dir <- file.path("replication-output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

endpoint <- "AUCinf"
mu_R <- setNames(1, endpoint)
mu_T <- setNames(1.05, endpoint)
sd_R <- setNames(0.30, endpoint)
sd_T <- setNames(0.30, endpoint)
lower <- setNames(0.80, endpoint)
upper <- setNames(1.25, endpoint)

# 1. The published Mielke-style reference calculation.
mielke <- SimTOST::sampleSize_Mielke(
  power = 0.80, Nmax = 200, m = 1, k = 1, rho = 0,
  sigma = 0.30, true.diff = log(1.05), equi.tol = log(1.25),
  design = "parallel", alpha = 0.05, nsim = 200, seed = 20240818
)

# 2. Main continuous sample-size interface.
continuous_ss <- SimTOST::sampleSize(
  power = 0.80, distribution = "Log Normal", alpha = 0.05,
  mu_list = list(R = mu_R, T = mu_T),
  sigma_list = list(R = sd_R, T = sd_T),
  list_comparator = list(T_vs_R = c("R", "T")),
  list_y_comparator = list(T_vs_R = endpoint),
  list_lequi.tol = list(T_vs_R = lower),
  list_uequi.tol = list(T_vs_R = upper),
  dtype = "parallel", ctype = "ROM", lower = 10, upper = 200,
  nsim = 200, seed = 20240818, ncores = 1
)

# 3. Fixed-sample-size power for the continuous outcome.
continuous_power <- SimTOST::simPower(
  n = 50, distribution = "Log Normal",
  mu_list = list(R = mu_R, T = mu_T),
  sigma_list = list(R = sd_R, T = sd_T),
  list_comparator = list(T_vs_R = c("R", "T")),
  list_y_comparator = list(T_vs_R = endpoint),
  list_lequi.tol = list(T_vs_R = lower),
  list_uequi.tol = list(T_vs_R = upper),
  dtype = "parallel", ctype = "ROM", nsim = 500, seed = 20240818
)

# 4. Fixed-sample-size power for a Poisson count outcome. The count API uses
# the same list-based comparison and margin arguments as sampleSize().
count_power <- SimTOST::simPower(
  n = 100, distribution = "Poisson",
  rate_list = list(TEST = setNames(0.21, endpoint),
                   REF = setNames(0.20, endpoint)),
  list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
  list_lequi.tol = list(TEST_vs_REF = setNames(0.80, endpoint)),
  list_uequi.tol = list(TEST_vs_REF = setNames(1.25, endpoint)),
  dtype = "parallel", nsim = 500, seed = 20240818, ncores = 1
)

# 5. Figure 2: reproducible fixed-sample-size power curve. The exact inputs,
# seed, and grid are recorded here so the manuscript figure can be regenerated
# without relying on an interactive session.
n_grid <- seq(20, 100, by = 20)
figure2 <- do.call(rbind, lapply(n_grid, function(n_i) {
  fit <- SimTOST::simPower(
    n = n_i, distribution = "Log Normal",
    mu_list = list(R = mu_R, T = mu_T),
    sigma_list = list(R = sd_R, T = sd_T),
    list_comparator = list(T_vs_R = c("R", "T")),
    list_y_comparator = list(T_vs_R = endpoint),
    list_lequi.tol = list(T_vs_R = lower),
    list_uequi.tol = list(T_vs_R = upper),
    dtype = "parallel", ctype = "ROM", nsim = 500,
    seed = 20240818 + n_i
  )
  data.frame(n_per_arm = n_i, n_total = 2 * n_i,
             power = fit$power, power_LCI = fit$power_LCI,
             power_UCI = fit$power_UCI)
}))

utils::write.csv(figure2, file.path(out_dir, "figure2_power_curve.csv"),
                 row.names = FALSE)
grDevices::pdf(file.path(out_dir, "figure2_power_curve.pdf"),
               width = 6.5, height = 4.5, useDingbats = FALSE)
graphics::plot(figure2$n_total, figure2$power, type = "o", pch = 19,
               col = "#0072B2", ylim = c(0, 1), xlab = "Total sample size",
               ylab = "Estimated power", main = "Figure 2: Power by sample size")
graphics::arrows(figure2$n_total, figure2$power_LCI,
                 figure2$n_total, figure2$power_UCI,
                 code = 3, angle = 90, length = 0.05, col = "#0072B2")
graphics::abline(h = 0.80, lty = 2, col = "#4B5563")
grDevices::dev.off()

saveRDS(list(mielke = mielke, continuous_ss = continuous_ss,
             continuous_power = continuous_power, count_power = count_power,
             figure2 = figure2), file.path(out_dir, "replication-results.rds"))

print(summary(mielke))
print(summary(continuous_power))
print(summary(count_power))
message("Replication outputs written to ", normalizePath(out_dir))
