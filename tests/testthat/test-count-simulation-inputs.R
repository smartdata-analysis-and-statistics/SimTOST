test_that("retained count outcomes reproduce supplied rates and correlations", {
  skip_if_not_installed("MASS")

  target_rates <- c(y1 = 0.50, y2 = 0.70)
  target_corr <- matrix(c(1, 0.50, 0.50, 1), nrow = 2,
                        dimnames = list(names(target_rates), names(target_rates)))

  for (distribution in c("poisson", "negative binomial")) {
    result <- simPower(
      n = 60,
      distribution = distribution,
      rate_list = list(T = target_rates, R = target_rates),
      exposure = 10,
      dispersion = 0.10,
      cor_mat = target_corr,
      list_comparator = list(T_vs_R = c("T", "R")),
      list_lequi.tol = list(T_vs_R = c(y1 = 0.80, y2 = 0.80)),
      list_uequi.tol = list(T_vs_R = c(y1 = 1.25, y2 = 1.25)),
      nsim = 80,
      seed = 2025,
      ncores = 1,
      keep_sim_data = TRUE
    )

    expect_true(is.data.frame(result$sim_data))
    expect_true(all(c("trial", "subject", "arm", "endpoint", "value",
                      "exposure") %in% names(result$sim_data)))

    observed <- result$sim_data[result$sim_data$arm == "T", ]
    observed_rate <- with(observed, tapply(value / exposure, endpoint, mean))
    expect_equal(as.numeric(observed_rate[names(target_rates)]),
                 unname(target_rates), tolerance = 0.06)

    wide <- reshape(
      observed[, c("trial", "subject", "endpoint", "value")],
      idvar = c("trial", "subject"), timevar = "endpoint",
      direction = "wide"
    )
    observed_corr <- stats::cor(wide$value.y1, wide$value.y2)
    expect_equal(observed_corr, target_corr["y1", "y2"], tolerance = 0.12)
  }
})

test_that("count distribution diagnostics preserve the requested arm order", {
  skip_if_not_installed("ggplot2")
  result <- simPower(
    n = 20, distribution = "poisson",
    rate_list = list(TEST = c(y1 = 0.5, y2 = 0.7),
                     REF = c(y1 = 0.5, y2 = 0.7)),
    exposure = 10,
    cor_mat = matrix(c(1, 0.5, 0.5, 1), nrow = 2,
                     dimnames = list(c("y1", "y2"), c("y1", "y2"))),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_lequi.tol = list(TEST_vs_REF = c(y1 = 0.8, y2 = 0.8)),
    list_uequi.tol = list(TEST_vs_REF = c(y1 = 1.25, y2 = 1.25)),
    nsim = 10, seed = 2026, ncores = 1, keep_sim_data = TRUE
  )

  rate_plot <- ggplot2::ggplot_build(
    plot_distribution(result, estimand = "rate", arms = c("TEST", "REF"))
  )
  correlation_plot <- ggplot2::ggplot_build(
    plot_distribution(result, estimand = "correlation", arms = c("TEST", "REF"))
  )
  expect_identical(unique(as.character(rate_plot$layout$layout$arm)), c("TEST", "REF"))
  expect_identical(unique(as.character(correlation_plot$layout$layout$arm)), c("TEST", "REF"))
})
