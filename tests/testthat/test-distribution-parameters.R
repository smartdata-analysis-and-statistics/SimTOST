test_that("plot_distribution supports single-arm parameter and multi-arm estimand plots", {
  skip_if_not_installed("MASS")
  result <- simPower(
    n = 15, distribution = "normal",
    mu_list = list(R = c(y1 = 1, y2 = 1), T = c(y1 = 1.1, y2 = 1.1)),
    sigma_list = list(R = c(y1 = .2, y2 = .2), T = c(y1 = .2, y2 = .2)),
    list_comparator = list(R_vs_T = c("R", "T")),
    list_lequi.tol = list(R_vs_T = c(y1 = -.3, y2 = -.3)),
    list_uequi.tol = list(R_vs_T = c(y1 = .3, y2 = .3)),
    ctype = "DOM", nsim = 20, ncores = 1, seed = 1,
    keep_sim_data = TRUE
  )

  expect_s3_class(plot_distribution(result, estimand = "mu", arms = "R",
                                    endpoints = "y1"), "ggplot")
  expect_s3_class(plot_distribution(result, estimand = "sigma", arms = "R",
                                    endpoints = "y1"), "ggplot")
  expect_s3_class(plot_distribution(result, estimand = "mu"), "ggplot")
  expect_s3_class(plot_distribution(result, estimand = "DOM",
                                    endpoints = "y1"), "ggplot")
  estimands <- SimTOST:::.estimand_data(result, endpoint = "y1")
  expect_equal(unique(estimands$truth), -0.1, tolerance = 1e-12)

  rom_result <- simPower(
    n = 15, distribution = "normal",
    mu_list = list(R = c(y1 = 1, y2 = 1), T = c(y1 = 1.1, y2 = 1.1)),
    sigma_list = list(R = c(y1 = .2, y2 = .2), T = c(y1 = .2, y2 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = .8, y2 = .8)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.25, y2 = 1.25)),
    ctype = "ROM", nsim = 20, ncores = 1, seed = 2,
    keep_sim_data = TRUE
  )
  rom_estimands <- SimTOST:::.estimand_data(rom_result, endpoint = "y1")
  expect_equal(unique(rom_estimands$truth), 1.1, tolerance = 1e-12)
})

test_that("plot_distribution validates parameter and estimand filters", {
  expect_error(
    plot_distribution(structure(list(), class = "simss"),
                      estimand = "DOM", arms = "R"),
    "Raw simulated data"
  )
})

test_that("plot_distribution uses the simplified estimand interface", {
  skip_if_not_installed("MASS")
  object <- simPower(
    n = 15, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1.1, y2 = 1.1), R = c(y1 = 1, y2 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2), R = c(y1 = .2, y2 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.5, y2 = -.5)),
    list_uequi.tol = list(T_vs_R = c(y1 = .5, y2 = .5)),
    nsim = 20, ncores = 1, seed = 11, keep_sim_data = TRUE
  )

  expect_s3_class(plot_distribution(object, estimand = "mu",
                                    arms = "T", endpoints = "y1"), "ggplot")
  expect_s3_class(plot_distribution(object, estimand = "mu",
                                    arms = "T", endpoints = "y1"), "ggplot")
  expect_s3_class(plot_distribution(object, estimand = "sigma"), "ggplot")
  expect_s3_class(plot_distribution(object, estimand = "dispersion"), "ggplot")
  expect_s3_class(plot_distribution(object, estimand = "DOM",
                                    arms = c("T", "R"), endpoints = "y2"), "ggplot")
  expect_s3_class(plot_distribution(object, estimand = "t_value",
                                    arms = "T"), "ggplot")
  expect_s3_class(plot_distribution(object, estimand = "correlation",
                                    endpoints = c("y1", "y2")), "ggplot")
  expect_warning(
    expect_s3_class(plot_distribution(object, estimand = "ROM"), "ggplot"),
    "incompatible.*ctype = 'DOM'.*Use plot_distribution.*estimand = \\\"dom\\\".*ctype")
  expect_error(plot_distribution(object, estimand = "mean", arms = "unknown"),
               "Unknown arm")
  expect_error(plot_distribution(object, estimand = "correlation", type = "qq"),
               "not available")
})

test_that("t-value plots retain all comparators selected by multiple arms", {
  object <- simPower(
    n = 12, distribution = "norm", ctype = "DOM",
    mu_list = list(T = c(y1 = 1.1), R1 = c(y1 = 1), R2 = c(y1 = 1.05)),
    sigma_list = list(T = c(y1 = .2), R1 = c(y1 = .2), R2 = c(y1 = .2)),
    list_comparator = list(T_vs_R1 = c("T", "R1"),
                           T_vs_R2 = c("T", "R2")),
    list_lequi.tol = list(T_vs_R1 = c(y1 = -.5), T_vs_R2 = c(y1 = -.5)),
    list_uequi.tol = list(T_vs_R1 = c(y1 = .5), T_vs_R2 = c(y1 = .5)),
    nsim = 10, ncores = 1, seed = 15, keep_sim_data = TRUE
  )
  plot <- plot_distribution(object, estimand = "t_value",
                             arms = c("T", "R1", "R2"))
  expect_setequal(unique(as.character(plot$data$Comparator)),
                  c("T_vs_R1", "T_vs_R2"))
})

test_that("plot_distribution displays count rate ratios with RR", {
  object <- simPower(
    n = 15, distribution = "Poisson",
    rate_list = list(T = c(y1 = .2, y2 = .2), R = c(y1 = .2, y2 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = .8, y2 = .8)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.25, y2 = 1.25)),
    nsim = 20, ncores = 1, seed = 12, keep_sim_data = TRUE
  )
  plot <- plot_distribution(object, estimand = "RR", endpoints = "y1")
  expect_s3_class(plot, "ggplot")
  expect_equal(unique(as.character(plot$data$endpoint)), "y1")
  expect_error(plot_distribution(object, estimand = "DOM"),
               "continuous")
})

test_that("plot_distribution accepts sampleSize simss objects", {
  object <- sampleSize(
    distribution = "norm",
    mu_list = list(T = c(y1 = 1.1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.5)),
    list_uequi.tol = list(T_vs_R = c(y1 = .5)),
    ctype = "DOM", nsim = 10, ncores = 1,
    lower = 2, upper = 20, keep_sim_data = TRUE
  )
  expect_s3_class(plot_distribution(object, estimand = "mu"), "ggplot")
  expect_s3_class(plot_distribution(object, estimand = "DOM"), "ggplot")
})

test_that("distribution levels use canonical R names and accept aliases", {
  common <- list(
    n = 10,
    mu_list = list(T = c(y1 = 1.1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = .8)),
    list_uequi.tol = list(T_vs_R = c(y1 = 1.25)),
    nsim = 10, seed = 42
  )
  canonical <- do.call(simPower, c(common, distribution = "lnorm"))
  legacy <- do.call(simPower, c(common, distribution = "log normal"))
  r_alias <- do.call(simPower, c(common, distribution = "lnorm"))
  expect_identical(canonical$distribution, "lnorm")
  expect_identical(legacy$distribution, "lnorm")
  expect_identical(r_alias$distribution, "lnorm")
  expect_error(do.call(simPower, c(common, distribution = "gamma")),
               "distribution")
})
