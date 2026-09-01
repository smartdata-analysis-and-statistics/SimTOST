test_that("sampleSize uses supplied endpoints per comparator", {
  mu <- list(T = c(y1 = 1, y2 = 1, y3 = 1),
             R = c(y1 = 1, y2 = 1, y3 = 1))
  sigma <- list(T = c(y1 = .2, y2 = .2, y3 = .2),
                R = c(y1 = .2, y2 = .2, y3 = .2))
  comparators <- list(T_vs_R = c("T", "R"))
  explicit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM", mu_list = mu,
    sigma_list = sigma, list_comparator = comparators,
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = NA, adjust = "no", lower = 2, upper = 2,
    optimization_method = "step-by-step", nsim = 5, ncores = 1, seed = 1
  )
  expect_identical(explicit$param$list_y_comparator[[1]], c("y1", "y2"))
  expect_equal(unname(explicit$param.d$k), 2L)
})

test_that("sampleSize infers all common endpoints when list_y_comparator is omitted", {
  mu <- list(T = c(y1 = 1, y2 = 1, y3 = 1),
             R = c(y1 = 1, y2 = 1, y3 = 1))
  sigma <- list(T = c(y1 = .2, y2 = .2, y3 = .2),
                R = c(y1 = .2, y2 = .2, y3 = .2))
  inferred <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM", mu_list = mu,
    sigma_list = sigma, list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2, y3 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2, y3 = .2)),
    k = NA, adjust = "no", lower = 2, upper = 2,
    optimization_method = "step-by-step", nsim = 5, ncores = 1, seed = 1
  )
  expect_identical(inferred$param$list_y_comparator[[1]], c("y1", "y2", "y3"))
  expect_equal(unname(inferred$param.d$k), 3L)
})

test_that("omitted endpoint selections warn when only a common subset is available", {
  fit <- NULL
  expect_warning(
    fit <- sampleSize(
      power = .5, distribution = "normal", ctype = "DOM",
      mu_list = list(T = c(y1 = 1, y2 = 1, y3 = 1),
                     R = c(y1 = 1, y2 = 1)),
      sigma_list = list(T = c(y1 = .2, y2 = .2, y3 = .2),
                        R = c(y1 = .2, y2 = .2)),
      list_comparator = list(T_vs_R = c("T", "R")),
      list_lequi.tol = list(T_vs_R = c(-.2, -.2)),
      list_uequi.tol = list(T_vs_R = c(.2, .2)),
      lower = 2, upper = 2, optimization_method = "step-by-step",
      nsim = 3, ncores = 1, seed = 10
    ),
    "effective endpoint count m is 2"
  )
  expect_identical(fit$param$list_y_comparator[[1]], c("y1", "y2"))
})

test_that("type1Error uses inferred or supplied endpoints for m and k", {
  common <- list(
    null = "lower", joint = TRUE, n = 10, distribution = "normal",
    ctype = "DOM", mu_list = list(T = c(y1 = 1, y2 = 1, y3 = 1),
                                   R = c(y1 = 1, y2 = 1, y3 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2, y3 = .2),
                      R = c(y1 = .2, y2 = .2, y3 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2, y3 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2, y3 = .2)),
    k = NA, nsim = 5, ncores = 1, seed = 2
  )
  explicit_args <- common
  explicit_args$list_y_comparator <- list(T_vs_R = c("y1", "y2"))
  explicit_args$list_lequi.tol <- list(T_vs_R = c(y1 = -.2, y2 = -.2))
  explicit_args$list_uequi.tol <- list(T_vs_R = c(y1 = .2, y2 = .2))
  explicit_args$k <- 2
  explicit_args$adjust <- "bon"
  explicit <- NULL
  expect_warning(explicit <- do.call(type1Error, explicit_args), "not necessary")
  inferred <- do.call(type1Error, common)
  expect_equal(nrow(explicit$table), 3L)
  expect_equal(nrow(inferred$table), 7L)
})

test_that("oversized k is warned about and capped at comparator m", {
  mu <- list(T = c(y1 = 1, y2 = 1, y3 = 1),
             R = c(y1 = 1, y2 = 1, y3 = 1))
  sigma <- list(T = c(y1 = .2, y2 = .2, y3 = .2),
                R = c(y1 = .2, y2 = .2, y3 = .2))
  fit <- NULL
  expect_warning(fit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM", mu_list = mu,
    sigma_list = sigma, list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = 3, adjust = "no", lower = 2, upper = 2,
    optimization_method = "step-by-step", nsim = 5, ncores = 1, seed = 3
  ), "larger than the number")
  expect_equal(unname(fit$param.d$k), 2L)

  args <- list(
    null = "lower", n = 10, distribution = "normal", ctype = "DOM",
    mu_list = mu, sigma_list = sigma,
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = 3, nsim = 5, ncores = 1, seed = 3
  )
  type1 <- NULL
  expect_warning(type1 <- do.call(type1Error, args), "larger than the number")
  expect_equal(unname(type1$result$param.d$k), 2L)
})

test_that("simPower applies endpoint-specific k validation", {
  args <- list(
    n = 10, distribution = "normal", ctype = "DOM",
    mu_list = list(T = c(y1 = 1, y2 = 1, y3 = 1),
                   R = c(y1 = 1, y2 = 1, y3 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2, y3 = .2),
                      R = c(y1 = .2, y2 = .2, y3 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = 3, adjust = "no", nsim = 5, ncores = 1, seed = 4
  )
  fit <- NULL
  expect_warning(fit <- do.call(simPower, args), "larger than the number")
  expect_equal(unname(fit$result$param.d$k), 2L)
  expect_equal(fit$result$param$list_y_comparator[[1]], c("y1", "y2"))
})

test_that("Bonferroni is retained but warns when k equals selected m", {
  mu <- list(T = c(y1 = 1, y2 = 1), R = c(y1 = 1, y2 = 1))
  sigma <- list(T = c(y1 = .2, y2 = .2), R = c(y1 = .2, y2 = .2))
  fit <- NULL
  expect_warning(fit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM", mu_list = mu,
    sigma_list = sigma, list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = 2, adjust = "bon", lower = 2, upper = 2,
    optimization_method = "step-by-step", nsim = 5, ncores = 1, seed = 5
  ), "not necessary")
  expect_equal(unname(fit$param.d$k), 2L)
  expect_identical(fit$param.d$adjust, "bon")
  expect_equal(SimTOST:::.endpoint_alpha(.05, 2, 2, "bon"), c(.025, .025))
})

test_that("t-adjustment uses m minus k plus one", {
  expect_equal(
    SimTOST:::.endpoint_alpha(.05, m = 3, k = 2, adjust = "t"),
    rep(.025, 3)
  )
  expect_equal(
    SimTOST:::.normalize_adjustment("partial-conjunction"), "t"
  )
  expect_warning(
    SimTOST:::.warn_adjustment_configuration(k = 3, m = 3, adjust = "t"),
    "has no effect"
  )
  printed <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM",
    mu_list = list(T = c(y1 = 1, y2 = 1), R = c(y1 = 1, y2 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2), R = c(y1 = .2, y2 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = 1, adjust = "t", lower = 2, upper = 2,
    optimization_method = "step-by-step", nsim = 5, ncores = 1, seed = 6
  )
  expect_output(print(printed), "t-adjustment")
  expect_error(
    SimTOST:::.endpoint_alpha(.05, m = 2, k = 3, adjust = "t"),
    "between 1 and 'm'"
  )
  expect_error(
    SimTOST:::.normalize_adjustment("not-an-adjustment"),
    "one of 'no', 'bon', 'sid', 'k', 't', or 'seq'"
  )
})

test_that("count k validation rejects invalid values", {
  expect_error(SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), k = 0, nsim = 5
  ), "positive integers")
  expect_error(SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), k = 1.5, nsim = 5
  ), "positive integers")
})

test_that("adjustment warnings distinguish co-primary and k-of-m decisions", {
  expect_warning(
    SimTOST:::.warn_adjustment_configuration(
      k = 2, m = 2, adjust = "k"
    ),
    "has no effect"
  )
  expect_warning(
    SimTOST:::.warn_adjustment_configuration(
      k = 2, m = 2, adjust = "sid"
    ),
    "not necessary"
  )
  expect_warning(
    SimTOST:::.warn_adjustment_configuration(
      k = 1, m = 2, adjust = "no"
    ),
    "does not provide multiplicity control"
  )
})

test_that("type_y is validated, aligned, and scoped to selected endpoints", {
  aligned <- SimTOST:::.prepare_type_y(
    type_y = c(y3 = 2, y1 = 1),
    all_endpoints = c("y1", "y2", "y3"),
    selected_endpoints = c("y3", "y1"), adjust = "seq"
  )
  expect_identical(unname(aligned$type_y), c(1L, 1L, 2L))
  expect_identical(
    SimTOST:::.type_y_for_endpoints(aligned$type_y, c("y3", "y1"), "seq"),
    c(2L, 1L)
  )
  expect_warning(
    SimTOST:::.prepare_type_y(c(1, 2), c("y1", "y2"),
                               c("y1", "y2"), "bon"),
    "used only with adjust = 'seq'"
  )
  expect_error(
    SimTOST:::.prepare_type_y(c(y1 = 3), c("y1"), c("y1"), "seq"),
    "only integer values 1"
  )
  expect_error(
    SimTOST:::.prepare_type_y(c(y1 = 1), c("y1", "y2"),
                               c("y1", "y2"), "seq"),
    "missing: y2"
  )
})

test_that("sequential adjustment follows named endpoint selections", {
  mu <- list(
    T = c(y1 = 1, y2 = 1, y3 = 1),
    R = c(y1 = 1, y2 = 1, y3 = 1),
    S = c(y1 = 1, y2 = 1, y3 = 1)
  )
  sigma <- list(
    T = c(y1 = .2, y2 = .2, y3 = .2),
    R = c(y1 = .2, y2 = .2, y3 = .2),
    S = c(y1 = .2, y2 = .2, y3 = .2)
  )
  fit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM",
    mu_list = mu, sigma_list = sigma,
    list_comparator = list(TR = c("T", "R"), TS = c("T", "S")),
    list_y_comparator = list(TR = c("y1", "y2"), TS = c("y1", "y3")),
    list_lequi.tol = list(TR = c(-.2, -.2), TS = c(-.2, -.2)),
    list_uequi.tol = list(TR = c(.2, .2), TS = c(.2, .2)),
    type_y = c(y1 = 1, y2 = 2, y3 = 2), k = 1, adjust = "seq",
    lower = 2, upper = 2, optimization_method = "step-by-step",
    nsim = 5, ncores = 1, seed = 11
  )
  expect_identical(unname(fit$param$type_y), c(1L, 2L, 2L))
  expect_equal(unname(fit$param.d$k), c(1L, 1L))
  expect_identical(fit$param$list_y_comparator[[2]], c("y1", "y3"))
})

test_that("sequential adjustment infers endpoints when endpoint selection is omitted", {
  fit <- sampleSize(
    power = .1, distribution = "normal", ctype = "DOM",
    mu_list = list(T = c(y1 = 1, y2 = 1, y3 = 1),
                   R = c(y1 = 1, y2 = 1, y3 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2, y3 = .2),
                      R = c(y1 = .2, y2 = .2, y3 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2, y3 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2, y3 = .2)),
    rho = .6, type_y = c(y1 = 2, y2 = 2, y3 = 1),
    k = c(T_vs_R = 1), adjust = "seq",
    lower = 2, upper = 2, optimization_method = "step-by-step",
    nsim = 3, ncores = 1, seed = 12
  )

  expect_identical(fit$param$list_y_comparator[[1]], c("y1", "y2", "y3"))
  expect_identical(unname(fit$param$type_y), c(2L, 2L, 1L))
})
