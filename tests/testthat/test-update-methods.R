test_that("update simss preserves stored settings and overrides supplied ones", {
  fit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM",
    mu_list = list(T = c(y1 = 1, y2 = 1), R = c(y1 = 1, y2 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2), R = c(y1 = .2, y2 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = 1, adjust = "no", lower = 2, upper = 2,
    optimization_method = "step-by-step", nsim = 5, ncores = 1, seed = 31
  )
  updated <- update(fit, nsim = 7)
  expect_s3_class(updated, "simss")
  expect_equal(updated$param.d$nsim, 7)
  expect_equal(updated$param.d$k, fit$param.d$k)
  expect_identical(updated$param$list_comparator, fit$param$list_comparator)

  call <- update(fit, nsim = 7, evaluate = FALSE)
  expect_type(call, "language")
  expect_identical(as.character(call[[1L]]), "sampleSize")
  expect_identical(as.numeric(call$nsim), 7)
})

test_that("update simss rebuilds covariance matrices from supplied standard deviations", {
  fit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM",
    mu_list = list(T = c(y1 = 10), R = c(y1 = 10)),
    sigma_list = list(T = c(y1 = 1), R = c(y1 = 1)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = "y1"),
    list_lequi.tol = list(T_vs_R = -.2),
    list_uequi.tol = list(T_vs_R = .2),
    lower = 2, upper = 2, optimization_method = "step-by-step",
    nsim = 5, ncores = 1, seed = 37
  )

  updated <- update(
    fit,
    mu_list = list(T = 10, R = 10),
    sigma_list = list(T = 3, R = 4)
  )

  expect_equal(as.numeric(updated$param$varcov$T), 9)
  expect_equal(as.numeric(updated$param$varcov$R), 16)

  call <- update(
    fit,
    sigma_list = list(T = 3, R = 4),
    evaluate = FALSE
  )
  expect_null(call$varcov_list)
  expect_identical(as.numeric(call$sigma_list$T), 3)
})

test_that("update rejects a supplied standard-deviation list with invalid dimensions", {
  fit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM",
    mu_list = list(T = c(y1 = 10), R = c(y1 = 10)),
    sigma_list = list(T = c(y1 = 1), R = c(y1 = 1)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = "y1"),
    list_lequi.tol = list(T_vs_R = -.2),
    list_uequi.tol = list(T_vs_R = .2),
    lower = 2, upper = 2, optimization_method = "step-by-step",
    nsim = 5, ncores = 1, seed = 38
  )

  expect_error(
    update(fit, sigma_list = list(T = c(y1 = 1, y2 = 2), R = c(y1 = 1))),
    "same length"
  )
})

test_that("update infers endpoint selections when comparators are replaced", {
  fit <- sampleSize(
    power = .5, distribution = "normal", ctype = "DOM",
    arm_names = c("T", "R1", "R2"),
    mu_list = list(T = c(y1 = 1), R1 = c(y1 = 1), R2 = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R1 = c(y1 = .2), R2 = c(y1 = .2)),
    list_comparator = list(T_vs_R1 = c("T", "R1")),
    list_y_comparator = list(T_vs_R1 = "y1"),
    list_lequi.tol = list(T_vs_R1 = -.2),
    list_uequi.tol = list(T_vs_R1 = .2),
    lower = 2, upper = 2, optimization_method = "step-by-step",
    nsim = 5, ncores = 1, seed = 39
  )

  updated <- update(
    fit,
    list_comparator = list(
      T_vs_R1 = c("T", "R1"),
      T_vs_R2 = c("T", "R2")
    ),
    list_lequi.tol = list(T_vs_R1 = -.2, T_vs_R2 = -.2),
    list_uequi.tol = list(T_vs_R1 = .2, T_vs_R2 = .2)
  )

  expect_length(updated$param$list_y_comparator, 2)
  expect_identical(
    unname(vapply(updated$param$list_y_comparator, `[[`, character(1), 1L)),
    c("y1", "y1")
  )

  call <- update(
    fit,
    list_comparator = list(T_vs_R1 = c("T", "R1"), T_vs_R2 = c("T", "R2")),
    evaluate = FALSE
  )
  expect_true(length(call$list_y_comparator) == 1L && is.na(call$list_y_comparator))
})

test_that("update simpower preserves the fixed sample size and settings", {
  fit <- simPower(
    n = 10, distribution = "normal", ctype = "DOM",
    mu_list = list(T = c(y1 = 1, y2 = 1), R = c(y1 = 1, y2 = 1)),
    sigma_list = list(T = c(y1 = .2, y2 = .2), R = c(y1 = .2, y2 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.2, y2 = -.2)),
    list_uequi.tol = list(T_vs_R = c(y1 = .2, y2 = .2)),
    k = 1, adjust = "no", nsim = 5, ncores = 1, seed = 32
  )
  updated <- update(fit, nsim = 7, alpha = .04)
  expect_s3_class(updated, "simpower")
  expect_equal(updated$n, fit$n)
  expect_equal(updated$nsim, 7)
  expect_equal(updated$result$param.d$alpha, .04)
  expect_equal(updated$result$param.d$k, fit$result$param.d$k)
})

test_that("update supports standalone and unified count objects", {
  count_power <- SimTOST:::power_count(
    n_per_arm = 20, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), k = 1, nsim = 5,
    seed = 33, ncores = 1
  )
  updated_power <- update(count_power, nsim = 7)
  expect_s3_class(updated_power, "countpower")
  expect_equal(updated_power$nsim, 7)
  expect_equal(updated_power$n_per_arm, count_power$n_per_arm)
  expect_identical(updated_power$parameters$rate_test,
                   count_power$parameters$rate_test)

  count_ss <- SimTOST:::sampleSize_count(
    power = .1, rate_test = c(y1 = 5, y2 = 5),
    rate_reference = c(y1 = 5, y2 = 5), margin_lower = .1,
    margin_upper = 10, k = 1, lower = 2, upper = 10,
    optimization_method = "fast", nsim = 5, seed = 34, ncores = 1
  )
  updated_ss <- update(count_ss, nsim = 7)
  expect_s3_class(updated_ss, "countss")
  expect_equal(updated_ss$nsim, 7)
  expect_equal(updated_ss$parameters$rate_test, count_ss$parameters$rate_test)

  unified <- simPower(
    n = 20, distribution = "poisson",
    rate_list = list(T = c(y1 = 5, y2 = 5), R = c(y1 = 5, y2 = 5)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = .1, y2 = .1)),
    list_uequi.tol = list(T_vs_R = c(y1 = 10, y2 = 10)),
    k = 1, nsim = 5, seed = 35, ncores = 1
  )
  updated_unified <- update(unified, nsim = 7)
  expect_s3_class(updated_unified, "simpower")
  expect_equal(updated_unified$nsim, 7)
  expect_equal(updated_unified$n, unified$n)

  unified_ss <- sampleSize(
    power = .1, distribution = "poisson",
    rate_list = list(T = c(y1 = 5, y2 = 5), R = c(y1 = 5, y2 = 5)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_y_comparator = list(T_vs_R = c("y1", "y2")),
    list_lequi.tol = list(T_vs_R = c(y1 = .1, y2 = .1)),
    list_uequi.tol = list(T_vs_R = c(y1 = 10, y2 = 10)),
    k = 1, lower = 2, upper = 2, nsim = 5, seed = 36, ncores = 1
  )
  updated_unified_ss <- update(unified_ss, nsim = 7)
  expect_true(inherits(updated_unified_ss, "simss"))
  expect_true(inherits(updated_unified_ss, "countss"))
  expect_equal(updated_unified_ss$nsim, 7)
})

test_that("update rejects unsupported objects", {
  expect_error(update(list(), nsim = 10), "call component")
})
