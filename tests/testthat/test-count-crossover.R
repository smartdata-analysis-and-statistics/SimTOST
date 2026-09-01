test_that("2x2 count power uses the paired crossover analysis", {
  poisson <- simPower(
    n = 30,
    distribution = "Poisson",
    rate_list = list(TEST = 0.21, REF = 0.20),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_lequi.tol = list(TEST_vs_REF = 0.80),
    list_uequi.tol = list(TEST_vs_REF = 1.25),
    exposure = 1,
    dtype = "2x2",
    sigmaB = 0.4,
    Eper = c(0.15, -0.15),
    Eco = c(0.02, 0.05),
    nsim = 100,
    seed = 123
  )

  negative_binomial <- simPower(
    n = 30,
    distribution = "Negative Binomial",
    rate_list = list(TEST = 0.21, REF = 0.20),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_lequi.tol = list(TEST_vs_REF = 0.80),
    list_uequi.tol = list(TEST_vs_REF = 1.25),
    exposure = 1,
    dispersion = 0.01,
    dtype = "2x2",
    sigmaB = 0.4,
    Eper = c(0.15, -0.15),
    Eco = c(0.02, 0.05),
    nsim = 100,
    seed = 123
  )

  expect_s3_class(poisson, "countpower")
  expect_s3_class(negative_binomial, "countpower")
  expect_true(is.finite(poisson$power))
  expect_true(is.finite(negative_binomial$power))
  expect_equal(poisson$n_total, 60)
  expect_equal(negative_binomial$n_total, 60)
})

test_that("2x2 count crossover validates sequence dropout", {
  expect_error(
    simPower(
      n = 20,
      distribution = "Poisson",
      rate_list = list(TEST = 0.2, REF = 0.2),
      list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
      list_lequi.tol = list(TEST_vs_REF = 0.8),
      list_uequi.tol = list(TEST_vs_REF = 1.25),
      dtype = "2x2", dropout = c(0, 1), nsim = 10
    ),
    "dropout"
  )
})

test_that("multi-endpoint count simulations use endpoint correlation", {
  endpoint_corr <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  result <- simPower(
    n = 40,
    distribution = "Poisson",
    rate_list = list(TEST = c(5, 6), REF = c(5, 6)),
    list_comparator = list(T_vs_R = c("TEST", "REF")),
    list_lequi.tol = list(T_vs_R = c(0.8, 0.8)),
    list_uequi.tol = list(T_vs_R = c(1.25, 1.25)),
    cor_mat = endpoint_corr,
    dtype = "parallel",
    nsim = 50,
    seed = 123
  )

  expect_s3_class(result, "simpower")
  expect_true(is.finite(result$power))
})

test_that("count sample-size search supports fast and sequential methods", {
  common <- list(
    power = 0.80,
    distribution = "Poisson",
    rate_list = list(TEST = 5.1, REF = 5),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_lequi.tol = list(TEST_vs_REF = 0.80),
    list_uequi.tol = list(TEST_vs_REF = 1.25),
    lower = 20, upper = 100, nsim = 50, seed = 123
  )
  fast <- do.call(sampleSize, c(common, list(optimization_method = "fast")))
  sequential <- do.call(sampleSize,
                         c(common, list(optimization_method = "step-by-step")))

  expect_equal(fast$n_per_arm, sequential$n_per_arm)
})

test_that("parallel and crossover count designs are distinct", {
  common <- list(
    n = 40,
    distribution = "Negative Binomial",
    rate_list = list(TEST = 0.21, REF = 0.20),
    list_comparator = list(TEST_vs_REF = c("TEST", "REF")),
    list_lequi.tol = list(TEST_vs_REF = 0.80),
    list_uequi.tol = list(TEST_vs_REF = 1.25),
    dispersion = 0.01,
    nsim = 100,
    seed = 321
  )
  parallel <- do.call(simPower, c(common, list(dtype = "parallel")))
  crossover <- do.call(simPower, c(common, list(
    dtype = "2x2", sigmaB = 0.4, Eper = c(0.2, -0.1), Eco = c(0.02, 0.05)
  )))

  expect_true(is.finite(parallel$power))
  expect_true(is.finite(crossover$power))
  expect_false(isTRUE(all.equal(parallel$power, crossover$power)))
})
