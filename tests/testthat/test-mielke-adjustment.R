test_that("sampleSize_Mielke supports the t-adjustment", {
  fit <- sampleSize_Mielke(
    power = 0.01, Nmax = 10, m = 3, k = 2, rho = 0,
    sigma = 0.3, true.diff = log(1.05), equi.tol = log(1.25),
    design = "parallel", alpha = 0.05, adjust = "t",
    seed = 40, nsim = 5
  )
  expect_s3_class(fit, "simss_mielke")
  expect_identical(fit$parameters$adjust, "t")
  expect_true(is.finite(fit$power))
  expect_equal(SimTOST:::.normalize_adjustment("partial-conjunction"), "t")
  expect_error(
    sampleSize_Mielke(
      power = 0.01, Nmax = 4, m = 3, k = 2, rho = 0,
      sigma = 0.3, true.diff = log(1.05), equi.tol = log(1.25),
      design = "parallel", alpha = 0.05, adjust = "not-valid",
      nsim = 2
    ),
    "one of 'no', 'bon', 'sid', 'k', 't', or 'seq'"
  )
})
