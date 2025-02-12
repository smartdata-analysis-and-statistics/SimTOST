test_that("power_cal produces valid results for parallel DOM design", {

  # Define dummy input parameters
  nsim <- 5  # Small number of simulations for testing
  n <- 50    # Sample size per arm

  # Create mock parameter lists
  param <- list(
    mu = list(
      SB2 = matrix(c(38703, 36862, 127), nrow = 1, dimnames = list(NULL, c("AUCinf", "AUClast", "Cmax"))),
      EUREF = matrix(c(39360, 37022, 126.2), nrow = 1, dimnames = list(NULL, c("AUCinf", "AUClast", "Cmax")))
    ),
    varcov = list(
      SB2 = matrix(c(123520996, 60902497.20, 112695.96,
                     60902497.20, 83411689.00,  92608.62,
                     112695.96,   92608.62,    285.61),
                   nrow = 3, byrow = TRUE, dimnames = list(c("AUCinf", "AUClast", "Cmax"), c("AUCinf", "AUClast", "Cmax"))),
      EUREF = matrix(c(152078224.0, 69537681.6, 132445.68,
                       69537681.6,  88322404.0, 100934.52,
                       132445.68,   100934.52,    320.41),
                     nrow = 3, byrow = TRUE, dimnames = list(c("AUCinf", "AUClast", "Cmax"), c("AUCinf", "AUClast", "Cmax")))
    ),
    sigmaB = NA,
    TAR_list = list(SB2 = 1, EUREF = 1),
    type_y = c(AUCinf = 1, AUClast = 1, Cmax = 1),
    weight_seq = c(AUCinf = 1/3, AUClast = 1/3, Cmax = 1/3),
    arm_names = c("SB2", "EUREF"),
    ynames_list = list(SB2 = c("AUCinf", "AUClast", "Cmax"),
                       EUREF = c("AUCinf", "AUClast", "Cmax")),
    list_comparator = list(EMA = c("SB2", "EUREF")),
    list_y_comparator = list(EMA = c("AUCinf", "AUClast", "Cmax")),
    list_lequi.tol = list(EMA = c(AUCinf = 0.8, AUClast = 0.8, Cmax = 0.8)),
    list_uequi.tol = list(EMA = c(AUCinf = 1.25, AUClast = 1.25, Cmax = 1.25)),
    Eper = c(0, 0),
    Eco = c(0, 0)
  )

  param.d <- list(
    dtype = "parallel",
    ctype = "ROM",
    lognorm = TRUE,
    vareq = TRUE,
    list_lequi.tol = list(EMA = c(AUCinf = 0.8, AUClast = 0.8, Cmax = 0.8)),
    list_uequi.tol = list(EMA = c(AUCinf = 1.25, AUClast = 1.25, Cmax = 1.25)),
    dropout = c("SB2" = 0.05, "EUREF" = 0.05),
    alpha = 0.05,
    adjust = "bon",
    k = 1  # Require one endpoint to meet equivalence
  )

  # Run the function
  result <- power_cal(n = n, nsim = nsim, param = param, param.d = param.d, seed = 1234, ncores = 1)

  # Check structure
  expect_type(result, "list")
  expect_true("power" %in% names(result), "Output should contain 'power'")
  expect_true("output.test" %in% names(result), "Output should contain 'output.test'")

  # Validate power
  expect_type(result$power, "double")
  expect_true(result$power >= 0 && result$power <= 1, "Power should be between 0 and 1")

  # Validate output.test
  expect_true(is.data.frame(result$output.test), "Output test should be a data frame")
  expect_gt(nrow(result$output.test), 0, "Output test should have rows")
  expect_gt(ncol(result$output.test), 0, "Output test should have columns")

  # Check if power is computed correctly
  expect_equal(result$power, sum(result$output.test[,1]) / nrow(result$output.test), tolerance = 1e-5)

  # Check for NA values in output.test
  expect_false(any(is.na(result$output.test)), "There should be no NA values in output.test")
})

test_that("power_cal produces valid results for 2x2 cross-over DOM design", {

  # Define dummy input parameters
  nsim <- 10  # Small number of simulations for testing
  n <- 50    # Sample size per arm

  param <- list(
    mu = list(
      "R" = matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("AUC", "Cmax"))),
      "T" = matrix(c(log(1.02), log(1.03)), nrow = 1, dimnames = list(NULL, c("AUC", "Cmax")))
    ),
    varcov = list(
      "R" = matrix(c(0.06250, 0.01875, 0.01875, 0.09000), nrow = 2, byrow = TRUE,
                   dimnames = list(c("AUC", "Cmax"), c("AUC", "Cmax"))),
      "T" = matrix(c(0.06250, 0.01875, 0.01875, 0.09000), nrow = 2, byrow = TRUE,
                   dimnames = list(c("AUC", "Cmax"), c("AUC", "Cmax")))
    ),
    TAR_list = list(
      "R" = 1,
      "T" = 1
    ),
    type_y = -1,
    weight_seq = c(AUC = 0.5, Cmax = 0.5),
    arm_names = c("R", "T"),
    ynames_list = list(
      R = c("AUC", "Cmax"),
      T = c("AUC", "Cmax")
    ),
    list_comparator = list(
      T_vs_R = c("R", "T")
    ),
    list_y_comparator = list(
      T_vs_R = c("AUC", "Cmax")
    ),
    list_lequi.tol = list(
      T_vs_R = c(AUC = -0.2231436, Cmax = -0.2231436)
    ),
    list_uequi.tol = list(
      T_vs_R = c(AUC = 0.2231436, Cmax = 0.2231436)
    ),
    sigmaB = NA,
    Eper = c(0, 0),
    Eco = c(0, 0)
  )

  param.d <- list(
    power = 0.8,
    alpha = 0.05,
    dtype = "2x2",
    ctype = "DOM",
    lognorm = FALSE,
    vareq = TRUE,
    k = 2,
    adjust = "no",
    dropout = c(0, 0),

    list_lequi.tol = list(
      T_vs_R = c(AUC = -0.2231436, Cmax = -0.2231436)
    ),

    list_uequi.tol = list(
      T_vs_R = c(AUC = 0.2231436, Cmax = 0.2231436)
    )
  )

  # Run the function
  result <- power_cal(n = n, nsim = nsim, param = param, param.d = param.d, seed = 1234, ncores = 1)

  # Check structure
  expect_type(result, "list")
  expect_true("power" %in% names(result), "Output should contain 'power'")
  expect_true("output.test" %in% names(result), "Output should contain 'output.test'")

  # Validate power
  expect_type(result$power, "double")
  expect_true(result$power >= 0 && result$power <= 1, "Power should be between 0 and 1")

  # Validate output.test
  expect_true(is.data.frame(result$output.test), "Output test should be a data frame")
  expect_gt(nrow(result$output.test), 0, "Output test should have rows")
  expect_gt(ncol(result$output.test), 0, "Output test should have columns")

  # Check if power is computed correctly
  expect_equal(result$power, sum(result$output.test[,1]) / nrow(result$output.test), tolerance = 1e-5)

  # Check for NA values in output.test
  expect_false(any(is.na(result$output.test)), "There should be no NA values in output.test")



})
