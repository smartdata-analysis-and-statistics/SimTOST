test_that("test_studies executes correctly and returns expected results", {

  # Define dummy input parameters
  nsim <- 5  # Small number of simulations for testing
  n <- 50    # Sample size per arm
  comp <- 1  # Index of comparator

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
    dropout = c(T1 = 0.05, T2 = 0.05),
    alpha = 0.05,
    adjust = "bon",
    k = 1  # Require one endpoint to meet equivalence
  )

  # Set a reproducible seed matrix
  set.seed(123)
  arm_seed <- matrix(sample(1:10000, nsim * 2, replace = TRUE), nrow = nsim, ncol = 2)
  colnames(arm_seed) <- param$arm_names

  # Run test_studies function
  result <- test_studies(nsim, n, comp, param, param.d, arm_seed, ncores = 1)

  # Expected output dimensions
  expect_true(is.matrix(result))   # Should return a matrix
  expect_equal(ncol(result), nsim) # Rows should match number of simulations
  expect_equal(nrow(result), 1 + length(param$list_y_comparator[[comp]]) * 5) # Based on output structure

  # Ensure totaly column contains only 0 or 1
  expect_true(all(result[1, ] %in% c(0, 1)))
})

