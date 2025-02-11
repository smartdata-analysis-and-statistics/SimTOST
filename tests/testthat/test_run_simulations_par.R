test_that("run_simulations_par runs correctly with standard settings", {

  # Define input parameters
  nsim <- 10  # Number of simulations
  n <- 50  # Sample size per arm
  muT <- c(100, 120)  # Mean vector for treatment arm
  muR <- c(98, 118)   # Mean vector for reference arm
  SigmaT <- diag(c(15, 10))  # Covariance matrix for treatment
  SigmaR <- diag(c(15, 10))  # Covariance matrix for reference
  lequi_tol <- c(95, 110)  # Lower equivalence bounds
  uequi_tol <- c(105, 130) # Upper equivalence bounds
  alpha <- c(0.05, 0.05)   # Significance level for each endpoint
  dropout <- c(0.05, 0.05) # 5% dropout per arm
  typey <- c(1, 2)  # First endpoint is primary, second is secondary
  adseq <- FALSE  # No hierarchical testing
  k <- 2         # Require equivalence for both endpoints
  arm_seed_T <- 1:nsim  # Different seeds for each simulation
  arm_seed_R <- (1:nsim) + 100  # Offset reference arm seeds
  ctype <- "DOM"  # Difference of Means test
  lognorm <- FALSE  # Assume normal distribution
  TART <- 0.5  # Treatment allocation proportion
  TARR <- 0.5  # Reference allocation proportion
  vareq <- TRUE  # Assume equal variance

  # Run function
  result <- run_simulations_par(nsim, n, muT, muR, SigmaT, SigmaR,
                                lequi_tol, uequi_tol, alpha, dropout,
                                typey, adseq, k, arm_seed_T, arm_seed_R,
                                ctype, lognorm, TART, TARR, vareq)

  # Check structure of output
  expect_silent(run_simulations_par(nsim, n, muT, muR, SigmaT, SigmaR,
                                    lequi_tol, uequi_tol, alpha, dropout,
                                    typey, adseq, k, arm_seed_T, arm_seed_R,
                                    ctype, lognorm, TART, TARR, vareq))

  # Check if result is a matrix
  expect_true(is.matrix(result))

  # Expected number of rows
  num_endpoints <- length(muT)
  expected_rows <- 1 + num_endpoints * 5  # totaly + 5 cols per endpoint
  expect_equal(nrow(result), expected_rows)

  # Expected number of columns
  expect_equal(ncol(result), nsim)

  # Ensure output contains valid numeric values
  expect_true(all(is.finite(result)))

  # Ensure totaly column (equivalence success) contains only 0 or 1
  expect_true(all(result[1, ] %in% c(0, 1)))

  # Ensure no NA values in the output
  expect_false(any(is.na(result)))
})
