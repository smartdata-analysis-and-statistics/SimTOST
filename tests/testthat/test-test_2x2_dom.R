test_that("test_2x2_dom runs correctly with non-hierarchical testing (adseq = FALSE)", {

  # Define input parameters
  n <- 50  # Sample size per arm
  muT <- c(0, 0)  # Mean vector for treatment arm
  muR <- c(log(1.02), log(1.03))   # Mean vector for reference arm
  SigmaW <- diag(c(0.25, 0.3))  # Covariance matrix for treatment
  SigmaB <- 2 * sqrt(max(diag(SigmaW)))
  lequi_tol <- c( log(0.80),  log(0.80))  # Lower equivalence bounds
  uequi_tol <- c( log(1.25),  log(1.25)) # Upper equivalence bounds
  alpha <- c(0.05, 0.05)   # Significance level for each endpoint
  dropout <- c(0.05, 0.05) # 5% dropout per arm
  typey <- -1  # No hierarchical testing
  adseq <- FALSE  # No hierarchical testing
  k <- 2         # Require equivalence for both endpoints
  arm_seed <- 123  # Seed for treatment group
  Eper <- Eco <- c(0, 0)

  # Run the function
  result <- test_2x2_dom(n, muT, muR, SigmaW, lequi_tol, uequi_tol,
                         alpha, SigmaB, dropout, Eper, Eco,
                         typey, adseq, k, arm_seed)

  # Check that result is a matrix
  expect_true(is.matrix(result))

  # Check that the number of columns matches expectations (depends on the function's return structure)
  expect_equal(ncol(result), 11) # totaly[overall], totaly[y0], totaly[y1], mu0[0], mu0[1], mu1[0], mu1[1], sd0[0], sd0[1], sd1[0], sd1[1]

  # Ensure output contains valid numeric values
  expect_true(all(is.finite(result)))

  # Ensure equivalence decision (`totaly`) is binary (0 or 1)
  expect_true(all(result[, 1] %in% c(0, 1)))

  # Ensure no NA values in the result
  expect_false(any(is.na(result)))
})
