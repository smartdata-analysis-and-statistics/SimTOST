test_that("test_par_dom runs correctly with non-hierarchical testing (adseq = FALSE)", {

  # Define input parameters
  n <- 50  # Sample size per arm
  muT <- c(100, 120)  # Mean vector for treatment arm
  muR <- c(98, 118)   # Mean vector for reference arm
  SigmaT <- diag(c(15, 10))  # Covariance matrix for treatment
  SigmaR <- diag(c(15, 10))  # Covariance matrix for reference
  lequi_tol <- c(95, 110)  # Lower equivalence bounds
  uequi_tol <- c(105, 130) # Upper equivalence bounds
  alpha <- c(0.05, 0.05)   # Significance level for each endpoint
  dropout <- c(0.05, 0.05) # 5% dropout per arm
  typey <- -1  # No hierarchical testing
  adseq <- FALSE  # No hierarchical testing
  k <- 2         # Require equivalence for both endpoints
  arm_seedT <- 123  # Seed for treatment group
  arm_seedR <- 456  # Seed for reference group
  TART <- 0.5  # Treatment allocation proportion
  TARR <- 0.5  # Reference allocation proportion
  vareq <- TRUE  # Assume equal variance

  # Run the function
  result <- test_par_dom(n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol,
                         alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR,
                         TART, TARR, vareq)

  # Check that the function runs without error
  expect_silent(test_par_dom(n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol,
                             alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR,
                             TART, TARR, vareq))

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
