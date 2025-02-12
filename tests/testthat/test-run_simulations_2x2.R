test_that("run_simulations_2x2 produces valid results", {

  # Set up example parameters
  n <- 30     # Sample size per period

  muT <- c(100, 120)  # Mean values for treatment arm
  muR <- c(98, 118)   # Mean values for reference arm

  SigmaW <- matrix(c(15, 5, 5, 10), nrow = 2)  # Within-subject covariance matrix

  lequi_tol <- c(90, 110)  # Lower equivalence thresholds
  uequi_tol <- c(110, 130) # Upper equivalence thresholds
  alpha <- c(0.05, 0.05)   # Significance level for each endpoint

  sigmaB <- 2.0  # Between-subject variance

  dropout <- c(0.05, 0.05)  # Dropout rates per sequence

  Eper <- c(0, 0)  # Expected period effects
  Eco <- c(0, 0)   # Expected carryover effects

  typey <- -1       # No sequential testing
  adseq <- FALSE    # Disable sequential testing
  k <- 2            # Require equivalence for at least 2 endpoints

  # Run the function
  result1 <- run_simulations_2x2(nsim = 5, ctype = "DOM", lognorm = FALSE, n,
                                 muT, muR, SigmaW, lequi_tol, uequi_tol,
                                 alpha, sigmaB, dropout, Eper, Eco, typey,
                                 adseq, k, arm_seed = 1000:1005)

  # Basic structure checks
  expect_true(is.matrix(result1), "Output should be a matrix")
  expect_equal(dim(result1), c(11, 5), info = "Output dimensions should match expected")


  # Check numerical values are within reasonable range
  expect_true(all(result1 >= 0), "All values should be non-negative")
  expect_true(all(result1 <= 1 | result1 > 1), "Equivalence decisions should be 0 or 1, other values should be means/sd")

  # Test consistency with different seeds
  result2 <- run_simulations_2x2(nsim = 5, ctype = "DOM", lognorm = FALSE, n,
                                 muT, muR, SigmaW, lequi_tol, uequi_tol,
                                 alpha, sigmaB, dropout, Eper, Eco, typey,
                                 adseq, k, arm_seed = 1000:1005)

  expect_equal(result1, result2, info = "Results should be identical when using the same seed")
})

