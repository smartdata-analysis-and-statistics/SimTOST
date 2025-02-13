test_that("Non-Hierarchical Equivalence Check (adseq = FALSE)", {

  # Define test parameters
  typey <- c(1, 2, 1)  # First and third endpoints are primary, second is secondary
  adseq <- FALSE        # No sequential adjustment (all endpoints evaluated independently)
  tbioq <- matrix(c(1, 0, 1), nrow = 1)  # Equivalence test results (1 = met, 0 = not met)
  k <- 2               # Require at least 2 endpoints to establish equivalence

  # Run function
  result <- check_equivalence(typey, adseq, tbioq, k)

  # Expected behavior:
  # Since `adseq = FALSE`, all endpoints are considered independently.
  # We need at least `k = 2` endpoints to pass equivalence.
  # Endpoints 1 and 3 pass (2 total), so equivalence should be established (`totaly = 1`).

  # Run the test
  expect_equal(result[1, 1], 1, info = "Equivalence should be established.")
})

test_that("Hierarchical Testing Blocks Secondary Endpoints if Primary Fails", {
  typey <- c(1, 2, 1)  # Two primary, one secondary
  adseq <- TRUE
  tbioq <- matrix(c(0, 1, 1), nrow = 1)  # First primary endpoint fails
  k <- 2

  result <- check_equivalence(typey, adseq, tbioq, k)

  expect_equal(result[1, 1], 0, info = "Equivalence should fail since a primary endpoint did not meet criteria.")
})

test_that("Hierarchical Testing All Endpoints Fail", {
  typey <- -1
  adseq <- FALSE
  tbioq <- matrix(c(0, 0, 0), nrow = 1)  # All fail
  k <- 2

  result <- check_equivalence(typey, adseq, tbioq, k)

  expect_equal(result[1, 1], 0, info = "Equivalence should not be established when all endpoints fail.")
})

test_that("Hierarchical Testing Of k Secondary Endpoints", {

  result <- check_equivalence(typey = c(1, 2, 1, 2, 2), adseq = TRUE, tbioq = matrix(c(1, 0, 1, 0, 1), nrow = 1), k = 1)
  expect_equal(result[1, 1], 1, info = "Equivalence should hold when all primary endpoints are equivalent and 1 secondary endpoint is equivalent.")

  result <- check_equivalence(typey = c(1, 2, 1, 2, 2), adseq = TRUE, tbioq = matrix(c(1, 0, 1, 1, 1), nrow = 1), k = 1)
  expect_equal(result[1, 1], 1, info = "Equivalence should hold when all primary endpoints are equivalent and 1 secondary endpoint is equivalent.")

  result <- check_equivalence(typey = c(1, 2, 1, 2, 2), adseq = TRUE, tbioq = matrix(c(1, 0, 1, 0, 1), nrow = 1), k = 2)
  expect_equal(result[1, 1], 0, info = "Equivalence should not hold (not enough secondary endpoints pass).")

  result <- check_equivalence(typey = c(1, 2, 1, 2, 2), adseq = TRUE, tbioq = matrix(c(0, 0, 1, 0, 1), nrow = 1), k = 2)
  expect_equal(result[1, 1], 0, info = "Equivalence should not hold (one of the primary endpoints does not pass).")

})
