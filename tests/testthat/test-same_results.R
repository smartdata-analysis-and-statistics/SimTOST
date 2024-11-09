#context("Same result")
testthat::test_that("Return expected result", {
  local_edition(3)
  mu_T <- c(AUCinf = 37162.0, AUClast = 35702.0, Cmax = 125.9)
  mu_R1 <- c(AUCinf = 37705.0, AUClast = 35930.0, Cmax = 125.1)
  mu_R2 <- c(AUCinf = 37702.8, AUClast = 35862.4, Cmax = 126.9)
  mu_list <- list("SB2" = mu_T,
                  "EUREF" = mu_R1,
                  "USREF" = mu_R2)

  # Save variance covariance matrix in a list
  sigma_T <- c(AUCinf = 11113.62172, AUClast = 9132.75342, Cmax = 16.89586)
  sigma_R1 <- c(AUCinf = 12332.41615, AUClast = 9398.42182, Cmax = 17.88151)
  sigma_R2 <- c(AUCinf = 12113.72, AUClast = 9098.42182, Cmax = 17.1586)

  sigma_list <- list("SB2" = sigma_T,
                     "EUREF" = sigma_R1,
                     "USREF" = sigma_R2)

  # Same treatment allocation rate
  TAR = c(1,1,1) # we assume same allocation rate in both arms

  # arms to be compared
  list_comparator <- list(EMA = c("SB2", "EUREF"),
                          FDA = c("SB2", "USREF"))

  # endpoint to be compared
  list_y_comparator <- list(EMA = c("AUCinf", "AUClast", "Cmax"),
                            FDA = c("AUCinf", "AUClast", "Cmax"))

  # Equivalence boundaries
  lequi.tol <- c(AUCinf = 0.8, AUClast = 0.8, Cmax = 0.8)
  uequi.tol <- c(AUCinf = 1.25, AUClast = 1.25, Cmax = 1.25)

  # Pass the user parameters into a list of parameters and calculate the sample size
  res_cal <- sampleSize(mu_list = mu_list, sigma_list = sigma_list,
                        power = 0.9, dtype = "parallel", ctype = "ROM",
                        vareq = T, lognorm = T, k = 3,
                        list_comparator = list_comparator,
                        list_y_comparator = list_y_comparator,
                        lequi.tol = lequi.tol, uequi.tol = uequi.tol,
                        ncores = 1)

  expect_equal(res_cal$response[["n_total"]], 150, tolerance = 6)
})
