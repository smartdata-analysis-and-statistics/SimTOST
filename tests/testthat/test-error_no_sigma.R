#context("Error no sigma")
test_that("error no sigma_list provided", {
  local_edition(3)
  mu_T <- c(37162.0, 35702.0)
  mu_R1 <- c(37705.0, 35930.0, 125.1)
  mu_R2 <- c(37702.8, 35862.4, 126.9)
  mu_list <- list(mu_T,mu_R1,mu_R2)

  sigma_list <- NA

  # Same treatment allocation rate
  TAR = c(1,1,1) # we assume same allocation rate in both arms

  expect_error(SimTOST::sampleSize( mu_list = mu_list,
                                      sigma_list = sigma_list,
                                      varcov_list = NA,
                                      power = 0.9,
                                      dtype = "parallel",
                                      ctype = "ROM",
                                      vareq = T,
                                      lognorm = T,
                                      k=1,
                                      list_comparator =list(c("T","R1"),c("T","R2")),
                                      arm_names=c("T","R_1","R2"),
                                      ncores=1),"No variance-covariance matrix provided, and a standard deviation list is also missing. Either a variance-covariance matrix or a standard deviation list is required.")


})
