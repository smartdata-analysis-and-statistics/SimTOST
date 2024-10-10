#context("Error no mu")
test_that("error no mu_list provided", {
  local_edition(3)
  # Save variance covariance matrix in a list
  sigma_T <- c(11113.62172, 9132.75342, 16.89586)
  sigma_R1 <- c(12332.41615, 9398.42182, 17.88151)
  sigma_R2 <- c(12113.72, 9098.42182, 17.1586)

  sigma_list <- list(sigma_T,sigma_R1,sigma_R2)

  mu_list<-NA
  # Same treatment allocation rate
  TAR = c(1,1,1) # we assume same allocation rate in both arms

  expect_error(simsamplesize::calopt( mu_list = mu_list,
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
                                      ncores=1),"mu_list must be provided")
})
