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

  expect_error(SimTOST::sampleSize( mu_list = mu_list,
                                      sigma_list = sigma_list,
                                      varcov_list = NA,
                                      power = 0.9,
                                      dtype = "parallel",
                                      ctype = "ROM",
                                      distribution = "Log Normal",
                                      vareq = T,
                                      k=1,
                                      list_comparator =list(c("T","R1"),c("T","R2")),
                                      arm_names=c("T","R_1","R2"),
                                      ncores=1),"mu_list must be provided")
})

test_that("simulation controls reject invalid values", {
  args <- list(
    n = 10, distribution = "norm",
    mu_list = list(T = c(y1 = 1), R = c(y1 = 1)),
    sigma_list = list(T = c(y1 = .2), R = c(y1 = .2)),
    list_comparator = list(T_vs_R = c("T", "R")),
    list_lequi.tol = list(T_vs_R = c(y1 = -.5)),
    list_uequi.tol = list(T_vs_R = c(y1 = .5)), nsim = 2
  )
  expect_error(do.call(simPower, c(args, alpha = 0)), "alpha")
  expect_error(do.call(simPower, c(args, dropout = 1)), "dropout")
  expect_error(do.call(simPower, c(args, rho = 1.1)), "rho")
})
