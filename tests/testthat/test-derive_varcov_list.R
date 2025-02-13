#context("Error lenght mu varcov tar")
test_that("differences across lenght of mu, varcov and tar", {

  mu_T <- c(37162.0, 35702.0, 125.1)
  mu_R1 <- c(37705.0, 35930.0, 125.1)
  mu_R2 <- c(37702.8, 35862.4, 126.9)
  mu_list <- list(mu_T,mu_R1,mu_R2)
  # Save variance covariance matrix in a list
  sigma_T <- c(11113.62172, 9132.75342, 16.89586)
  sigma_R1 <- c(12332.41615, 9398.42182, 17.88151)
  sigma_R2 <- c(12113.72, 9098.42182, 17.1586)

  sigma_list <- list(sigma_T,sigma_R1,sigma_R2)

  expect_type(derive_varcov_list(mu_list = mu_list, sigma_list = sigma_list), "list")
  expect_error(derive_varcov_list(mu_list = mu_list[[-1]], sigma_list = sigma_list))
})
