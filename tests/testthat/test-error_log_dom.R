# #context("Error log dom")
# test_that("error when lognorm =True and ctype=DOM", {
#   local_edition(3)
#   mu_T <- c(37162.0, 35702.0, 125.1)
#   mu_R1 <- c(37705.0, 35930.0, 125.1)
#   mu_R2 <- c(37702.8, 35862.4, 126.9)
#   mu_list <- list(mu_T,mu_R1,mu_R2)
#   # Save variance covariance matrix in a list
#   sigma_T <- c(11113.62172, 9132.75342, 16.89586)
#   sigma_R1 <- c(12332.41615, 9398.42182, 17.88151)
#   sigma_R2 <- c(12113.72, 9098.42182, 17.1586)
#
#   sigma_list <- list(sigma_T,sigma_R1,sigma_R2)
#
#   # Same treatment allocation rate
#   TAR = c(1,1,1) # we assume same allocation rate in both arms
#
#   expect_error(SimTOST::sampleSize( mu_list = mu_list,
#                                       sigma_list = sigma_list,
#                                       varcov_list = NA,
#                                       TAR=TAR,
#                                       power = 0.9,
#                                       dtype = "parallel",
#                                       ctype = "DOM",
#                                       vareq = T,
#                                       lognorm = T,
#                                       k=1,
#                                       list_comparator =list(c("T","R1"),c("T","R2")),
#                                       arm_names=c("T","R_1","R2"),
#                                       ncores=1),"No test available for DOM when variables are log normal distributed")
# })
