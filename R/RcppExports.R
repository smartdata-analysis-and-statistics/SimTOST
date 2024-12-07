# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' ptv
#'
#' Description : to get the p-value with fixed degrees of freedom across the endpoints
#'
#' @param x mat(vector) random variable
#' @param df double degrees of freedom
#' @param lower boolean, calculate the probability of the region where the random variable is less than or equal to x
#' @return mat(vector) with the cumulative distribution function (p-values)
#' @export
ptv <- function(x, df, lower) {
    .Call(`_SimTOST_ptv`, x, df, lower)
}

#' @title Calculate p-values using t-distribution with Variable Degrees of Freedom
#'
#' @description This function computes the cumulative distribution function (p-values) for a given random variable `x`
#' and corresponding degrees of freedom `df` using the t-distribution. The function can compute the lower or upper
#' tail probabilities depending on the value of the `lower` parameter.
#'
#' @param x arma::mat (vector) - A matrix or vector of random variable values for which the p-values will be calculated.
#' @param df arma::mat (vector) - A matrix or vector of degrees of freedom for the t-distribution, matching the size of `x`.
#' @param lower bool - If `TRUE`, calculates the lower-tail probability (P(T <= x)); if `FALSE`, calculates the upper-tail probability.
#' @return arma::mat (vector) - A matrix containing the computed cumulative distribution function (p-values) for each element in `x`.
#' The result is returned as a 1xN matrix, where N is the number of elements in `x`.
#' @export
ptvdf <- function(x, df, lower) {
    .Call(`_SimTOST_ptvdf`, x, df, lower)
}

#' test_2x2_dom_cpp
#'
#' Description : to simulate a 2x2 design and return the p-tost for difference of means (DOM)
#'
#' @param n integer number of subjects per sequence
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaW matrix  within subject covar-variance matrix across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param sigmaB double between subject variance (assumed same for all endpoints)
#' @param dropout vector of size 2 with dropout proportion per sequence (0,1)
#' @param Eper vector of size 2 with period effect on period (0,1)
#' @param Eco vector of size 2 with carry over effect of arm c(Reference, Treatment).
#' @param typey vector with positions of primary endpoints
#' @param adseq boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seed seed for the simulation
#' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
#' @export
test_2x2_dom <- function(n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed) {
    .Call(`_SimTOST_test_2x2_dom`, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed)
}

#' test_2x2_rom_cpp
#'
#' Description : to simulate a 2x2 design and return the p-tost for ratio of means (ROM)
#'
#' @param n integer number of subjects per sequence
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaW matrix  within subject covar-variance matrix across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param sigmaB double between subject variance (assumed same for all endpoints)
#' @param dropout vector of size 2 with dropout proportion per sequence (0,1)
#' @param Eper vector of size 2 with period effect on period (0,1)
#' @param Eco vector of size 2 with carry over effect of arm c(Reference, Treatment).
#' @param typey vector with positions of primary endpoints
#' @param adseq boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seed seed for the simulation
#' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
#' @export
test_2x2_rom <- function(n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed) {
    .Call(`_SimTOST_test_2x2_rom`, n, muT, muR, SigmaW, lequi_tol, uequi_tol, alpha, sigmaB, dropout, Eper, Eco, typey, adseq, k, arm_seed)
}

#' test_par_dom_cpp
#'
#' Description : to simulate a parallel design and return the p-tost for differnece of means (DOM)
#'
#' @param n integer number of subjects per arm
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaT matrix covar-variance matrix on treatment arm across endpoints
#' @param SigmaR matrix covar-variance matrix on reference arm across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param dropout vector of size 2 with dropout proportion per arm (T,R)
#' @param typey vector with positions of primary endpoints
#' @param adseq  boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seedT integer seed for the simulation on treatment arm
#' @param arm_seedR integer seed for the simulation on reference arm
#' @param TART double treatment allocation rate for the treatment arm
#' @param TARR double treatment allocation rate for the reference arm
#' @param vareq boolean assumed equivalence variance between arms for the t-test
#' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
#' @export
test_par_dom <- function(n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq) {
    .Call(`_SimTOST_test_par_dom`, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq)
}

#' test_par_rom_cpp
#'
#' Description : to simulate a parallel design and return the p-tost for ratio of means (ROM)
#'
#' @param n integer number of subjects per arm
#' @param muT vector mean of endpoints on treatment arm
#' @param muR vector mean of endpoints on reference arm
#' @param SigmaT matrix covar-variance matrix on treatment arm across endpoints
#' @param SigmaR matrix covar-variance matrix on reference arm across endpoints
#' @param lequi_tol vector  lower equivalence tolerance band across endpoints
#' @param uequi_tol vector  upper equivalence tolerance band across endpoints
#' @param alpha vector alpha value across endpoints
#' @param dropout vector of size 2 with dropout proportion per arm (T,R)
#' @param typey vector with positions of primary endpoints
#' @param adseq boolean is used a sequential adjustment?
#' @param k integer minimum number of equivalent endpoints
#' @param arm_seedT integer seed for the simulation on treatment arm
#' @param arm_seedR integer seed for the simulation on reference arm
#' @param TART double treatment allocation rate for the treatment arm
#' @param TARR double treatment allocation rate for the reference arm
#' @param vareq boolean assumed equivalence variance between arms for the t-test
#' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
#' @export
test_par_rom <- function(n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq) {
    .Call(`_SimTOST_test_par_rom`, n, muT, muR, SigmaT, SigmaR, lequi_tol, uequi_tol, alpha, dropout, typey, adseq, k, arm_seedT, arm_seedR, TART, TARR, vareq)
}

