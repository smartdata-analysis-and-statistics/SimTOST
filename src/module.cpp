#include <iostream>
#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

using namespace Rcpp;
using namespace arma;

//' @title Compute p-values for a t-distribution with Fixed Degrees of Freedom
//'
//' @description Computes p-values for a given set of random variables under a t-distribution with fixed degrees of freedom.
//'
//' @param x A numeric matrix (or vector) representing the random variables.
//' @param df A double specifying the degrees of freedom.
//' @param lower A logical value indicating whether to compute the lower-tail probability (\code{P(T <= x)}). If \code{FALSE}, the function returns the upper-tail probability (\code{P(T > x)}).
//'
//' @return A numeric matrix containing the computed cumulative distribution function (CDF) values (p-values).
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::mat ptv(arma::mat x, const double df, const bool lower) {
  Rcpp::NumericVector x_rcpp= as<NumericVector>(wrap(x));
  int n = x_rcpp.size();
  Rcpp::NumericVector y(n);
  y = Rcpp::pt(x_rcpp, df, lower, false);
  arma::vec y_arma = as<arma::vec>(wrap(y));
  mat y_mat = reshape(y_arma, 1, n);
  return y_mat;
}

//' @title Calculate p-values using t-distribution with Variable Degrees of Freedom
//'
//' @description This function computes the cumulative distribution function (p-values) for a given random variable \code{x} and corresponding degrees of freedom \code{df} using the t-distribution. The function can compute the lower or upper tail probabilities depending on the value of the \code{lower} argument.
//'
//' @param x arma::mat (vector) - A matrix or vector of random variable values for which the p-values will be calculated.
//' @param df arma::mat (vector) - A matrix or vector of degrees of freedom for the t-distribution, matching the size of \code{x}.
//' @param lower bool - If \code{TRUE}, calculates the lower-tail probability (P(T <= x)); if \code{FALSE}, calculates the upper-tail probability.
//'
//' @return arma::mat (vector) - A matrix containing the computed cumulative distribution function (p-values) for each element in \code{x}. The result is returned as a 1xN matrix, where N is the number of elements in \code{x}.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::mat ptvdf(arma::mat x, arma::mat df, const bool lower) {
  std::size_t n = x.n_elem;  // Get the number of elements in x
  arma::vec y(n);  // Create an Armadillo vector for the result

  for (std::size_t i = 0; i < n; ++i) {
    // Directly pass individual elements of x and df to Rcpp::pt()
    y[i] = R::pt(x(i), df(i), lower, false);  // Use R's pt function directly
    }

  // Reshape y into a 1xN matrix and return
  return arma::reshape(y, 1, n);
}

//' @title Check Equivalence for Multiple Endpoints
//'
//' @description This function evaluates whether equivalence criteria are met based on a predefined set of endpoints. It first checks whether all primary endpoints satisfy equivalence (if sequential testing is enabled). Then, it determines whether the required number of endpoints (\code{k}) meet the equivalence threshold. The function returns a binary decision indicating whether overall equivalence is established.
//'
//' @param typey An integer vector specifying the hierarchy of each endpoint, where \code{1} denotes a primary endpoint and \code{2} denotes a secondary endpoint.
//' @param adseq A boolean flag indicating whether sequential testing is enabled. If set to \code{TRUE}, all primary endpoints must pass equivalence before secondary endpoints are evaluated. If set to \code{FALSE}, primary and secondary endpoints are assessed independently.
//' @param tbioq A matrix containing the equivalence test results for each endpoint, where \code{1} indicates that equivalence is met and \code{0} indicates that equivalence is not met.
//' @param k An integer specifying the minimum number of endpoints required for overall equivalence.
//'
//' @details When sequential testing is enabled (\code{adseq = TRUE}), all primary endpoints must meet equivalence before secondary endpoints are considered. If sequential testing is disabled (\code{adseq = FALSE}), all endpoints are evaluated simultaneously without hierarchical constraints. The function then determines whether at least \code{k} endpoints meet the equivalence criteria. If the conditions are satisfied, the final equivalence decision (\code{totaly}) is \code{1}; otherwise, it is \code{0}.
//'
//' @return Returns a (1 × 1 matrix) containing a binary equivalence decision. A value of \code{1} indicates that equivalence is established, while \code{0} indicates that equivalence is not established.
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::mat check_equivalence(const arma::uvec& typey,
                            const bool adseq,
                            const arma::mat& tbioq,
                            const int k) {

  // Initialize final equivalence decision
  arma::mat totaly(1, 1, arma::fill::zeros);

  // If no sequential testing is required, evaluate equivalence criteria directly
  if (!adseq || typey.empty() || any(typey < 0)) {
    if (k < 0) {
      totaly(0, 0) = (accu(tbioq) == tbioq.n_cols) ? 1 : 0; // All endpoints must pass
    } else {
      totaly(0, 0) = (accu(tbioq) >= k) ? 1 : 0; // At least k endpoints must pass
    }
    return totaly;
  }

  // Count the number of primary endpoints
  int num_primary_endpoints = accu(typey == 1);

  // If no primary endpoints exist, assume primary endpoint test is automatically passed
  bool primary_test_passed = (num_primary_endpoints == 0);
  int num_primary_successes = arma::accu(tbioq.cols(arma::find(typey == 1)));

  // Ensure all primary endpoints pass if sequential testing is enabled
  if (num_primary_endpoints > 0) {
    primary_test_passed = (num_primary_successes == num_primary_endpoints);
  }

  // Ensure at least k secondary endpoints pass
  int num_secondary_successes = arma::accu(tbioq.cols(arma::find(typey == 2)));
  bool secondary_test_passed = (num_secondary_successes >= k);

  // Final decision
  totaly(0, 0) = (secondary_test_passed && primary_test_passed) ? 1 : 0;

  return totaly;
}


//' @title Simulate a 2x2 Crossover Design and Compute Difference of Means (DOM)
//'
//' @description Simulates a two-sequence, two-period (2x2) crossover design and evaluate equivalence for the difference of means (DOM).
//'
//' @param n integer number of subjects per sequence
//' @param muT vector mean of endpoints on treatment arm
//' @param muR vector mean of endpoints on reference arm
//' @param SigmaW matrix  within subject covar-variance matrix across endpoints
//' @param lequi_tol vector  lower equivalence tolerance band across endpoints
//' @param uequi_tol vector  upper equivalence tolerance band across endpoints
//' @param alpha vector alpha value across endpoints
//' @param sigmaB double between subject variance (assumed same for all endpoints)
//' @param dropout vector of size 2 with dropout proportion per sequence (0,1)
//' @param Eper vector of size 2 with period effect on period (0,1)
//' @param Eco vector of size 2 with carry over effect of arm c(Reference, Treatment).
//' @param typey vector with positions of primary endpoints
//' @param adseq boolean is used a sequential adjustment?
//' @param k integer minimum number of equivalent endpoints
//' @param arm_seed seed for the simulation
//'
//' @return A numeric matrix containing the simulated hypothesis test results.
//' The first column represents the overall equivalence decision, where 1 indicates
//' success and 0 indicates failure. The subsequent columns contain the hypothesis
//' test results for each endpoint, followed by mean estimates for the reference and
//' treatment groups, and standard deviations for the reference and treatment groups.
//'
//' @export
// [[Rcpp::export]]
arma::mat test_2x2_dom(const int n,
                       const arma::vec& muT,
                       const arma::vec& muR,
                       const arma::mat& SigmaW,
                       const arma::rowvec& lequi_tol,
                       const arma::rowvec& uequi_tol,
                       const arma::rowvec& alpha,
                       const double sigmaB,
                       const arma::vec& dropout,
                       const arma::vec& Eper,
                       const arma::vec& Eco,
                       const arma::uvec& typey,
                       const bool adseq,
                       const int k,
                       const int arm_seed){

  // Set random seed
  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(arm_seed);

  // Transform drop out
   int n0i = ceil(n/2);
   int n1i = n - n0i;
   int n0 = ceil((1 - dropout[0])* n0i);
   if (n0 < 2) n0 = 2;

   int n1 = ceil((1 - dropout[1])*n1i);
   if (n1 < 2) n1 = 2;

   int nt = n0 + n1;

   ivec pid = rep_each(seq_len(nt), 2);
   vec pid_u = rep_each(rnorm(nt, 0, sigmaB), 2);
   ivec seq0 = rep(0, n0*2);
   ivec seq1 = rep(1, n1*2);
   ivec seq = join_cols(seq0,seq1);
   ivec per = rep(IntegerVector::create(0, 1), nt );


   // Generate data based on the provided covariance matrix and means
   mat result(nt*2, muT.size());
   mat diff(nt, muT.size());
   ivec trt(nt * 2);
   ivec s0p0(nt * 2);
   vec mean_val;
   int j = 0;

   for (int i = 0; i < nt * 2; ++i) {
     trt[i] = (seq[i] == 1) ? 1 - per[i] : per[i];
     if (seq[i] == 0 && per[i] == 0) {
       mean_val = muR + Eper[0] + pid_u[i];
       s0p0[i] = 1;
     } else if (seq[i] == 1 && per[i] == 0) {
       mean_val = muT + Eper[0] + pid_u[i];
     } else if (seq[i] == 0 && per[i] == 1) {
       mean_val = muT + Eper[1] + Eco[0] + pid_u[i];
     } else {//seq ==1 per==1
       mean_val = muR + Eper[1] + Eco[1] + pid_u[i];
     }
     result.row(i) = arma::mvnrnd(mean_val,SigmaW,1).t();
     if ( per[i] == 1){
       diff.row(j) = result.row(i)-result.row(i-1);
       if (seq[i] == 0){
         diff.row(j) = - diff.row(j);
       }
       j++;
     }
   }

   uvec trt0 = find(trt == 0); // filter reference observations
   uvec trt1 = find(trt == 1); // filter treatment observations
   uvec fs0p0 = find(s0p0 == 1); // filter period ==0 sequence == 0

   mat mut0 = mean(result.rows(trt0),0);
   mat mut1 = mean(result.rows(trt1),0);
   mat sdt0 = stddev(result.rows(trt0),0,0);
   mat sdt1 = stddev(result.rows(trt1),0,0);
   mat sdb = stddev(result.rows(fs0p0),0,0);
   mat cor01 = cor(result.rows(trt0),result.rows(trt1));
   mat sdw = stddev(diff,0,0)/sqrt(2.0);
   mat sde = sdw * sqrt(2.0/nt);
   mat tlb = (mut0 - mut1 - lequi_tol) /sde;
   mat tub = (mut0 - mut1 - uequi_tol) /sde;

   // Calculate p-value
   double df = 1.0*(nt - 2);
   mat plb = ptv(tlb,df,false);
   mat pub = ptv(tub,df,true);
   mat ptost = max(plb, pub);

   ptost.replace(datum::nan, 0); // in case of NA values due to no sd in some studies
   mat alpha0 = conv_to<mat>::from(alpha);
   mat tbioq  = conv_to<mat>::from((ptost < alpha0));

   // Call the check_equivalence function to determine if equivalence is established
   arma::mat totaly = check_equivalence(typey, adseq, tbioq, k);

   mat response0 = join_rows<mat>(totaly,tbioq);
   mat response1 = join_rows<mat>(mut0,mut1);
   mat response2 = join_rows<mat>(sdw,sdb);
   mat response3 = join_rows<mat>(response0,response1);

   return join_rows<mat>(response3,response2);
}

//' @title Simulate a 2x2 Crossover Design and Compute Ratio of Means (ROM)
//'
//' @description Simulates a two-sequence, two-period (2x2) crossover design and evaluate equivalence for the ratio of means (ROM).
//'
//' @param n integer number of subjects per sequence
//' @param muT vector mean of endpoints on treatment arm
//' @param muR vector mean of endpoints on reference arm
//' @param SigmaW matrix  within subject covar-variance matrix across endpoints
//' @param lequi_tol vector  lower equivalence tolerance band across endpoints
//' @param uequi_tol vector  upper equivalence tolerance band across endpoints
//' @param alpha vector alpha value across endpoints
//' @param sigmaB double between subject variance (assumed same for all endpoints)
//' @param dropout vector of size 2 with dropout proportion per sequence (0,1)
//' @param Eper vector of size 2 with period effect on period (0,1)
//' @param Eco vector of size 2 with carry over effect of arm c(Reference, Treatment).
//' @param typey vector with positions of primary endpoints
//' @param adseq boolean is used a sequential adjustment?
//' @param k integer minimum number of equivalent endpoints
//' @param arm_seed seed for the simulation
//'
//' @return A numeric matrix containing the simulated hypothesis test results. The first column represents the overall equivalence decision, where 1 indicates success and 0 indicates failure. The subsequent columns contain the hypothesis test results for each endpoint, followed by mean estimates for the reference and treatment groups, and standard deviations for the reference and treatment groups.
//'
//' @export
// [[Rcpp::export]]
arma::mat test_2x2_rom(const int n,
                       const arma::vec& muT,
                       const arma::vec& muR,
                       const arma::mat& SigmaW,
                       const arma::rowvec& lequi_tol,
                       const arma::rowvec& uequi_tol,
                       const arma::rowvec& alpha,
                       const double sigmaB,
                       const arma::vec& dropout,
                       const arma::vec& Eper,
                       const arma::vec& Eco,
                       const arma::uvec& typey,
                       const bool adseq,
                       const int k,
                       const int arm_seed){

  // Set random seed
  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(arm_seed);


  // Power test based on Hauschke et al, 1999
   // Transform drop out

   int n0i = ceil(n/2);
   int n0 = ceil((1 - dropout[0])* n0i);
   if (n0 < 2) n0 = 2;

   int n1 = ceil((1 - dropout[1])* (n-n0i));
   if (n1 < 2) n1 = 2;

   int nt = n0 + n1;

   ivec pid = rep_each(seq_len(nt),2);
   vec pid_u = rep_each(rnorm(nt, 0, sigmaB),2);
   ivec seq0 = rep(0, n0*2);
   ivec seq1 = rep(1, n1*2);
   ivec seq = join_cols(seq0,seq1);
   ivec per = rep(IntegerVector::create(0, 1), nt );

   // Generate data based on the provided covariance matrix and means
   mat result(nt*2, muT.size());
   ivec trt(nt * 2);
   ivec s1p0(nt * 2);
   ivec s1p1(nt * 2);
   ivec s0p1(nt * 2);
   ivec s0p0(nt * 2);
   vec mean_val;


   for (int i = 0; i < nt * 2; ++i) {
     trt[i] = (seq[i] == 1) ? 1 - per[i] : per[i];
     if (seq[i] == 0 && per[i] == 0) {
       mean_val = muR + Eper[0] + pid_u[i];
       s0p0[i] = 1;
     } else if (seq[i] == 1 && per[i] == 0) {
       mean_val = muT + Eper[0] + pid_u[i];
       s1p0[i] = 1;
     } else if (seq[i] == 0 && per[i] == 1) {
       mean_val = muT + Eper[1] + Eco[0] + pid_u[i];
       s0p1[i] = 1;
     } else {//seq ==1 per==1
       mean_val = muR + Eper[1] + Eco[1] + pid_u[i];
       s1p1[i] = 1;
     }
     result.row(i) = arma::mvnrnd(mean_val,SigmaW,1).t();
   }

   uvec fs0p0 = find(s0p0 == 1); // filter sequence == 0 period 0
   uvec fs0p1 = find(s0p1 == 1);
   uvec fs1p0 = find(s1p0 == 1);
   uvec fs1p1 = find(s1p1 == 1);

   mat mus0p0 = mean(result.rows(fs0p0),0);
   mat mus0p1 = mean(result.rows(fs0p1),0);
   mat mus1p0 = mean(result.rows(fs1p0),0);
   mat mus1p1 = mean(result.rows(fs1p1),0);

   mat vars0p0 = var(result.rows(fs0p0),0,0);
   mat vars0p1 = var(result.rows(fs0p1),0,0);
   mat vars1p0 = var(result.rows(fs1p0),0,0);
   mat vars1p1 = var(result.rows(fs1p1),0,0);

   mat covs0 = cov(result.rows(fs0p0),result.rows(fs0p1));
   mat covs1 = cov(result.rows(fs1p0),result.rows(fs1p1));

   mat mut0 = (mus0p0 + mus1p1)/2.0;   // mean reference
   mat mut1 = (mus0p1 + mus1p0)/2.0;    // mean treatment
   mat vart0 = (vars0p0*(n0-1.0) + vars1p1*(n1-1.0))/(n0 + n1 - 2.0); // var reference;
   mat vart1 = (vars0p1*(n0-1.0) + vars1p0*(n1-1.0))/(n0 + n1 - 2.0); // var treatment;
   mat covt0t1 = (covs0*(n0-1.0) + covs1*(n1-1.0))/(n0 + n1 - 2.0); // cov reference-treatment;

   mat sdlb = sqrt(vart0 - 2.0*lequi_tol%covt0t1 + arma::pow(lequi_tol, 2)%vart1);
   mat sdub = sqrt(vart0 - 2.0*uequi_tol%covt0t1 + arma::pow(uequi_tol, 2)%vart1);

   mat tlb = (mut0 - lequi_tol%mut1) / (sdlb/2.0*sqrt(1.0/n0 + 1.0/n1));
   mat tub = (mut0 - uequi_tol%mut1) / (sdub/2.0*sqrt(1.0/n0 + 1.0/n1));

   // Calculate p-value
   double df = 1.0*(n0 + n1 - 2);
   mat plb = ptv(tlb,df,false);
   mat pub = ptv(tub,df,true);
   mat ptost = max(plb, pub);

   ptost.replace(datum::nan, 0); // in case of NA values due to no sd in some studies
   mat alpha0 = conv_to<mat>::from(alpha);
   mat tbioq  = conv_to<mat>::from((ptost < alpha0));

   // Call the check_equivalence function to determine if equivalence is established
   arma::mat totaly = check_equivalence(typey, adseq, tbioq, k);

   mat response0 = join_rows<mat>(totaly,tbioq);
   mat response1 = join_rows<mat>(mut0,mut1);
   mat response2 = join_rows<mat>(vart0,vart1);
   mat response3 = join_rows<mat>(response2,covt0t1);
   mat response4 = join_rows<mat>(response0,response1);

   return join_rows<mat>(response4,response3);
}

//' @title Simulate a Parallel Design and Test Difference of Means (DOM)
//'
//' @description
//' Simulates a parallel-group design and performs equivalence testing using the difference of means (DOM) approach.
//' This function evaluates whether the treatment and reference groups are equivalent based on predefined
//' equivalence margins and hypothesis testing criteria.
//'
//' @param n integer number of subjects per arm
//' @param muT vector mean of endpoints on treatment arm
//' @param muR vector mean of endpoints on reference arm
//' @param SigmaT matrix covar-variance matrix on treatment arm across endpoints
//' @param SigmaR matrix covar-variance matrix on reference arm across endpoints
//' @param lequi_tol vector  lower equivalence tolerance band across endpoints
//' @param uequi_tol vector  upper equivalence tolerance band across endpoints
//' @param alpha vector alpha value across endpoints
//' @param dropout vector of size 2 with dropout proportion per arm (T,R)
//' @param typey vector with positions of primary endpoints
//' @param adseq  boolean is used a sequential adjustment?
//' @param k integer minimum number of equivalent endpoints
//' @param arm_seedT integer seed for the simulation on treatment arm
//' @param arm_seedR integer seed for the simulation on reference arm
//' @param TART double treatment allocation rate for the treatment arm
//' @param TARR double treatment allocation rate for the reference arm
//' @param vareq boolean assumed equivalence variance between arms for the t-test
//'
//' @return A numeric matrix containing the simulated hypothesis test results.
//' The first column represents the overall equivalence decision, where 1 indicates
//' success and 0 indicates failure. The subsequent columns contain the hypothesis
//' test results for each endpoint, followed by mean estimates for the reference and
//' treatment groups, and standard deviations for the reference and treatment groups.
//'
//' @details
//' The function simulates a parallel-group study design and evaluates equivalence
//' using the difference of means (DOM) approach. It accounts for dropout rates and
//' treatment allocation proportions while generating simulated data based on the
//' specified covariance structure. The test statistics are computed, and a final
//' equivalence decision is made based on the predefined number of required significant
//' endpoints (\code{k}). If sequential testing (\code{adseq}) is enabled, primary endpoints
//' must establish equivalence before secondary endpoints are evaluated.
//' When \code{vareq = TRUE}, the test assumes equal variances between groups and
//' applies Schuirmann's two one-sided tests (TOST).
//'
//' @export
// [[Rcpp::export]]
arma::mat test_par_dom(const int n,
                       const arma::vec& muT,
                       const arma::vec& muR,
                       const arma::mat& SigmaT,
                       const arma::mat& SigmaR,
                       const arma::rowvec& lequi_tol,
                       const arma::rowvec& uequi_tol,
                       const arma::rowvec& alpha,
                       const arma::vec& dropout,
                       const arma::uvec& typey,
                       const bool adseq,
                       const int k,
                       const int arm_seedT,
                       const int arm_seedR,
                       const double TART,
                       const double TARR,
                       const bool vareq){

  // Transform drop out
   int n0i = ceil(n*TART);
   int n1i = ceil(n*TARR);

   int n0 = ceil((1 - dropout[0])* n0i);
   if (n0 < 2) n0 = 2;

   int n1 = ceil((1 - dropout[1])*n1i);
   if (n1 < 2) n1 = 2;


   RNGScope scope;
   Environment base_env("package:base");
   Function set_seed = base_env["set.seed"];

   // Generate data based on the provided covariance matrix and means
   set_seed(arm_seedT);
   mat yT = arma::mvnrnd(muT,SigmaT,n0).t();

   set_seed(arm_seedR);
   mat yR = arma::mvnrnd(muR,SigmaR,n1).t();

   mat mu0 = mean(yT,0);
   mat mu1 = mean(yR,0);
   mat sd0 = stddev(yT,0,0);
   mat sd1 = stddev(yR,0,0);
   mat sde;
   mat df;
   mat tlb;
   mat tub;
   mat plb;
   mat pub;

   if(vareq == true){
     // Schuirmann’s test
     sde = arma::pow(((n0 - 1)*arma::pow(sd0, 2) + (n1 - 1)*arma::pow(sd1, 2))/(n0 + n1 - 2.0)*(1.0/n0 + 1.0/n1),0.5);
     df = datum::nan;
     tlb = (mu0 - mu1 - lequi_tol)/sde ;
     tub = (mu0 - mu1 - uequi_tol)/sde ;
     plb = ptv(tlb, 1.0*(n0 + n1 - 2),false);
     pub = ptv(tub, 1.0*(n0 + n1 - 2),true);
   }else{
     sde = sqrt(arma::pow(sd0, 2)/n0 + arma::pow(sd1, 2)/n1);
     df = (arma::pow(sd0, 2)/n0 + arma::pow(sd1, 2)/n1)%(arma::pow(sd0, 2)/n0 + arma::pow(sd1, 2)/n1)/(arma::pow(arma::pow(sd0, 2)/n0,2)/(n0 - 1.0) + arma::pow(arma::pow(sd1, 2)/n1,2)/(n1 - 1.0) );

     tlb = (mu0 - mu1 - lequi_tol)/sde;
     tub = (mu0 - mu1 - uequi_tol)/sde;
     plb = ptvdf(tlb, df, false);
     pub = ptvdf(tub, df, true);

   }


   // Calculate p-value

   mat ptost = max(plb, pub);

   ptost.replace(datum::nan, 0); // in case of NA values due to no sd in some studies
   mat alpha0 = conv_to<mat>::from(alpha);
   mat tbioq  = conv_to<mat>::from((ptost < alpha0));

   // Call the check_equivalence function to determine if equivalence is established
   arma::mat totaly = check_equivalence(typey, adseq, tbioq, k);

   // Combine results into a response matrix
   arma::mat response0 = join_rows(totaly, tbioq);
   arma::mat response1 = join_rows(mu0, mu1);
   arma::mat response2 = join_rows(sd0, sd1);
   arma::mat response3 = join_rows(response0, response1);

   return join_rows(response3, response2);
}

//' @title Simulate a Parallel Design and Test Ratio of Means (ROM)
//'
//' @description
//' Simulates a parallel-group design and performs equivalence testing using the ratio of means (ROM) approach.
//' This function evaluates whether the treatment and reference groups are equivalent based on predefined
//' equivalence margins and hypothesis testing criteria.
//'
//' @param n integer number of subjects per arm
//' @param muT vector mean of endpoints on treatment arm
//' @param muR vector mean of endpoints on reference arm
//' @param SigmaT matrix covar-variance matrix on treatment arm across endpoints
//' @param SigmaR matrix covar-variance matrix on reference arm across endpoints
//' @param lequi_tol vector  lower equivalence tolerance band across endpoints
//' @param uequi_tol vector  upper equivalence tolerance band across endpoints
//' @param alpha vector alpha value across endpoints
//' @param dropout vector of size 2 with dropout proportion per arm (T,R)
//' @param typey vector with positions of primary endpoints
//' @param adseq boolean is used a sequential adjustment?
//' @param k integer minimum number of equivalent endpoints
//' @param arm_seedT integer seed for the simulation on treatment arm
//' @param arm_seedR integer seed for the simulation on reference arm
//' @param TART double treatment allocation rate for the treatment arm
//' @param TARR double treatment allocation rate for the reference arm
//' @param vareq Boolean. If \code{TRUE}, assumes equal variance between arms and applies Schuirmann's two one-sided tests (TOST) for equivalence using a pooled variance.
//'
//' @return A numeric matrix containing the simulated hypothesis test results.
//' The first column represents the overall equivalence decision, where 1 indicates
//' success and 0 indicates failure. The subsequent columns contain the hypothesis
//' test results for each endpoint, followed by mean estimates for the reference and
//' treatment groups, and standard deviations for the reference and treatment groups.
//'
//' @details
//' The function simulates a parallel-group study design and evaluates equivalence
//' using the ratio of means (ROM) approach. It accounts for dropout rates and
//' treatment allocation proportions while generating simulated data based on the
//' specified covariance structure. The test statistics are computed, and a final
//' equivalence decision is made based on the predefined number of required significant
//' endpoints (\code{k}). If sequential testing (\code{adseq}) is enabled, primary endpoints
//' must establish equivalence before secondary endpoints are evaluated.
//' When \code{vareq = TRUE}, the test assumes equal variances between groups and
//' applies Schuirmann's two one-sided tests (TOST).
//'
//' @export
// [[Rcpp::export]]
arma::mat test_par_rom(const int n,
                       const arma::vec& muT,
                       const arma::vec& muR,
                       const arma::mat& SigmaT,
                       const arma::mat& SigmaR,
                       const arma::rowvec& lequi_tol,
                       const arma::rowvec& uequi_tol,
                       const arma::rowvec& alpha,
                       const arma::vec& dropout,
                       const arma::uvec& typey,
                       const bool adseq,
                       const int k,
                       const int arm_seedT,
                       const int arm_seedR,
                       const double TART,
                       const double TARR,
                       const bool vareq){
  // Transform drop out
   int n0i = ceil(n*TART);
   int n1i = ceil(n*TARR);

   int n0 = ceil((1 - dropout[0])* n0i);
   if (n0 < 2) n0 = 2;

   int n1 = ceil((1 - dropout[1])*n1i);
   if (n1 < 2) n1 = 2;

   RNGScope scope;
   Environment base_env("package:base");
   Function set_seed = base_env["set.seed"];

   // Generate data based on the provided covariance matrix and means
   set_seed(arm_seedT);
   mat yT = arma::mvnrnd(muT,SigmaT,n0).t();

   set_seed(arm_seedR);
   mat yR = arma::mvnrnd(muR,SigmaR,n1).t();

   mat mu0 = mean(yT,0);
   mat mu1 = mean(yR,0);
   mat sd0 = stddev(yT,0,0);
   mat sd1 = stddev(yR,0,0);
   mat sde_l;
   mat sde_u;
   int df = n0 + n1 - 2;

   if(vareq == true){
     // Schuirmann’s test
     sde_l = arma::pow(((n0 - 1)*arma::pow(sd0, 2) + (n1 - 1)*arma::pow(sd1, 2))/(n0 + n1 - 2.0)%(1.0/n0 + arma::pow(lequi_tol, 2)/n1),0.5);
     sde_u = arma::pow(((n0 - 1)*arma::pow(sd0, 2) + (n1 - 1)*arma::pow(sd1, 2))/(n0 + n1 - 2.0)%(1.0/n0 + arma::pow(uequi_tol, 2)/n1),0.5);

   }else{
     sde_l = arma::pow(arma::pow(sd0, 2)/n0 + arma::pow(lequi_tol, 2)%arma::pow(sd1, 2)/n1,0.5);
     sde_u = arma::pow(arma::pow(sd0, 2)/n0 + arma::pow(uequi_tol, 2)%arma::pow(sd1, 2)/n1,0.5);
   }

   mat tlb = (mu0 - mu1%lequi_tol)/sde_l ;
   mat tub = (mu0 - mu1%uequi_tol)/sde_u ;

   // Calculate p-value
   mat plb = ptv(tlb,df,false);
   mat pub = ptv(tub,df,true);
   mat ptost = max(plb, pub);

   ptost.replace(datum::nan, 0); // in case of NA values due to no sd in some studies
   mat alpha0 = conv_to<mat>::from(alpha);
   mat tbioq  = conv_to<mat>::from((ptost < alpha0));

   // Call the check_equivalence function to determine if equivalence is established
   arma::mat totaly = check_equivalence(typey, adseq, tbioq, k);

   // Combine results into a response matrix
   arma::mat response0 = join_rows(totaly, tbioq);
   arma::mat response1 = join_rows(mu0, mu1);
   arma::mat response2 = join_rows(sd0, sd1);
   arma::mat response3 = join_rows(response0, response1);

   return join_rows(response3, response2);
}


//' @title Run Simulations for a Parallel Design with Difference of Means (DOM) test
//'
//' @description
//' This function simulates a parallel-group trial across multiple iterations.
//' It evaluates equivalence across multiple endpoints using the
//' Difference of Means (DOM) test.
//'
//' @param nsim Integer. The number of simulations to run.
//' @param n Integer. The sample size per arm (before dropout).
//' @param muT arma::vec. Mean vector for the treatment arm.
//' @param muR arma::vec. Mean vector for the reference arm.
//' @param SigmaT arma::mat. Covariance matrix for the treatment arm.
//' @param SigmaR arma::mat. Covariance matrix for the reference arm.
//' @param lequi_tol arma::rowvec. Lower equivalence thresholds for each endpoint.
//' @param uequi_tol arma::rowvec. Upper equivalence thresholds for each endpoint.
//' @param alpha arma::rowvec. Significance level for each endpoint.
//' @param dropout arma::vec. Dropout rates for each arm (T, R).
//' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
//' @param adseq Boolean. If \code{TRUE}, applies sequential (hierarchical) testing.
//' @param k Integer. Minimum number of endpoints required for equivalence.
//' @param arm_seed_T arma::ivec. Random seed vector for the treatment group (one per simulation).
//' @param arm_seed_R arma::ivec. Random seed vector for the reference group (one per simulation).
//' @param TART Double. Treatment allocation ratio (proportion of subjects in treatment arm).
//' @param TARR Double. Reference allocation ratio (proportion of subjects in reference arm).
//' @param vareq Boolean. If \code{TRUE}, assumes equal variances across treatment and reference groups.
//'
//' @details
//' Equivalence testing uses either the Difference of Means (DOM) test,
//' applying predefined equivalence thresholds and significance levels. When hierarchical testing (\code{adseq})
//' is enabled, all primary endpoints must demonstrate equivalence before secondary endpoints are evaluated.
//' Dropout rates are incorporated into the sample size calculation to ensure proper adjustment.
//' Randomization is controlled through separate random seeds for the treatment and reference groups,
//' enhancing reproducibility.
//'
//' @return
//' The function returns an arma::mat storing simulation results row-wise for consistency
//' with R's output format. The first row (\code{totaly}) contains the overall equivalence decision
//' (1 for success, 0 for failure). The subsequent rows include equivalence deicisons for each endpoint,
//' mean estimates for both treatment and reference groups, and corresponding standard deviations.
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
arma::mat run_simulations_par_dom(const int nsim,
                                  const int n,
                               const arma::vec& muT,
                               const arma::vec& muR,
                               const arma::mat& SigmaT,
                               const arma::mat& SigmaR,
                               const arma::rowvec& lequi_tol,
                               const arma::rowvec& uequi_tol,
                               const arma::rowvec& alpha,
                               const arma::vec& dropout,
                               const arma::uvec& typey,
                               const bool adseq,
                               const int k,
                               const arma::ivec& arm_seed_T,
                               const arma::ivec& arm_seed_R,
                               const double TART,
                               const double TARR,
                               const bool vareq) {

   // **Determine number of endpoints**
   int num_endpoints = muT.n_elem; // Assuming muT and muR have the same number of elements

   // **Define the number of columns in result matrix**
   int num_cols = 1 + num_endpoints * 5; // totaly + 5 columns per endpoint

   // **Initialize result matrix**
   arma::mat results(nsim, num_cols, arma::fill::zeros);

   for (int i = 0; i < nsim; i++) {
     arma::mat outtest = test_par_dom(n, muT, muR,
                                      SigmaT, SigmaR,
                                      lequi_tol, uequi_tol,
                                      alpha, dropout,
                                      typey, adseq, k,
                                      arm_seed_T(i), arm_seed_R(i),
                                      TART, TARR, vareq);


     // **Store results in the output matrix**
     results.row(i) = outtest;
   }

   // Transpose results to match R's output format
   return results.t(); // Transpose before returning
}

//' @title Run Simulations for a Parallel Design with Ratio of Means (ROM) test
//'
//' @description
//' This function simulates a parallel-group trial across multiple iterations.
//' It evaluates equivalence across multiple endpoints using the
//' Ratio of Means (ROM) test.
//'
//' @param nsim Integer. The number of simulations to run.
//' @param n Integer. The sample size per arm (before dropout).
//' @param muT arma::vec. Mean vector for the treatment arm.
//' @param muR arma::vec. Mean vector for the reference arm.
//' @param SigmaT arma::mat. Covariance matrix for the treatment arm.
//' @param SigmaR arma::mat. Covariance matrix for the reference arm.
//' @param lequi_tol arma::rowvec. Lower equivalence thresholds for each endpoint.
//' @param uequi_tol arma::rowvec. Upper equivalence thresholds for each endpoint.
//' @param alpha arma::rowvec. Significance level for each endpoint.
//' @param dropout arma::vec. Dropout rates for each arm (T, R).
//' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
//' @param adseq Boolean. If \code{TRUE}, applies sequential (hierarchical) testing.
//' @param k Integer. Minimum number of endpoints required for equivalence.
//' @param arm_seed_T arma::ivec. Random seed vector for the treatment group (one per simulation).
//' @param arm_seed_R arma::ivec. Random seed vector for the reference group (one per simulation).
//' @param TART Double. Treatment allocation ratio (proportion of subjects in treatment arm).
//' @param TARR Double. Reference allocation ratio (proportion of subjects in reference arm).
//' @param vareq Boolean. If \code{TRUE}, assumes equal variances across treatment and reference groups.
//'
//' @details
//' Equivalence testing uses either the Ratio of Means (ROM) test,
//' applying predefined equivalence thresholds and significance levels. When hierarchical testing (\code{adseq})
//' is enabled, all primary endpoints must demonstrate equivalence before secondary endpoints are evaluated.
//' Dropout rates are incorporated into the sample size calculation to ensure proper adjustment.
//' Randomization is controlled through separate random seeds for the treatment and reference groups,
//' enhancing reproducibility.
//'
//' @return
//' The function returns an arma::mat storing simulation results row-wise for consistency
//' with R's output format. The first row (\code{totaly}) contains the overall equivalence decision
//' (1 for success, 0 for failure). The subsequent rows include equivalence decisions for each endpoint,
//' mean estimates for both treatment and reference groups, and corresponding standard deviations.
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
arma::mat run_simulations_par_rom(const int nsim,
                                   const int n,
                                   const arma::vec& muT,
                                   const arma::vec& muR,
                                   const arma::mat& SigmaT,
                                   const arma::mat& SigmaR,
                                   const arma::rowvec& lequi_tol,
                                   const arma::rowvec& uequi_tol,
                                   const arma::rowvec& alpha,
                                   const arma::vec& dropout,
                                   const arma::uvec& typey,
                                   const bool adseq,
                                   const int k,
                                   const arma::ivec& arm_seed_T,
                                   const arma::ivec& arm_seed_R,
                                   const double TART,
                                   const double TARR,
                                   const bool vareq) {

   // **Determine number of endpoints**
   int num_endpoints = muT.n_elem; // Assuming muT and muR have the same number of elements

   // **Define the number of columns in result matrix**
   int num_cols = 1 + num_endpoints * 5; // totaly + 5 columns per endpoint

   // **Initialize result matrix**
   arma::mat results(nsim, num_cols, arma::fill::zeros);

   for (int i = 0; i < nsim; i++) {
     arma::mat outtest = test_par_rom(n, muT, muR,
                                      SigmaT, SigmaR,
                                      lequi_tol, uequi_tol,
                                      alpha, dropout,
                                      typey, adseq, k,
                                      arm_seed_T(i), arm_seed_R(i),
                                      TART, TARR, vareq);


     // **Store results in the output matrix**
     results.row(i) = outtest;
   }

   // Transpose results to match R's output format
   return results.t(); // Transpose before returning
 }



//' @title Run Simulations for a 2x2 Crossover Design with Difference of Means (DOM) test
//'
//' @description
//' This function simulates a 2x2 crossover trial across multiple iterations.
//' It evaluates equivalence across multiple endpoints using the
//' Difference of Means (DOM) test.
//'
//' @param nsim Integer. The number of simulations to run.
//' @param n Integer. The sample size per period.
//' @param muT Numeric vector. Mean outcomes for the active treatment.
//' @param muR Numeric vector. Mean outcomes for the reference treatment.
//' @param SigmaW Numeric matrix. Within-subject covariance matrix for endpoints.
//' @param lequi_tol Numeric vector. Lower equivalence thresholds for each endpoint.
//' @param uequi_tol Numeric vector. Upper equivalence thresholds for each endpoint.
//' @param alpha Numeric vector. Significance levels for hypothesis testing across endpoints.
//' @param sigmaB Numeric. Between-subject variance for the crossover model.
//' @param dropout Numeric vector of size 2. Dropout rates for each sequence.
//' @param Eper Numeric vector. Expected period effects for each sequence.
//' @param Eco Numeric vector. Expected carryover effects for each sequence.
//' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
//' @param adseq Logical. If \code{TRUE}, applies sequential (hierarchical) testing.
//' @param k Integer. Minimum number of endpoints required for equivalence.
//' @param arm_seed Integer vector. Random seed for each simulation.
//'
//' @details
//' This function evaluates equivalence using the Difference of Means (DOM) test.
//' Equivalence is determined based on predefined lower (\code{lequi_tol}) and upper (\code{uequi_tol}) equivalence thresholds,
//' and hypothesis testing is conducted at the specified significance level (\code{alpha}).
//' If \code{adseq} is \code{TRUE}, primary endpoints must establish equivalence before secondary endpoints are evaluated.
//' The sample size per period is adjusted based on dropout rates, ensuring valid study conclusions.
//' The simulation incorporates within-subject correlation using \code{SigmaW} and accounts for between-subject variance with \code{sigmaB}.
//' Expected period effects (\code{Eper}) and carryover effects (\code{Eco}) are included in the model.
//' A fixed random seed (\code{arm_seed}) is used to ensure reproducibility across simulations.
//'
//' @return
//' A numeric matrix where each column stores simulation results:
//' The first row (\code{totaly}) represents the overall equivalence decision (1 = success, 0 = failure).
//' Subsequent rows contain equivalence decisions per endpoint,
//' mean estimates for the treatment group, mean estimates for the reference group,
//' standard deviations for treatment, and standard deviations for reference.
//'
//' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
arma::mat run_simulations_2x2_dom(const int nsim,
                                  const int n,
                               const arma::vec& muT,
                               const arma::vec& muR,
                               const arma::mat& SigmaW,
                               const arma::rowvec& lequi_tol,
                               const arma::rowvec& uequi_tol,
                               const arma::rowvec& alpha,
                               const double sigmaB,
                               const arma::vec& dropout,
                               const arma::vec& Eper,
                               const arma::vec& Eco,
                               const arma::uvec& typey,
                               const bool adseq,
                               const int k,
                               const arma::ivec& arm_seed){

   // **Determine number of endpoints**
   int num_endpoints = muT.n_elem; // Assuming muT and muR have the same number of elements

   // **Define the number of columns in result matrix**
   int num_cols = 1 + num_endpoints * 5; // totaly + 5 columns per endpoint

   // **Initialize result matrix**
   arma::mat results(nsim, num_cols, arma::fill::zeros);

   for (int i = 0; i < nsim; i++) {
     arma::mat outtest = test_2x2_dom(n, muT, muR, SigmaW, lequi_tol, uequi_tol,
                                      alpha, sigmaB, dropout, Eper, Eco, typey,
                                      adseq, k, arm_seed(i));

     // **Store results in the output matrix**
     results.row(i) = outtest;
   }

   // Transpose results to match R's output format
   return results.t(); // Transpose before returning
 }

//' @title Run Simulations for a 2x2 Crossover Design with Ratio of Means (ROM) test
//'
//' @description
//' This function simulates a 2x2 crossover trial across multiple iterations.
//' It evaluates equivalence across multiple endpoints using the
//' Ratio of Means (ROM) test.
//'
//' @param nsim Integer. The number of simulations to run.
//' @param n Integer. The sample size per period.
//' @param muT Numeric vector. Mean outcomes for the active treatment.
//' @param muR Numeric vector. Mean outcomes for the reference treatment.
//' @param SigmaW Numeric matrix. Within-subject covariance matrix for endpoints.
//' @param lequi_tol Numeric vector. Lower equivalence thresholds for each endpoint.
//' @param uequi_tol Numeric vector. Upper equivalence thresholds for each endpoint.
//' @param alpha Numeric vector. Significance levels for hypothesis testing across endpoints.
//' @param sigmaB Numeric. Between-subject variance for the crossover model.
//' @param dropout Numeric vector of size 2. Dropout rates for each sequence.
//' @param Eper Numeric vector. Expected period effects for each sequence.
//' @param Eco Numeric vector. Expected carryover effects for each sequence.
//' @param typey Integer vector indicating the classification of each endpoint, where \code{1} corresponds to a primary endpoint and \code{2} corresponds to a secondary endpoint.
//' @param adseq Logical. If \code{TRUE}, applies sequential (hierarchical) testing.
//' @param k Integer. Minimum number of endpoints required for equivalence.
//' @param arm_seed Integer vector. Random seed for each simulation.
//'
//' @details
//' This function evaluates equivalence using the Ratio of Means (ROM) test.
//' Equivalence is determined based on predefined lower \code{lequi_tol} and upper \code{uequi_tol} equivalence thresholds,
//' and hypothesis testing is conducted at the specified significance level \code{alpha}.
//' If \code{adseq} is \code{TRUE}, primary endpoints must establish equivalence before secondary endpoints are evaluated.
//' The sample size per period is adjusted based on dropout rates, ensuring valid study conclusions.
//' The simulation incorporates within-subject correlation using \code{SigmaW} and accounts for between-subject variance with \code{sigmaB}.
//' Expected period effects \code{Eper} and carryover effects \code{Eco} are included in the model.
//' A fixed random seed \code{arm_seed} is used to ensure reproducibility across simulations.//'
//'
//' @return
//' A numeric matrix where each column stores simulation results:
//' The first row (\code{totaly}) represents the overall equivalence decision (1 = success, 0 = failure).
//' Subsequent rows contain equivalence decisions per endpoint,
//' mean estimates for the treatment group, mean estimates for the reference group,
//' standard deviations for treatment, and standard deviations for reference.
//'
//'  @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
//' @export
// [[Rcpp::export]]
 arma::mat run_simulations_2x2_rom(const int nsim,
                                   const int n,
                                   const arma::vec& muT,
                                   const arma::vec& muR,
                                   const arma::mat& SigmaW,
                                   const arma::rowvec& lequi_tol,
                                   const arma::rowvec& uequi_tol,
                                   const arma::rowvec& alpha,
                                   const double sigmaB,
                                   const arma::vec& dropout,
                                   const arma::vec& Eper,
                                   const arma::vec& Eco,
                                   const arma::uvec& typey,
                                   const bool adseq,
                                   const int k,
                                   const arma::ivec& arm_seed){

   // **Determine number of endpoints**
   int num_endpoints = muT.n_elem; // Assuming muT and muR have the same number of elements

   // **Define the number of columns in result matrix**
   int num_cols = 1 + num_endpoints * 5; // totaly + 5 columns per endpoint

   // **Initialize result matrix**
   arma::mat results(nsim, num_cols, arma::fill::zeros);

   for (int i = 0; i < nsim; i++) {
     arma::mat outtest = test_2x2_rom(n, muT, muR, SigmaW, lequi_tol, uequi_tol,
                                      alpha, sigmaB, dropout, Eper, Eco, typey,
                                      adseq, k, arm_seed(i));

     // **Store results in the output matrix**
     results.row(i) = outtest;
   }

   // Transpose results to match R's output format
   return results.t(); // Transpose before returning
 }

RCPP_MODULE(test)
{
   function("rcpp_test_2x2_dom", &test_2x2_dom);
   function("rcpp_test_2x2_rom", &test_2x2_rom);
   function("rcpp_test_par_dom", &test_par_dom);
   function("rcpp_test_par_rom", &test_par_rom);
}
