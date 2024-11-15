#include <iostream>
#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

using namespace Rcpp;
using namespace arma;


//' ptv
//'
//' Description : to get the p-value with fixed degrees of freedom across the endpoints
//'
//' @param x mat(vector) random variable
//' @param df double degrees of freedom
//' @param lower boolean, calculate the probability of the region where the random variable is less than or equal to x
//' @return mat(vector) with the cumulative distribution function (p-values)
//' @export
// [[Rcpp::export]]
arma::mat ptv(arma::mat x, double df, bool lower) {
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
//' @description This function computes the cumulative distribution function (p-values) for a given random variable `x`
//' and corresponding degrees of freedom `df` using the t-distribution. The function can compute the lower or upper
//' tail probabilities depending on the value of the `lower` parameter.
//'
//' @param x arma::mat (vector) - A matrix or vector of random variable values for which the p-values will be calculated.
//' @param df arma::mat (vector) - A matrix or vector of degrees of freedom for the t-distribution, matching the size of `x`.
//' @param lower bool - If `TRUE`, calculates the lower-tail probability (P(T <= x)); if `FALSE`, calculates the upper-tail probability.
//' @return arma::mat (vector) - A matrix containing the computed cumulative distribution function (p-values) for each element in `x`.
//' The result is returned as a 1xN matrix, where N is the number of elements in `x`.
//' @export
// [[Rcpp::export]]
arma::mat ptvdf(arma::mat x, arma::mat df, bool lower) {
  std::size_t n = x.n_elem;  // Get the number of elements in x
  arma::vec y(n);  // Create an Armadillo vector for the result

  for (std::size_t i = 0; i < n; ++i) {
    // Directly pass individual elements of x and df to Rcpp::pt()
    y[i] = R::pt(x(i), df(i), lower, false);  // Use R's pt function directly
  }

  // Reshape y into a 1xN matrix and return
  return arma::reshape(y, 1, n);
}

//' test_2x2_dom_cpp
//'
//' Description : to simulate a 2x2 design and return the p-tost for difference of means (DOM)
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
//' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
//' @export
// [[Rcpp::export]]
arma::mat test_2x2_dom(int n, arma::vec muT, arma::vec muR,
                           arma::mat SigmaW, arma::rowvec lequi_tol, arma::rowvec uequi_tol,
                           arma::rowvec alpha, double sigmaB, arma::vec dropout,
                           arma::vec Eper, arma::vec Eco, arma::uvec typey, bool adseq, int k, int arm_seed){

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


  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(arm_seed);

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

    // primary endpoints in case of sequencial adjustment
    bool sumpe;
    double sumtypey = accu(tbioq.cols(typey));
    double lentypey = typey.n_cols;
    if(adseq == TRUE){
      sumpe =  sumtypey == lentypey;
    }else{
      sumpe = TRUE;
    }

    bool sumt = accu(tbioq) >= k;

    mat totaly(1,1);

    if(sumt&sumpe){
      totaly(0, 0) = 1;
    }else{
      totaly(0, 0) = 0;
    }


    mat response0 = join_rows<mat>(totaly,tbioq);
    mat response1 = join_rows<mat>(mut0,mut1);
    mat response2 = join_rows<mat>(sdw,sdb);
    mat response3 = join_rows<mat>(response0,response1);

    return join_rows<mat>(response3,response2);

}


//' test_2x2_rom_cpp
//'
//' Description : to simulate a 2x2 design and return the p-tost for ratio of means (ROM)
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
//' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
//' @export
// [[Rcpp::export]]
arma::mat test_2x2_rom(int n, arma::vec muT, arma::vec muR,
                             arma::mat SigmaW, arma::rowvec lequi_tol, arma::rowvec uequi_tol,
                             arma::rowvec alpha, double sigmaB, arma::vec dropout,
                             arma::vec Eper, arma::vec Eco, arma::uvec typey, bool adseq, int k, int arm_seed){

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

  RNGScope scope;
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(arm_seed);

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

  // primary endpoints in case of sequencial adjustment
  bool sumpe;
  double sumtypey = accu(tbioq.cols(typey));
  double lentypey = typey.n_cols;
  if(adseq == TRUE){
    sumpe =  sumtypey == lentypey;
  }else{
    sumpe = TRUE;
  }


  bool sumt = accu(tbioq) >= k;

  mat totaly(1,1);

  if(sumt&sumpe){
    totaly(0, 0) = 1;
  }else{
    totaly(0, 0) = 0;
  }



  mat response0 = join_rows<mat>(totaly,tbioq);
  mat response1 = join_rows<mat>(mut0,mut1);
  mat response2 = join_rows<mat>(vart0,vart1);
  mat response3 = join_rows<mat>(response2,covt0t1);
  mat response4 = join_rows<mat>(response0,response1);

  return join_rows<mat>(response4,response3);

}

//' test_par_dom_cpp
//'
//' Description : to simulate a parallel design and return the p-tost for differnece of means (DOM)
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
//' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
//' @export
// [[Rcpp::export]]
arma::mat test_par_dom(int n, arma::vec muT, arma::vec muR,
                           arma::mat SigmaT, arma::mat SigmaR,
                           arma::rowvec lequi_tol, arma::rowvec uequi_tol,
                           arma::rowvec alpha, arma::vec dropout,
                           arma::uvec typey, bool adseq, int k,
                           int arm_seedT, int arm_seedR,
                           double TART, double TARR,
                           bool vareq){


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

 // primary endpoints in case of sequencial adjustment
 bool sumpe;
 double sumtypey = accu(tbioq.cols(typey));
 double lentypey = typey.n_cols;
 if(adseq == TRUE){
   sumpe =  sumtypey == lentypey;
 }else{
   sumpe = TRUE;
 }


 bool sumt = accu(tbioq) >= k;

 mat totaly(1,1);

 if(sumt&sumpe){
   totaly(0, 0) = 1;
 }else{
   totaly(0, 0) = 0;
 }



  mat response0 = join_rows<mat>(totaly,tbioq);
  mat response1 = join_rows<mat>(mu0,mu1);
  mat response2 = join_rows<mat>(sd0,sd1);
  mat response3 = join_rows<mat>(response0,response1);

  return join_rows<mat>(response3,response2);

}

//' test_par_rom_cpp
//'
//' Description : to simulate a parallel design and return the p-tost for ratio of means (ROM)
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
//' @param vareq boolean assumed equivalence variance between arms for the t-test
//' @return mat(vector) with ptost and other simulated statistics such as mean (mu) and standard deviation(std) per sequence (0,1)-endpoint
//' @export
// [[Rcpp::export]]

arma::mat test_par_rom(int n, arma::vec muT, arma::vec muR,
                           arma::mat SigmaT, arma::mat SigmaR,
                           arma::rowvec lequi_tol, arma::rowvec uequi_tol,
                           arma::rowvec alpha,  arma::vec dropout,
                           arma::uvec typey, bool adseq,int k,
                           int arm_seedT, int arm_seedR,
                           double TART, double TARR,
                           bool vareq){


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

  // primary endpoints in case of sequencial adjustment
  bool sumpe;
  double sumtypey = accu(tbioq.cols(typey));
  double lentypey = typey.n_cols;
  if(adseq == TRUE){
    sumpe =  sumtypey == lentypey;
  }else{
    sumpe = TRUE;
  }

  bool sumt = accu(tbioq) >= k;

  mat totaly(1,1);

  if(sumt&sumpe){
    totaly(0, 0) = 1;
  }else{
    totaly(0, 0) = 0;
  }


  mat response0 = join_rows<mat>(totaly,tbioq);
  mat response1 = join_rows<mat>(mu0,mu1);
  mat response2 = join_rows<mat>(sd0,sd1);
  mat response3 = join_rows<mat>(response0,response1);

  return join_rows<mat>(response3,response2);

}


RCPP_MODULE(test)
 {
   function("rcpp_test_2x2_dom", &test_2x2_dom);
   function("rcpp_test_2x2_rom", &test_2x2_rom);
   function("rcpp_test_par_dom", &test_par_dom);
   function("rcpp_test_par_rom", &test_par_rom);
 }

