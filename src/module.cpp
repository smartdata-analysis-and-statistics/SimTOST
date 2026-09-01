#include <iostream>
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

using namespace Rcpp;
using namespace arma;

// Draw rows from N(mu, Sigma) using an explicit Cholesky factor. This avoids
// relying on arma::mvnrnd's internal dimension handling in repeated calls from
// the simulation loops.
arma::mat draw_mvnorm_rows(const arma::vec& mu, const arma::mat& Sigma,
                           const int n) {
  if (mu.n_elem == 0 || Sigma.n_rows != mu.n_elem ||
      Sigma.n_cols != mu.n_elem) {
    Rcpp::stop("Mean vector and covariance matrix dimensions do not agree.");
  }
  if (!Sigma.is_finite() || !Sigma.is_symmetric(1e-10)) {
    Rcpp::stop("Covariance matrix must be finite and symmetric.");
  }
  arma::mat chol_sigma;
  bool ok = arma::chol(chol_sigma, Sigma);
  if (!ok) {
    Rcpp::stop("Covariance matrix must be positive definite.");
  }
  arma::mat out = arma::randn(n, mu.n_elem) * chol_sigma;
  out.each_row() += mu.t();
  return out;
}

double draw_count_value(const int model, const double dispersion,
                        const double mu) {
  if (model == 0) return R::rpois(mu);
  if (model == 1) return ::Rf_rnbinom_mu(1.0 / dispersion, mu);
  Rcpp::stop("'model' must be 0 (Poisson) or 1 (negative-binomial).");
  return NA_REAL;
}

// Paired log-link analysis for a balanced 2x2 crossover.  The subject-level
// treatment contrast is formed from the two observed period counts.  Averaging
// the contrasts from the two sequences removes the period effect; the
// carry-over correction uses Eco = (reference carry-over, treatment
// carry-over).  Because the contrast is within subject, the log-normal
// subject effect generated with sigmaB cancels from the estimand and its
// residual sampling variation is retained through the empirical contrast
// variance.
bool crossover_count_equivalent(const int n_per_arm,
                                const double rate_test,
                                const double rate_reference,
                                const double exposure_test,
                                const double exposure_reference,
                                const double margin_lower,
                                const double margin_upper,
                                const int model,
                                const double dispersion_test,
                                const double dispersion_reference,
                                const double alpha,
                                const double sigmaB,
                                const arma::vec& Eper,
                                const arma::vec& Eco,
                                const arma::vec& dropout) {
  double sum[2] = {0.0, 0.0};
  double sumsq[2] = {0.0, 0.0};
  int n_complete[2] = {0, 0};

  for (int sequence = 0; sequence < 2; ++sequence) {
    const int n_subjects = R::rbinom(n_per_arm, 1.0 - dropout[sequence]);
    n_complete[sequence] = n_subjects;
    for (int subject = 0; subject < n_subjects; ++subject) {
      const double subject_effect = R::rnorm(0.0, sigmaB);
      double y_test = 0.0;
      double y_reference = 0.0;

      for (int period = 0; period < 2; ++period) {
        const bool receives_test = (sequence == 0 && period == 0) ||
                                    (sequence == 1 && period == 1);
        // Eco is ordered as reference carry-over, treatment carry-over.
        // In sequence T-R, the period-2 reference observation carries T;
        // in sequence R-T, the period-2 treatment observation carries R.
        const double carry = period == 0 ? 0.0 :
          (sequence == 0 ? Eco[1] : Eco[0]);
        const double log_rate =
          (receives_test ? std::log(rate_test) : std::log(rate_reference)) +
          Eper[period] + carry + subject_effect;
        const double exposure = receives_test ? exposure_test : exposure_reference;
        const double dispersion = receives_test ? dispersion_test : dispersion_reference;
        const double count = draw_count_value(
          model, dispersion, exposure * std::exp(log_rate));
        if (receives_test) y_test = count;
        else y_reference = count;
      }

      const double subject_contrast =
        std::log((y_test + 0.5) / exposure_test) -
        std::log((y_reference + 0.5) / exposure_reference);
      sum[sequence] += subject_contrast;
      sumsq[sequence] += subject_contrast * subject_contrast;
    }
  }

  if (n_complete[0] < 2 || n_complete[1] < 2) return false;
  const double mean0 = sum[0] / n_complete[0];
  const double mean1 = sum[1] / n_complete[1];
  const double var0 = std::max(0.0,
    (sumsq[0] - n_complete[0] * mean0 * mean0) /
      (n_complete[0] - 1.0));
  const double var1 = std::max(0.0,
    (sumsq[1] - n_complete[1] * mean1 * mean1) /
      (n_complete[1] - 1.0));

  // Period effects cancel in the average of the two sequence contrasts.
  // The following term removes the expected carry-over bias.
  const double estimate = 0.5 * (mean0 + mean1) +
    0.5 * (Eco[1] - Eco[0]);
  const double variance = 0.25 *
    (var0 / n_complete[0] + var1 / n_complete[1]);
  const double se = std::sqrt(std::max(variance, 1e-12));
  const double z = R::qnorm5(1.0 - alpha, 0.0, 1.0, 1, 0);

  return (estimate - std::log(margin_lower)) / se > z &&
    (std::log(margin_upper) - estimate) / se > z;
}

// Apply the same hierarchy used by the continuous simulation kernels. Without
// sequential testing, at least k endpoints must pass. With sequential testing,
// all primary endpoints form a gate and the secondary family must contain at
// least k passing endpoints. If no secondary endpoints exist, the primary
// family is treated as a co-primary intersection.
bool count_endpoint_decision(const arma::ivec& typey,
                             const bool adseq,
                             const std::vector<int>& endpoint_passed,
                             const int k) {
  int endpoint_successes = 0;
  for (const int passed : endpoint_passed) endpoint_successes += passed;
  if (!adseq || typey.n_elem == 0 || arma::any(typey < 0))
    return endpoint_successes >= k;

  int n_primary = 0;
  int n_secondary = 0;
  int primary_successes = 0;
  int secondary_successes = 0;
  for (unsigned int j = 0; j < endpoint_passed.size(); ++j) {
    if (typey[j] == 1) {
      ++n_primary;
      primary_successes += endpoint_passed[j];
    } else if (typey[j] == 2) {
      ++n_secondary;
      secondary_successes += endpoint_passed[j];
    }
  }
  if (n_primary == 0) return secondary_successes >= k;
  if (n_secondary == 0) return primary_successes == n_primary;
  return primary_successes == n_primary && secondary_successes >= k;
}

// [[Rcpp::export]]
Rcpp::NumericVector count_power_cpp(const int n_per_arm,
                                    const double rate_test,
                                    const double rate_reference,
                                    const double exposure_test,
                                    const double exposure_reference,
                                    const double margin_lower,
                                    const double margin_upper,
                                    const int model,
                                    const double dispersion_test,
                                    const double dispersion_reference,
                                    const double alpha,
                                    const int nsim,
                                    const int design,
                                    const double sigmaB,
                                    const arma::vec& Eper,
                                    const arma::vec& Eco,
                                    const arma::vec& dropout) {
  if (n_per_arm < 2 || rate_test <= 0 || rate_reference <= 0 ||
      exposure_test <= 0 || exposure_reference <= 0 ||
      margin_lower <= 0 || margin_lower >= margin_upper ||
      alpha <= 0 || alpha >= 1 || nsim < 1 ||
      (model == 1 && (dispersion_test <= 0 || dispersion_reference <= 0)) || sigmaB < 0 ||
      Eper.n_elem != 2 || Eco.n_elem != 2 || dropout.n_elem != 2 ||
      any(dropout < 0) || any(dropout >= 1) ||
      (design != 0 && design != 1)) {
    Rcpp::stop("Invalid count simulation parameters.");
  }

  RNGScope scope;
  const int subjects_per_treatment = design == 1 ? 2 * n_per_arm : n_per_arm;
  const double total_exposure_test = subjects_per_treatment * exposure_test;
  const double total_exposure_reference = subjects_per_treatment * exposure_reference;
  const double z = R::qnorm5(1.0 - alpha, 0.0, 1.0, 1, 0);
  const double lower = std::log(margin_lower);
  const double upper = std::log(margin_upper);
  int successes = 0;
  int endpoint_successes = 0;

  for (int i = 0; i < nsim; ++i) {
    double y_test;
    double y_reference;
    if (design == 0) {
      y_test = draw_count_value(model, dispersion_test,
                                total_exposure_test * rate_test);
      y_reference = draw_count_value(model, dispersion_reference,
                                     total_exposure_reference * rate_reference);
    } else {
      const bool passed = crossover_count_equivalent(
        n_per_arm, rate_test, rate_reference, exposure_test,
        exposure_reference, margin_lower, margin_upper, model,
        dispersion_test, dispersion_reference, alpha, sigmaB, Eper, Eco,
        dropout);
      if (passed) ++successes;
      if (passed) ++endpoint_successes;
      continue;
    }

    const double test_cc = y_test + 0.5;
    const double reference_cc = y_reference + 0.5;
    const double log_rr = std::log((test_cc / total_exposure_test) /
                                   (reference_cc / total_exposure_reference));
    const double se = std::sqrt(1.0 / test_cc + 1.0 / reference_cc);
    if ((log_rr - lower) / se > z &&
        (upper - log_rr) / se > z) {
      ++successes;
      ++endpoint_successes;
    }
  }

  Rcpp::NumericVector out(4);
  out[0] = static_cast<double>(successes) / nsim;
  out[1] = successes;
  out[2] = endpoint_successes;
  out[3] = successes;
  out.names() = Rcpp::CharacterVector::create(
    "power", "successes", "endpoint_success_1", "comparison_success_1");
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector count_power_multi_cpp(
    const int n_per_arm,
    const Rcpp::NumericVector& rate_test,
    const Rcpp::NumericVector& rate_reference,
    const Rcpp::NumericVector& exposure_test,
    const Rcpp::NumericVector& exposure_reference,
    const Rcpp::NumericVector& margin_lower,
    const Rcpp::NumericVector& margin_upper,
    const int model,
    const Rcpp::NumericVector& dispersion_test,
    const Rcpp::NumericVector& dispersion_reference,
    const Rcpp::NumericVector& alpha,
    const int nsim,
    const int design,
    const arma::mat& endpoint_corr,
    const arma::ivec& typey,
    const bool adseq,
    const int k,
    const double sigmaB,
    const arma::vec& Eper,
    const arma::vec& Eco,
    const arma::vec& dropout) {
  const int m = rate_test.size();
  if (m < 2 || rate_reference.size() != m || exposure_test.size() != m ||
      exposure_reference.size() != m || dispersion_test.size() != m ||
      dispersion_reference.size() != m ||
      margin_lower.size() != m || margin_upper.size() != m ||
      alpha.size() != m || endpoint_corr.n_rows != m ||
      endpoint_corr.n_cols != m || typey.n_elem != static_cast<unsigned int>(m) ||
      k < 1 || k > m) {
      Rcpp::stop("Multi-endpoint count inputs have incompatible dimensions.");
  }
  if (adseq && (arma::any(typey < 1) || arma::any(typey > 2)))
    Rcpp::stop("'typey' must contain only 1 (primary) or 2 (secondary) when sequential testing is enabled.");
  if ((design != 0 && design != 1) || sigmaB < 0 || Eper.n_elem != 2 ||
      Eco.n_elem != 2 || dropout.n_elem != 2 || any(dropout < 0) ||
      any(dropout >= 1)) {
    Rcpp::stop("'design' must be 0 (parallel) or 1 (2x2).");
  }
  for (int j = 0; j < m; ++j) {
    if (!R_finite(exposure_test[j]) || !R_finite(exposure_reference[j]) ||
        exposure_test[j] <= 0 || exposure_reference[j] <= 0 ||
        (model == 1 && (!R_finite(dispersion_test[j]) ||
                        !R_finite(dispersion_reference[j]) ||
                        dispersion_test[j] <= 0 ||
                        dispersion_reference[j] <= 0))) {
      Rcpp::stop("Exposure and dispersion must be valid for every endpoint.");
    }
  }
  if (!endpoint_corr.is_finite() || !endpoint_corr.is_symmetric(1e-10) ||
      arma::any(arma::abs(endpoint_corr.diag() - 1.0) > 1e-10))
    Rcpp::stop("'endpoint_corr' must be a finite correlation matrix.");
  arma::mat chol_corr;
  if (!arma::chol(chol_corr, endpoint_corr))
    Rcpp::stop("'endpoint_corr' must be positive definite.");
  RNGScope scope;
  const int subjects_per_treatment = design == 1 ? 2 * n_per_arm : n_per_arm;
  int successes = 0;
  std::vector<int> endpoint_successes(m, 0);
  int comparison_successes = 0;
  for (int i = 0; i < nsim; ++i) {
    std::vector<int> endpoint_passed(m, 0);
    if (design == 0) {
      std::vector<double> y_test(m), y_reference(m);
      arma::vec latent_test = chol_corr.t() * arma::randn<arma::vec>(m);
      arma::vec latent_reference = chol_corr.t() * arma::randn<arma::vec>(m);
      for (int j = 0; j < m; ++j) {
        const double p_test = std::min(1.0 - 1e-12,
                                       std::max(1e-12,
                                         R::pnorm5(latent_test[j], 0.0, 1.0, 1, 0)));
        const double p_reference = std::min(1.0 - 1e-12,
                                            std::max(1e-12,
                                              R::pnorm5(latent_reference[j], 0.0, 1.0, 1, 0)));
        const double total_exposure_test =
          subjects_per_treatment * exposure_test[j];
        const double total_exposure_reference =
          subjects_per_treatment * exposure_reference[j];
        y_test[j] = model == 0 ? R::qpois(p_test, total_exposure_test * rate_test[j], 1, 0) :
          R::qnbinom(p_test, 1.0 / dispersion_test[j],
                     (1.0 / dispersion_test[j]) /
                       ((1.0 / dispersion_test[j]) +
                        total_exposure_test * rate_test[j]), 1, 0);
        y_reference[j] = model == 0 ?
          R::qpois(p_reference, total_exposure_reference * rate_reference[j], 1, 0) :
          R::qnbinom(p_reference, 1.0 / dispersion_reference[j],
                     (1.0 / dispersion_reference[j]) /
                       ((1.0 / dispersion_reference[j]) +
                        total_exposure_reference * rate_reference[j]), 1, 0);
      }
      for (int j = 0; j < m; ++j) {
        const double test_cc = y_test[j] + 0.5;
        const double reference_cc = y_reference[j] + 0.5;
        const double log_rr = std::log((test_cc / (subjects_per_treatment * exposure_test[j])) /
                                       (reference_cc / (subjects_per_treatment * exposure_reference[j])));
        const double se = std::sqrt(1.0 / test_cc + 1.0 / reference_cc);
        const double z = R::qnorm5(1.0 - alpha[j], 0.0, 1.0, 1, 0);
        if ((log_rr - std::log(margin_lower[j])) / se > z &&
            (std::log(margin_upper[j]) - log_rr) / se > z) {
          endpoint_passed[j] = 1;
          ++endpoint_successes[j];
        }
      }
    } else {
      std::vector<double> sum0(m, 0.0), sum1(m, 0.0);
      std::vector<double> sumsq0(m, 0.0), sumsq1(m, 0.0);
      int n_complete[2] = {0, 0};
      for (int sequence = 0; sequence < 2; ++sequence) {
        n_complete[sequence] = R::rbinom(n_per_arm, 1.0 - dropout[sequence]);
        for (int subject = 0; subject < n_complete[sequence]; ++subject) {
          arma::vec subject_effect = sigmaB *
            (chol_corr.t() * arma::randn<arma::vec>(m));
          std::vector<double> y_test(m, 0.0), y_reference(m, 0.0);
          for (int period = 0; period < 2; ++period) {
            const bool receives_test = (sequence == 0 && period == 0) ||
                                        (sequence == 1 && period == 1);
            const double carry = period == 0 ? 0.0 :
              (sequence == 0 ? Eco[1] : Eco[0]);
            arma::vec latent = chol_corr.t() * arma::randn<arma::vec>(m);
            for (int j = 0; j < m; ++j) {
              const double p = std::min(1.0 - 1e-12,
                                        std::max(1e-12,
                                          R::pnorm5(latent[j], 0.0, 1.0, 1, 0)));
              const double rate = receives_test ? rate_test[j] : rate_reference[j];
              const double exposure = receives_test ? exposure_test[j] : exposure_reference[j];
              const double dispersion = receives_test ? dispersion_test[j] : dispersion_reference[j];
              const double mu = exposure * std::exp(std::log(rate) +
                Eper[period] + carry + subject_effect[j]);
              const double count = model == 0 ? R::qpois(p, mu, 1, 0) :
                R::qnbinom(p, 1.0 / dispersion,
                           (1.0 / dispersion) / ((1.0 / dispersion) + mu), 1, 0);
              if (receives_test) y_test[j] = count;
              else y_reference[j] = count;
            }
          }
          for (int j = 0; j < m; ++j) {
            const double contrast = std::log((y_test[j] + 0.5) / exposure_test[j]) -
              std::log((y_reference[j] + 0.5) / exposure_reference[j]);
            if (sequence == 0) {
              sum0[j] += contrast;
              sumsq0[j] += contrast * contrast;
            } else {
              sum1[j] += contrast;
              sumsq1[j] += contrast * contrast;
            }
          }
        }
      }
      for (int j = 0; j < m; ++j) {
        if (n_complete[0] < 2 || n_complete[1] < 2) continue;
        const double mean0 = sum0[j] / n_complete[0];
        const double mean1 = sum1[j] / n_complete[1];
        const double var0 = std::max(0.0,
          (sumsq0[j] - n_complete[0] * mean0 * mean0) /
            (n_complete[0] - 1.0));
        const double var1 = std::max(0.0,
          (sumsq1[j] - n_complete[1] * mean1 * mean1) /
            (n_complete[1] - 1.0));
        const double estimate = 0.5 * (mean0 + mean1) +
          0.5 * (Eco[1] - Eco[0]);
        const double se = std::sqrt(std::max(0.25 *
          (var0 / n_complete[0] + var1 / n_complete[1]), 1e-12));
        const double z = R::qnorm5(1.0 - alpha[j], 0.0, 1.0, 1, 0);
        if ((estimate - std::log(margin_lower[j])) / se > z &&
            (std::log(margin_upper[j]) - estimate) / se > z) {
          endpoint_passed[j] = 1;
          ++endpoint_successes[j];
        }
      }
    }
    if (count_endpoint_decision(typey, adseq, endpoint_passed, k)) {
      ++successes;
      ++comparison_successes;
    }
  }
  Rcpp::NumericVector out(3 + m);
  out[0] = static_cast<double>(successes) / nsim;
  out[1] = successes;
  for (int j = 0; j < m; ++j) out[2 + j] = endpoint_successes[j];
  out[2 + m] = comparison_successes;
  Rcpp::CharacterVector names(3 + m);
  names[0] = "power";
  names[1] = "successes";
  for (int j = 0; j < m; ++j)
    names[2 + j] = "endpoint_success_" + std::to_string(j + 1);
  names[2 + m] = "comparison_success_1";
  out.names() = names;
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector count_power_joint_cpp(
    const int n_per_arm,
    const Rcpp::NumericMatrix& rates,
    const Rcpp::NumericMatrix& exposure,
    const Rcpp::NumericMatrix& margin_lower,
    const Rcpp::NumericMatrix& margin_upper,
    const int model,
    const Rcpp::NumericMatrix& dispersion,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericMatrix& endpoint_corr,
    const Rcpp::IntegerMatrix& comparisons,
    const arma::ivec& typey,
    const bool adseq,
    const int k,
    const int nsim) {
  const int n_arms = rates.nrow();
  const int m = rates.ncol();
  const int n_comparisons = comparisons.nrow();
  if (n_per_arm < 2 || n_arms < 2 || m < 1 || n_comparisons < 1 ||
      rates.nrow() != exposure.nrow() || rates.ncol() != exposure.ncol() ||
      rates.nrow() != dispersion.nrow() || rates.ncol() != dispersion.ncol() ||
      margin_lower.nrow() != n_comparisons ||
      margin_upper.nrow() != n_comparisons ||
      margin_lower.ncol() != m || margin_upper.ncol() != m ||
      endpoint_corr.nrow() != m || endpoint_corr.ncol() != m ||
      comparisons.ncol() != 2 || typey.n_elem != static_cast<unsigned int>(m) ||
      k < 1 || k > m || nsim < 1 ||
      model < 0 || model > 1) {
    Rcpp::stop("Invalid joint count simulation dimensions or parameters.");
  }
  if (adseq && (arma::any(typey < 1) || arma::any(typey > 2)))
    Rcpp::stop("'typey' must contain only 1 (primary) or 2 (secondary) when sequential testing is enabled.");
  for (int i = 0; i < n_comparisons; ++i) {
    if (comparisons(i, 0) < 0 || comparisons(i, 0) >= n_arms ||
        comparisons(i, 1) < 0 || comparisons(i, 1) >= n_arms ||
        comparisons(i, 0) == comparisons(i, 1)) {
      Rcpp::stop("Comparison indices must refer to two distinct arms.");
    }
  }
  arma::mat corr = Rcpp::as<arma::mat>(endpoint_corr);
  if (!corr.is_finite() || !corr.is_symmetric(1e-10)) {
    Rcpp::stop("'endpoint_corr' must be finite and symmetric.");
  }
  arma::mat chol_corr;
  if (!arma::chol(chol_corr, corr)) {
    Rcpp::stop("'endpoint_corr' must be positive definite.");
  }
  for (int j = 0; j < m; ++j) {
    if (alpha[j] <= 0 || alpha[j] >= 1) {
      Rcpp::stop("Alpha must be valid for every endpoint.");
    }
    for (int comparison = 0; comparison < n_comparisons; ++comparison) {
      if (margin_lower(comparison, j) <= 0 ||
          margin_lower(comparison, j) >= margin_upper(comparison, j)) {
        Rcpp::stop("Margins must be valid for every comparison and endpoint.");
      }
    }
  }
  for (int i = 0; i < n_arms; ++i) {
    for (int j = 0; j < m; ++j) {
      if (rates(i, j) <= 0 || !std::isfinite(rates(i, j))) {
        Rcpp::stop("All arm-specific rates must be positive and finite.");
      }
      if (exposure(i, j) <= 0 || !std::isfinite(exposure(i, j)) ||
          (model == 1 && (dispersion(i, j) <= 0 ||
                          !std::isfinite(dispersion(i, j))))) {
        Rcpp::stop("All arm-specific exposure and dispersion values must be valid.");
      }
    }
  }

  RNGScope scope;
  int successes = 0;
  std::vector<int> endpoint_successes(n_comparisons * m, 0);
  std::vector<int> comparison_successes(n_comparisons, 0);
  arma::vec latent(m);
  arma::vec uniforms(m);
  std::vector< std::vector<double> > counts(
      n_arms, std::vector<double>(m, 0.0));

  for (int simulation = 0; simulation < nsim; ++simulation) {
    for (int arm = 0; arm < n_arms; ++arm) {
      latent = arma::randn<arma::vec>(m);
      latent = chol_corr.t() * latent;
      for (int endpoint = 0; endpoint < m; ++endpoint) {
        const double p = std::min(1.0 - 1e-12,
                                  std::max(1e-12,
                                           R::pnorm5(latent[endpoint], 0.0,
                                                     1.0, 1, 0)));
        const double mu = n_per_arm * exposure(arm, endpoint) * rates(arm, endpoint);
        if (model == 0) {
          counts[arm][endpoint] = R::qpois(p, mu, 1, 0);
        } else {
          const double size = 1.0 / dispersion(arm, endpoint);
          const double probability = size / (size + mu);
          counts[arm][endpoint] = R::qnbinom(p, size, probability, 1, 0);
        }
      }
    }

    bool all_comparisons_pass = true;
    for (int comparison = 0; comparison < n_comparisons; ++comparison) {
      const int test_arm = comparisons(comparison, 0);
      const int reference_arm = comparisons(comparison, 1);
      std::vector<int> endpoint_passed(m, 0);
      for (int endpoint = 0; endpoint < m; ++endpoint) {
        const double test_cc = counts[test_arm][endpoint] + 0.5;
        const double reference_cc = counts[reference_arm][endpoint] + 0.5;
        const double test_exposure = n_per_arm * exposure(test_arm, endpoint);
        const double reference_exposure = n_per_arm * exposure(reference_arm, endpoint);
        const double log_rr = std::log((test_cc / test_exposure) /
                                       (reference_cc / reference_exposure));
        const double se = std::sqrt(1.0 / test_cc + 1.0 / reference_cc);
        const double z = R::qnorm5(1.0 - alpha[endpoint], 0.0, 1.0, 1, 0);
        if ((log_rr - std::log(margin_lower(comparison, endpoint))) / se > z &&
            (std::log(margin_upper(comparison, endpoint)) - log_rr) / se > z) {
          endpoint_passed[endpoint] = 1;
          ++endpoint_successes[comparison * m + endpoint];
        }
      }
      if (count_endpoint_decision(typey, adseq, endpoint_passed, k)) {
        ++comparison_successes[comparison];
      } else {
        all_comparisons_pass = false;
      }
    }
    if (all_comparisons_pass) ++successes;
  }

  const int n_values = 2 + n_comparisons * m + n_comparisons;
  Rcpp::NumericVector out(n_values);
  out[0] = static_cast<double>(successes) / nsim;
  out[1] = successes;
  Rcpp::CharacterVector names(n_values);
  names[0] = "power";
  names[1] = "successes";
  int index = 2;
  for (int comparison = 0; comparison < n_comparisons; ++comparison) {
    for (int endpoint = 0; endpoint < m; ++endpoint) {
      out[index] = endpoint_successes[comparison * m + endpoint];
      names[index] = "comparison_endpoint_success_" +
        std::to_string(comparison + 1) + "_" + std::to_string(endpoint + 1);
      ++index;
    }
  }
  for (int comparison = 0; comparison < n_comparisons; ++comparison) {
    out[index] = comparison_successes[comparison];
    names[index] = "comparison_success_" + std::to_string(comparison + 1);
    ++index;
  }
  out.names() = names;
  return out;
}

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
arma::mat check_equivalence(const arma::ivec& typey,
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

  // Ensure at least k secondary endpoints pass. If the hierarchy contains no
  // secondary endpoints, all primary endpoints are co-primary and their
  // intersection is the final decision.
  int num_secondary_endpoints = accu(typey == 2);
  int num_secondary_successes = num_secondary_endpoints > 0 ?
    arma::accu(tbioq.cols(arma::find(typey == 2))) : 0;
  bool secondary_test_passed = num_secondary_endpoints == 0 ?
    (num_primary_successes == num_primary_endpoints) :
    (num_secondary_successes >= k);

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
                       const arma::ivec& typey,
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
                       const arma::ivec& typey,
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
                       const arma::ivec& typey,
                       const bool adseq,
                       const int k,
                       const int arm_seedT,
                       const int arm_seedR,
                       const double TART,
                       const double TARR,
                       const bool vareq){

   const arma::uword m = muT.n_elem;
   if (muR.n_elem != m || SigmaT.n_rows != m || SigmaT.n_cols != m ||
       SigmaR.n_rows != m || SigmaR.n_cols != m ||
       lequi_tol.n_elem != m || uequi_tol.n_elem != m ||
       alpha.n_elem != m) {
     Rcpp::stop("test_par_dom received inconsistent dimensions: muT=%d, muR=%d, SigmaT=%dx%d, SigmaR=%dx%d, lower=%d, upper=%d, alpha=%d.",
                static_cast<int>(muT.n_elem), static_cast<int>(muR.n_elem),
                static_cast<int>(SigmaT.n_rows), static_cast<int>(SigmaT.n_cols),
                static_cast<int>(SigmaR.n_rows), static_cast<int>(SigmaR.n_cols),
                static_cast<int>(lequi_tol.n_elem), static_cast<int>(uequi_tol.n_elem),
                static_cast<int>(alpha.n_elem));
   }

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
   mat yT = draw_mvnorm_rows(muT, SigmaT, n0);

   set_seed(arm_seedR);
   mat yR = draw_mvnorm_rows(muR, SigmaR, n1);

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
                       const arma::ivec& typey,
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
   mat yT = draw_mvnorm_rows(muT, SigmaT, n0);

   set_seed(arm_seedR);
   mat yR = draw_mvnorm_rows(muR, SigmaR, n1);

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
                               const arma::ivec& typey,
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
                                   const arma::ivec& typey,
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
                               const arma::ivec& typey,
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
                                   const arma::ivec& typey,
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
