#include "gprd.h"

inline arma::vec ones_vec(arma::uword n) {
    return arma::ones<arma::vec>(n);
}

inline arma::vec zeros_vec(arma::uword n) {
    return arma::zeros<arma::vec>(n);
}

// [[Rcpp::export(.gprd_predict)]]
arma::mat predict(const arma::vec& outcome,
                  const arma::vec& training_cases,
                  const arma::vec& test_cases,
                  const arma::mat& B, // beta prior covariance
                  const arma::vec& b, // beta prior mean
                  const double cov_sigma,
                  const double cov_ell,
                  const double obs_sigma) {
    // Bookkeeping variables
    arma::uword n = outcome.n_elem;
    arma::uword m = test_cases.n_elem;

    // Covariance matrices
    arma::mat noise = obs_sigma * arma::eye(n, n);
    arma::mat Ky    = covSEiso(training_cases, cov_sigma, cov_ell) + noise;
    arma::mat Ky_i  = Ky.i();
    arma::mat Kstar = covSEiso(training_cases, test_cases, cov_sigma, cov_ell);
    arma::mat Kstar_t   = Kstar.t();
    arma::mat Kstarstar = covSEiso(test_cases, cov_sigma, cov_ell);

    // Add ones for intercept
    arma::mat H = arma::join_horiz(ones_vec(n), training_cases);
    arma::mat Hstar = arma::join_horiz(ones_vec(m), test_cases);

    // Objects we will reuse
    // The notation is from Rasmussen & Williams 2006, eq. 2.41,
    // but note that we use rows as observations instead of columns
    arma::mat H_Ky_i = H.t() * Ky_i; // 2 x n
    arma::mat R = Hstar.t() - (H_Ky_i * Kstar); // 2 x m
    arma::mat B_i = B.i();
    arma::mat beta_bar_start = B_i + (H_Ky_i * H); // 2 x 2
    // then beta_bar will be 2 x 1
    arma::vec beta_bar = beta_bar_start.i() * ((H_Ky_i * outcome) + (B_i * b));

    // Results
    arma::vec mX = outcome - (H * beta_bar); // n x 1
    arma::vec predictive_mean = Hstar * beta_bar + Kstar_t * Ky_i * mX; // m x 1
    arma::mat cov_fstar = Kstarstar - (Kstar_t * Ky_i * Kstar); // m x m
    arma::mat tmp = B_i + (H_Ky_i * H); // 2 x 2
    arma::mat predictive_cov = cov_fstar + R.t() * tmp.i() * R; // m x m
    arma::vec predictive_var = predictive_cov.diag(); // m x 1
    arma::vec predictive_sds = arma::sqrt(predictive_var); // m x 1
    arma::mat res(m, 3);
    res.col(0) = predictive_mean;
    res.col(1) = predictive_var;
    res.col(2) = predictive_sds;
    return res;
}
