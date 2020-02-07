#include "gprd.h"

// [[Rcpp::export(.gprd_predict)]]
Rcpp::List predict(const arma::vec& outcome,
                   const arma::mat& training_cases,
                   const arma::mat& test_cases,
                   const arma::mat& B, // beta prior covariance
                   const arma::vec& b, // beta prior mean
                   const double cov_sigma,
                   const arma::vec& cov_ell,
                   const double obs_sigma,
                   const arma::mat& Ky,
                   const arma::mat& Ky_i,
                   const arma::vec& beta_bar) {
    // Bookkeeping variables
    arma::uword n = outcome.n_elem;
    arma::uword m = test_cases.n_rows;

    // Covariance matrices
    arma::mat covX  = training_cases;
    covX.shed_col(0);
    arma::mat Kstar = covSEard(covX, test_cases, cov_sigma, cov_ell);
    arma::mat Kstar_t   = Kstar.t();
    arma::mat Kstarstar = covSEard(test_cases, cov_sigma, cov_ell);

    // Add ones for intercept
    arma::mat H = training_cases;
    arma::mat Hstar = arma::join_horiz(ones_vec(m), test_cases);

    // Objects we will reuse
    // The notation is from Rasmussen & Williams 2006, eq. 2.41,
    // but note that we use rows as observations instead of columns
    arma::mat H_Ky_i = H.t() * Ky_i; // 2 x n
    arma::mat R = Hstar.t() - (H_Ky_i * Kstar); // 2 x m
    arma::mat B_i = B.i();

    // Results
    arma::vec mX = outcome - (H * beta_bar); // n x 1
    arma::vec predictive_mean = Hstar * beta_bar + Kstar_t * Ky_i * mX; // m x 1
    arma::mat cov_fstar = Kstarstar - (Kstar_t * Ky_i * Kstar); // m x m
    arma::mat tmp = B_i + (H_Ky_i * H); // 2 x 2
    arma::mat predictive_cov = cov_fstar + R.t() * tmp.i() * R; // m x m
    arma::vec predictive_var = predictive_cov.diag(); // m x 1
    Rcpp::List hyperparameters = Rcpp::List::create(
            Rcpp::Named("beta_prior_cov")    = B,
            Rcpp::Named("beta_prior_mean")   = b,
            Rcpp::Named("scale_factor")      = cov_sigma,
            Rcpp::Named("length_scale")      = cov_ell,
            Rcpp::Named("outcome_sigma")     = obs_sigma
    );
    Rcpp::List res = Rcpp::List::create(
            Rcpp::Named("predictive_mean")   = predictive_mean,
            Rcpp::Named("predictive_var")    = predictive_var,
            Rcpp::Named("beta_bar")          = beta_bar,
            Rcpp::Named("training_input")    = H,
            Rcpp::Named("training_outcomes") = outcome,
            Rcpp::Named("test_input")        = Hstar,
            Rcpp::Named("prior_covariance")  = Ky,
            Rcpp::Named("prior_cov_inverse") = Ky_i,
            Rcpp::Named("hyperparameters")   = hyperparameters
    );
    res.attr("class") = "gpr";
    return res;
}
