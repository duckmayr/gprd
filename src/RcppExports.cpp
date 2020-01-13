// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// predict
arma::mat predict(const arma::vec& outcome, const arma::vec& training_cases, const arma::vec& test_cases, const arma::mat& B, const arma::vec& b, const double cov_sigma, const double cov_ell, const double obs_sigma);
RcppExport SEXP _gprd_predict(SEXP outcomeSEXP, SEXP training_casesSEXP, SEXP test_casesSEXP, SEXP BSEXP, SEXP bSEXP, SEXP cov_sigmaSEXP, SEXP cov_ellSEXP, SEXP obs_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type training_cases(training_casesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type test_cases(test_casesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type cov_sigma(cov_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type cov_ell(cov_ellSEXP);
    Rcpp::traits::input_parameter< const double >::type obs_sigma(obs_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(predict(outcome, training_cases, test_cases, B, b, cov_sigma, cov_ell, obs_sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gprd_predict", (DL_FUNC) &_gprd_predict, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_gprd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
