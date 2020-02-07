#ifndef GPRD_H
#define GPRD_H

#include <RcppArmadillo.h>

arma::mat covSEiso(const arma::vec& x, const double sigma, const double ell);
arma::mat covSEiso(const arma::vec& x1, const arma::vec& x2,
                   const double sigma, const double ell);

arma::mat covSEard(const arma::mat& x, const double sigma_f,
                   const arma::vec& ell);
arma::mat covSEard(const arma::mat& x1, const arma::mat& x2,
                   const double sigma_f,
                   const arma::vec& ell);

inline arma::vec ones_vec(arma::uword n) {
    return arma::ones<arma::vec>(n);
}

#endif
