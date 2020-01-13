#ifndef GPRD_H
#define GPRD_H

#include <RcppArmadillo.h>

arma::mat covSEiso(const arma::vec& x, const double sigma, const double ell);
arma::mat covSEiso(const arma::vec& x1, const arma::vec& x2,
                   const double sigma, const double ell);

#endif
