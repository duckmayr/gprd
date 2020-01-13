#include <RcppArmadillo.h>

arma::mat covSEiso(const arma::vec& x, const double sigma, const double ell) {
    arma::uword n = x.n_elem;
    double ssq = sigma * sigma;
    double z   = -0.5 * (1.0 / (ell * ell));
    arma::mat res(n, n);
    for ( arma::uword i = 0; i < n; ++i ) {
        res(i, i) = ssq;
    }
    for ( arma::uword j = 1; j < n; ++j ) {
        for ( arma::uword i = 0; i < (n-1); ++i ) {
            double diff = x[i] - x[j];
            double elem = ssq * std::exp(z * diff * diff);
            res(i, j) = elem;
            res(j, i) = elem;
        }
    }
    return res;
}

arma::mat covSEiso(const arma::vec& x1, const arma::vec& x2,
                   const double sigma, const double ell) {
    arma::uword n = x1.n_elem;
    arma::uword m = x2.n_elem;
    double ssq = sigma * sigma;
    double z   = -0.5 * (1.0 / (ell * ell));
    arma::mat res(n, m);
    for ( arma::uword j = 0; j < m; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            double diff = x1[i] - x2[j];
            res(i, j) = ssq * std::exp(z * diff * diff);
        }
    }
    return res;
}
