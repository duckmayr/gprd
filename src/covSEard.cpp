#include <RcppArmadillo.h>

arma::mat covSEard(const arma::mat& x, const double sigma_f,
                   const arma::vec& ell) {
    arma::uword n = x.n_rows;
    arma::uword p = ell.n_elem;
    double ssq  = sigma_f * sigma_f;
    arma::vec z = -0.5 * (1.0 / (ell % ell));
    arma::mat res(n, n);
    for ( arma::uword i = 0; i < n; ++i ) {
        res(i, i) = ssq;
    }
    for ( arma::uword j = 1; j < n; ++j ) {
        for ( arma::uword i = 0; i < (n-1); ++i ) {
            double elem = 0.0;
            for ( arma::uword k = 0; k < p; ++k ) {
                double diff = x(i, k) - x(j, k);
                elem += (z[k] * diff * diff);
            }
            elem = ssq * std::exp(elem);
            res(i, j) = elem;
            res(j, i) = elem;
        }
    }
    return res;
}

arma::mat covSEard(const arma::mat& x1, const arma::mat& x2,
                   const double sigma_f,
                   const arma::vec& ell) {
    arma::uword n = x1.n_rows;
    arma::uword m = x2.n_rows;
    arma::uword p = ell.n_elem;
    double ssq = sigma_f * sigma_f;
    arma::vec z = -0.5 * (1.0 / (ell % ell));
    arma::mat res(n, m);
    for ( arma::uword j = 0; j < m; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            double elem = 0.0;
            for ( arma::uword k = 0; k < p; ++k ) {
                double diff = x1(i, k) - x2(j, k);
                elem += (z[k] * diff * diff);
            }
            res(i, j) = ssq * std::exp(elem);
        }
    }
    return res;
}
