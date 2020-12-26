#include <RcppArmadillo.h>

// [[Rcpp::export(.squared_distance)]]
arma::cube squared_distance(const arma::mat& x) {
    arma::uword n = x.n_rows;
    arma::uword m = x.n_cols;
    arma::cube res(n, n, m);
    for ( arma::uword s = 0; s < m; ++s ) {
        for ( arma::uword j = 0; j < n; ++j ) {
            res(j, j, s) = 0.0;
            for ( arma::uword i = j+1; i < n; ++ i ) {
                double diff = x(i, s) - x(j, s);
                double z = diff * diff;
                res(i, j, s) = z;
                res(j, i, s) = z;
            }
        }
    }
    return res;
}

// [[Rcpp::export(.gradient)]]
arma::vec gradient(const arma::vec& hypers,
                   const arma::cube& K0,
                   const arma::mat& Q0,
                   const arma::vec& M) {
    /* K0 here is squared distance on each dimension.
     * K1 will be K / s_f^2.
     */
    // Bookkeeping variables
    arma::uword n = K0.n_rows;
    arma::uword p = K0.n_slices;
    // Hyperparameters
    double sy    = hypers[0];
    double sy_sq = sy * sy;
    double sf    = hypers[1];
    double sf_sq = sf * sf;
    arma::vec ell    = hypers.subvec(2, p+1);
    arma::vec ell_sq = ell % ell;
    arma::vec ell_3  = ell % ell_sq;
    // Generate various stages of computation of K
    arma::mat K1(n, n);
    for ( arma::uword j = 0; j < n; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            arma::vec tunnel = K0.tube(i, j);
            K1(i, j) = std::exp(-0.5 * arma::sum(tunnel / ell_sq));
        }
    }
    arma::mat K     = sf_sq * K1;
    arma::mat I     = arma::eye(n, n);
    arma::mat Qinv  = arma::inv_sympd(Q0 + K + (sy_sq * I));
    arma::mat QinvM = Qinv * M;
    arma::vec res(hypers.n_elem);
    arma::mat tmp = 2.0 * sy * I;
    res[0] = 0.5 * (arma::as_scalar(QinvM.t() * tmp * QinvM)
                        - arma::trace(Qinv));
    tmp = 2.0 * sf * K1;
    res[1] = 0.5 * (arma::as_scalar(QinvM.t() * tmp * QinvM)
                        - arma::trace(Qinv * tmp));
    for ( arma::uword k = 0; k < p; ++k ) {
        tmp =  K % (K0.slice(k) / ell_3[k]);
        res[k+2] = 0.5 * (arma::as_scalar(QinvM.t() * tmp * QinvM)
                              - arma::trace(Qinv * tmp));
    }
    return res;
}

// [[Rcpp::export(.log_marginal_likelihood)]]
double log_marginal_likelihood(const arma::vec& hypers,
                               const arma::cube& K0,
                               const arma::mat& Q0,
                               const arma::vec& M) {
    arma::uword n = K0.n_rows;
    arma::uword p = K0.n_slices;
    arma::mat K(n, n);
    double sf_sq = hypers[1];
    arma::vec ell_sq = hypers.subvec(2, p+1);
    for ( arma::uword j = 0; j < n; ++j ) {
        K(j, j) = sf_sq;
        for ( arma::uword i = j+1; i < n; ++i ) {
            arma::vec tunnel = K0.tube(i, j);
            double z = sf_sq * std::exp(-0.5 * arma::sum(tunnel / ell_sq));
            K(i, j) = z;
            K(j, i) = z;
        }
    }
    arma::mat Q  = Q0 + K + (hypers[0] * arma::eye(n, n));
    arma::mat U  = arma::chol(Q);
    double ld    =  2.0 * arma::sum(arma::log(U.diag()));
    U            = arma::inv(arma::trimatu(U));
    double mTQim = arma::as_scalar(M.t() * U * U.t() * M);
    return -0.5 * ( mTQim + ld + n * std::log(M_2PI) );
}
