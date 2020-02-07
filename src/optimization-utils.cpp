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
    arma::uword n = K0.n_rows;
    arma::uword p = K0.n_slices;
    arma::cube K1(n, n, p);
    arma::mat K2(n, n); // Getting expensive memory allocation wise
    arma::mat K(n, n);  // Faster to redo calculations?
    arma::vec ell_sq = hypers.subvec(2, p+1) % hypers.subvec(2, p+1);
    double sf = hypers[1];
    double sf_sq = sf * sf;
    for ( arma::uword j = 0; j < n; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            arma::vec tunnel = K0.tube(i, j);
            tunnel /= ell_sq;
            K1.tube(i, j) = tunnel;
            double z = arma::sum(tunnel);
            K2(i, j) = z;
            K(i, j)  = sf_sq * std::exp(-0.5 * z);
        }
    }
    arma::mat I = arma::eye(n, n);
    arma::mat Qinv = arma::inv(Q0 + K + (hypers[0] * I));
    arma::mat QinvM = Qinv * M;
    arma::vec res(hypers.n_elem);
    arma::mat tmp = 2.0 * hypers[0] * I;
    res[0] = 0.5 * (arma::as_scalar(QinvM.t() * QinvM) - arma::trace(Qinv));
    tmp = 2.0 * sf * arma::exp(-0.5 * K2);
    res[1] = 0.5 * (arma::as_scalar(QinvM.t() * tmp * QinvM)
                        - arma::trace(Qinv * tmp));
    for ( arma::uword k = 0; k < p; ++k ) {
        tmp = K % K1.slice(k);
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
    double sf = hypers[1];
    double sf_sq = sf * sf;
    arma::vec ell_sq = hypers.subvec(2, p+1) % hypers.subvec(2, p+1);
    for ( arma::uword j = 0; j < n; ++j ) {
        K(j, j) = sf_sq;
        for ( arma::uword i = j+1; i < n; ++i ) {
            arma::vec tunnel = K0.tube(i, j);
            double z = sf_sq * std::exp(-0.5 * arma::sum(tunnel / ell_sq));
            K(i, j) = z;
            K(j, i) = z;
        }
    }
    arma::mat I = arma::eye(n, n);
    arma::mat Q = Q0 + K + (hypers[0] * I);
    arma::mat U;
    bool able_to_decompose = arma::chol(U, Q);
    if ( !able_to_decompose ) {
        Rcpp::stop("Optimization unstable; try another optimization method.");
    }
    arma::mat w = arma::solve(arma::trimatu(U),
                              arma::solve(arma::trimatl(U.t()), M));
    double tmp = arma::as_scalar(M.t() * w);
    double ld;
    double sign;
    arma::log_det(ld, sign, Q);
    return -0.5 * ( tmp + ld + n * std::log(M_2PI) );
}
