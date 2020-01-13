#' Gaussian Process Regression Discontinuity
#'
#' @param x A numeric vector giving the observations of the explanatory variable
#' @param y A numeric vector giving the observations of the outcome
#' @param cutoff A numeric vector of length one giving the discontinutiy cutoff;
#'     the default is 0.
#' @param b_l A numeric vector of length two giving the prior means for the
#'     intercept and slope to the left of the cutoff
#' @param b_r A numeric vector of length two giving the prior means for the
#'     intercept and slope to the right of the cutoff
#' @param B_l A numeric matrix with 2 rows and 2 columns giving the prior
#'     covariance for the intercept and slope to the left of the cutoff
#' @param B_r A numeric matrix with 2 rows and 2 columns giving the prior
#'     covariance for the intercept and slope to the right of the cutoff
#' @param output_sigma A numeric vector of length one; a hyperparameter for
#'     the observation noise
#' @param cov_sigma A numeric vector of length one; a hyperparameter for the
#'     signal noise in the unobserved function
#' @param cov_ell A numeric vector of length one; a hyperparameter for the
#'     characteristic length scale for the unobserved function
#' @param ci_width A numeric vector of length one between 0 and 1 giving
#'     the width of the confidence interval for tau
#'
#' @return A named vector of length three with elements
#'     \describe{
#'         \item{tau}{The treatment effect estimate}
#'         \item{tau_low}{The lower bound of the confidence interval}
#'         \item{tau_high}{The upper bound of the confidence interval}
#'     }
#'
#' @importFrom stats qnorm
#'
#' @export
gprd <- function(x, y, cutoff = 0, # data
                 b_l = c(0, 0), b_r = c(0, 0), # prior means
                 B_l = diag(10, nrow = 2), B_r = diag(10, nrow = 2), # prior cov
                 output_sigma = 1, cov_sigma = 1, cov_ell = 1, # other hypers
                 ci_width = 0.95) {
    ## Sanity checks
    if ( !is.vector(x) | !is.numeric(x) ) {
        stop("x should be a numeric vector.")
    }
    if ( !is.vector(y) | !is.numeric(y) ) {
        stop("x should be a numeric vector.")
    }
    if ( length(x) != length(y) ) {
        stop("x and y should be of the same length.")
    }
    ## Find which cases belong to treatment, control, or are missing
    omit_cases  <- which(is.na(x) | is.na(y))
    left_cases  <- which(x < cutoff)
    right_cases <- which(x > cutoff)
    ## Farm out to C++ functions for computation
    left_preds  <- .gprd_predict(y[setdiff(left_cases, omit_cases)],
                                 x[setdiff(left_cases, omit_cases)],
                                 cutoff, # only need to predict at the cutoff
                                 B_l, b_l, cov_sigma, cov_ell, output_sigma)
    right_preds <- .gprd_predict(y[setdiff(right_cases, omit_cases)],
                                 x[setdiff(right_cases, omit_cases)],
                                 cutoff, # only need to predict at the cutoff
                                 B_r, b_r, cov_sigma, cov_ell, output_sigma)
    ## Then get the estimate and the CI
    q_high   <- 1 - ((1 - ci_width) / 2)
    q_low    <- 1 - q_high
    tau_mean <- right_preds[1, 1] - left_preds[1, 1]
    tau_sd   <- sqrt(right_preds[1, 2] + left_preds[1, 2])
    tau_low  <- qnorm(p = q_low,  mean = tau_mean, sd = tau_sd)
    tau_high <- qnorm(p = q_high, mean = tau_mean, sd = tau_sd)
    return(c(tau = tau_mean, tau_low = tau_low, tau_high = tau_high))
}
