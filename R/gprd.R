#' Gaussian Process Regression for Regression Discontinuity
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
#' @return A list of length five with class \code{gprd} and elements
#'     \describe{
#'         \item{tau_mean}{The treatment effect estimate}
#'         \item{tau_ci}{A named numeric ector of length two giving the
#'                       requested confidence interval (the names are the
#'                       confidence interval quantiles)}
#'         \item{cutoff}{A numeric vector of length one giving the cutoff}
#'         \item{f_l}{An object of class \code{gpr} with information on
#'                    the posterior over the function to the left of the cutoff}
#'         \item{f_r}{An object of class \code{gpr} with information on
#'                    the posterior over the function to the right of the cutoff}
#'     }
#'
#'     An object of class \code{gpr} is a list of length eight with the following
#' elements:
#'
#' \describe{
#'     \item{predictive_mean}{The mean of the predictive distribution at test
#'                            cases \code{test_input}; in this case, the cutoff}
#'     \item{predictive_var}{The variance of the predictive distribution at test
#'                           cases \code{test_input}; in this case, the cutoff}
#'     \item{beta_bar}{The mean of the posterior of the linear mean coefficients}
#'     \item{training_input}{The input variable obervations for training cases}
#'     \item{training_outcomes}{The outcome observations for training cases}
#'     \item{test_input}{The input variable observations predictions were
#'                       provided for}
#'     \item{prior_covariance}{The covariance function evaluated for the
#'                             training inputs, with observation noise added}
#'     \item{prior_cov_inverse}{The inverse of \code{prior_covariance};
#'                              having this available is useful for generating
#'                              other predictions}
#'     \item{hyperparameters}{The hyperparameters used for the GP regression,
#'                            provided in a list of length five:
#'                            \code{beta_prior_cov}, giving the
#'                            beta prior covariance;
#'                            \code{beta_prior_mean}, giving the beta prior
#'                            mean;
#'                            \code{scale_factor}, giving the scale factor for
#'                            the GP prior;
#'                            \code{length_scale}, giving the length scale for
#'                            the GP prior;
#'                            and \code{outcome_sigma}, giving the outcome
#'                            observation noise (the likelihood hyper).}
#' }
#'
#' @importFrom stats qnorm
#' @seealso plot.gprd
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
    f_l  <- .gprd(y[setdiff(left_cases, omit_cases)],
                  x[setdiff(left_cases, omit_cases)],
                  cutoff, # only need to predict at the cutoff
                  B_l, b_l, cov_sigma, cov_ell, output_sigma)
    f_r <- .gprd(y[setdiff(right_cases, omit_cases)],
                 x[setdiff(right_cases, omit_cases)],
                 cutoff, # only need to predict at the cutoff
                 B_r, b_r, cov_sigma, cov_ell, output_sigma)
    ## Then get the estimate and the CI
    q_high   <- 1 - ((1 - ci_width) / 2)
    q_low    <- 1 - q_high
    tau_mean <- f_r$predictive_mean[1,1] - f_l$predictive_mean[1,1]
    tau_sd   <- sqrt(f_r$predictive_var[1,1] + f_l$predictive_var[1,1])
    tau_low  <- qnorm(p = q_low,  mean = tau_mean, sd = tau_sd)
    tau_high <- qnorm(p = q_high, mean = tau_mean, sd = tau_sd)
    res <- list(tau_mean = tau_mean, tau_ci = c(tau_low, tau_high),
                cutoff = cutoff, f_l = f_l, f_r = f_r)
    class(res) <- c("gprd", class(res))
    return(res)
}
