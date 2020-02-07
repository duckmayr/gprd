#' Gaussian Process Regression for Regression Discontinuity
#'
#' @param x A numeric vector giving the observations of the explanatory variable
#' @param y A numeric vector giving the observations of the outcome
#' @param cutoff A numeric vector of length one giving the discontinutiy cutoff;
#'     the default is 0.
#' @param estimator A character vector of length one giving the estimation
#'     strategy; should be one of "global" (one unknown function for all
#'     observations, with dummy variable added to indicate treatment) or
#'     "local" (different unknown functions for observations to the left of
#'     the cutoff and to the right). The default is "local".
#' @param hypers A list giving the hyperparameters of the model. If
#'     \code{estimator = "local"}, should be a list of length two, where each
#'     element is itself a list giving the hyperparameters for the functions to
#'     each side of the cutoff. The list(s) of hypers should have elements
#'     "b" (a numeric vector giving the coefficients' prior mean),
#'     "B" (a numeric matrix giving the coefficients' prior covariance),
#'     "sigma_y" (a numeric vector of length one giving the output variance),
#'     "sigma_f" (a numeric vector of length one giving the function variance),
#'     and "ell" (a numeric vector giving the length scale(s)).
#' @param ci_width A numeric vector of length one between 0 and 1 giving
#'     the width of the confidence interval for the treatment effect.
#'
#' @return A list of length five with class \code{gprd} and elements
#'     \describe{
#'         \item{tau_mean}{The treatment effect estimate}
#'         \item{tau_ci}{A named numeric ector of length two giving the
#'                       requested confidence interval (the names are the
#'                       confidence interval quantiles)}
#'         \item{cutoff}{A numeric vector of length one giving the cutoff}
#'         \item{estimator}{A character vector of length one giving the
#'                          estimation strategy}
#'         \item{f}{A list where each element is an object of class \code{gpr}
#'                  with information on the posterior over the relevant
#'                  function; if \code{estimator = "global"}, the list will be
#'                  of length one, giving the posterior over the global mapping
#'                  from \code{x} to \code{y}, while if
#'                  \code{estimator = "local"}, the list will be of length two
#'                  with elements \code{f_l} (giving the posterior over the
#'                  function to the left of the cutoff) and \code{f_r} (giving
#'                  the posterior over the function to the right of the cutoff)}
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
#'                            \code{sigma_f}, giving the scale factor for the
#'                            GP prior;
#'                            \code{ell}, giving the length scale for the GP
#'                            prior;
#'                            and \code{sigma_y}, giving the outcome observation
#'                            noise (the likelihood hyper).}
#' }
#'
#' @importFrom stats qnorm
#' @seealso plot.gprd
#'
#' @export
gprd <- function(x, y, cutoff = 0, estimator = "local",
                 hypers = list(list(b = c(0, 0), B = diag(10, nrow = 2),
                                    sigma_y = 1, sigma_f = 1, ell = 1),
                               list(b = c(0, 0), B = diag(10, nrow = 2),
                                    sigma_y = 1, sigma_f = 1, ell = 1)),
                 ci_width = 0.95) {
    ## Sanity checks
    if ( !is.vector(x) | !is.numeric(x) ) {
        stop("x should be a numeric vector.")
    }
    if ( !is.vector(y) | !is.numeric(y) ) {
        stop("y should be a numeric vector.")
    }
    if ( length(x) != length(y) ) {
        stop("x and y should be of the same length.")
    }
    if ( ci_width > 1 | ci_width < 0 ) {
        stop("ci_width must be between 0 and 1.")
    }
    ## Find which cases belong to treatment, control, or are missing
    omit_cases  <- which(is.na(x) | is.na(y))
    left_cases  <- which(x < cutoff)
    right_cases <- which(x > cutoff)
    ## Farm out to C++ functions for computation
    if ( estimator == "local" ) {
        f_l  <- .gprd(y[setdiff(left_cases, omit_cases)],
                      matrix(x[setdiff(left_cases, omit_cases)]),
                      matrix(cutoff), # only need to predict at the cutoff
                      hypers[[1]][["B"]],
                      hypers[[1]][["b"]],
                      hypers[[1]][["sigma_f"]],
                      hypers[[1]][["ell"]],
                      hypers[[1]][["sigma_y"]])
        f_r <- .gprd(y[setdiff(right_cases, omit_cases)],
                     matrix(x[setdiff(right_cases, omit_cases)]),
                     matrix(cutoff), # only need to predict at the cutoff
                     hypers[[2]][["B"]],
                     hypers[[2]][["b"]],
                     hypers[[2]][["sigma_f"]],
                     hypers[[2]][["ell"]],
                     hypers[[2]][["sigma_y"]])
        f_res <- list(f_l = f_l, f_r = f_r)
    } else if ( estimator == "global" ) {
        X     <- x[setdiff(seq_along(x), omit_cases)]
        X     <- cbind(X, X > cutoff)
        cases <- matrix(c(cutoff, cutoff, 1, 0), ncol = 2)
        f     <- .gprd(y[setdiff(seq_along(y), omit_cases)],
                       X,
                       cases,
                       hypers[["B"]],
                       hypers[["b"]],
                       hypers[["sigma_f"]],
                       hypers[["ell"]],
                       hypers[["sigma_y"]])
        f_res <- list(f)
    } else {
        stop("estimator should be one of global or local.")
    }
    ## Then get the estimate and the CI
    q_high   <- 1 - ((1 - ci_width) / 2)
    q_low    <- 1 - q_high
    if ( estimator == "local" ) {
        tau_mean <- f_r$predictive_mean[1,1] - f_l$predictive_mean[1,1]
        tau_sd   <- sqrt(f_r$predictive_var[1,1] + f_l$predictive_var[1,1])
    } else {
        tau_mean <- f$predictive_mean[1,1] - f$predictive_mean[2,1]
        tau_sd   <- sqrt(f$predictive_var[1,1] + f$predictive_var[2,1])
    }
    tau_low  <- qnorm(p = q_low,  mean = tau_mean, sd = tau_sd)
    tau_high <- qnorm(p = q_high, mean = tau_mean, sd = tau_sd)
    tau_ci   <- c(tau_low, tau_high)
    names(tau_ci) <- c(q_low, q_high)
    res <- list(tau_mean = tau_mean, tau_ci = tau_ci,
                cutoff = cutoff, estimator = estimator, f = f_res)
    class(res) <- c("gprd", class(res))
    return(res)
}
