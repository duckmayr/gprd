#' Plot Gaussian Process Regression Discontinuity
#'
#' @param x A numeric vector giving the observations of the explanatory variable
#' @param y A numeric vector giving the observations of the outcome
#' @param cutoff A numeric vector of length one giving the discontinutiy cutoff;
#'     the default is 0
#' @param from A numeric vector of length one giving the lowest x value
#'     to plot prediction for
#' @param to A numeric vector of length one giving the highest x value
#'     to plot prediction for
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
#' @param data_color A color to plot the data points in
#' @param line_color A color to plot the predictive mean line in
#' @param ci_color A color to plot the CI polygon in
#' @param main_title A character vector of length one giving a main title for
#'     the plot
#' @param xlab A character vector of length one giving a title for the x axis
#' @param ylab A character vector of length one giving a title for the y axis
#' @param ... Other arguments passed to \code{\link[graphics]{plot}}
#'
#' @return Returns \code{NULL} invisibly
#'
#' @importFrom stats qnorm
#' @importFrom graphics plot
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics polygon
#'
#' @export
plotrd <- function(x, y, cutoff = 0, from = min(x), to = max(x), # data
                   b_l = c(0, 0), b_r = c(0, 0), # beta prior means
                   B_l = diag(10, nrow = 2), # beta prior cov
                   B_r = diag(10, nrow = 2), # beta prior cov
                   output_sigma = 1, cov_sigma = 1, cov_ell = 1, # other hypers
                   ci_width = 0.95,
                   data_color = "#5f5f5f5f",
                   line_color = "black",
                   ci_color = "#87878787",
                   main_title = "", xlab = "", ylab = "", ...) {
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
    ## Make prediction case vectors
    left_test_cases  <- seq(from = from, to = cutoff, by = 0.01)
    right_test_cases <- seq(from = cutoff, to = to, by = 0.01)
    ## Farm out to C++ functions for computation
    left_preds  <- .gprd_predict(y[setdiff(left_cases, omit_cases)],
                                 x[setdiff(left_cases, omit_cases)],
                                 left_test_cases,
                                 B_l, b_l, cov_sigma, cov_ell, output_sigma)
    right_preds <- .gprd_predict(y[setdiff(right_cases, omit_cases)],
                                 x[setdiff(right_cases, omit_cases)],
                                 right_test_cases,
                                 B_r, b_r, cov_sigma, cov_ell, output_sigma)
    ## Setup vectors to plot
    q_high     <- 1 - ((1 - ci_width) / 2)
    q_low      <- 1 - q_high
    left_mean  <- left_preds[ , 1]
    left_low   <- qnorm(q_low,  mean = left_mean, sd = left_preds[ , 3])
    left_high  <- qnorm(q_high, mean = left_mean, sd = left_preds[ , 3])
    right_mean <- right_preds[ , 1]
    right_low  <- qnorm(q_low,  mean = right_mean, sd = right_preds[ , 3])
    right_high <- qnorm(q_high, mean = right_mean, sd = right_preds[ , 3])
    ## Set up plot
    margins <- "if"(main_title == "", c(3, 3, 1, 1), c(3, 3, 3, 1)) + 0.1
    opar <- par(mar = margins)
    on.exit(par(opar))
    lim <- range(c(left_low, right_low, left_high, right_high, y), na.rm = TRUE)
    plot(x = c(left_test_cases, right_test_cases), y = c(left_mean, right_mean),
         ylim = lim, type = "n", ...,
         main = main_title, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(side = 1, tick = FALSE, line = -0.75)
    axis(side = 2, tick = FALSE, line = -0.75)
    mtext(side = 1, text = xlab, line =  1.5)
    mtext(side = 2, text = ylab, line =  1.5)
    polygon(x = c(left_test_cases, rev(left_test_cases)),
            y = c(left_low, rev(left_high)),
            border = NA, col = ci_color)
    polygon(x = c(right_test_cases, rev(right_test_cases)),
            y = c(right_low, rev(right_high)),
            border = NA, col = ci_color)
    points(x, y, pch = 19, col = data_color)
    lines(left_test_cases, left_mean, col = line_color)
    lines(right_test_cases, right_mean, col = line_color)
    return(invisible(NULL))
}
