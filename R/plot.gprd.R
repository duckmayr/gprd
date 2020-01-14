#' Plot Gaussian Process Regression for Regression Discontinuity
#'
#' @param x A object of class \code{\link{gprd}}
#' @param from A numeric vector of length one giving the lowest x value
#'     to plot prediction for; if \code{NULL}, only values in the training
#'     data and at the cutoff are plotted
#' @param to A numeric vector of length one giving the highest x value
#'     to plot prediction for; if \code{NULL}, only values in the training
#'     data and at the cutoff are plotted
#' @param n_points An integer vector of length one giving the number of
#'     prediction points to plot; if \code{NULL} and \code{from} and \code{to}
#'     are given, \code{n_points = length(seq(from = from, to = to, by = 0.01))}
#' @param ci_width A numeric vector of length one between 0 and 1 giving
#'     the width of the confidence interval for tau
#' @param data_color A color to plot the data points in
#' @param line_color A color to plot the predictive mean line in
#' @param ci_color A color to plot the CI polygon in
#' @param plot_cutoff A logical vector of length one; if \code{TRUE}, a dashed
#'     vertical line (\code{lty = 2}) marks the cutoff; default is \code{TRUE}
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
#' @seealso gprd
#'
#' @export
plot.gprd <- function(x,
                      from = min(x$f_l$training_input),
                      to = max(x$f_r$training_input), # data
                      n_points = NULL, ci_width = 0.95,
                      data_color = "#1c1c1c1c",
                      line_color = "black",
                      ci_color = "#87878787",
                      plot_cutoff = TRUE,
                      main_title = "", xlab = "", ylab = "", ...) {
    ## Sanity checks
    is_gprd <- inherits(x, "gprd")
    has_both_fs <- all(c("f_r", "f_l") %in% names(x))
    both_fs_are_gpr_objects <- all(sapply(x[c("f_r", "f_l")], inherits, "gpr"))
    if ( !is_gprd | !has_both_fs | !both_fs_are_gpr_objects ) {
        stop("x should be a gprd object returned from gprd().")
    }
    if ( ci_width > 1 | ci_width < 0 ) {
        stop("ci_width must be between 0 and 1.")
    }
    ## Predict at new points
    if ( !is.null(from) & !is.null(to) ) {
        if ( is.null(n_points) ) {
            left_test  <- seq(from = from, to = x$cutoff, by = 0.01)
            right_test <- seq(from = x$cutoff, to = to, by = 0.01)
        } else {
            left_test  <- seq(from = from, to = x$cutoff, length.out = n_points)
            right_test <- seq(from = x$cutoff, to = to, length.out = n_points)
        }
        ## Make prediction case vectors
        ## Farm out to C++ functions for computation
        left_preds  <- .gprd_predict(x$f_l$training_outcomes,
                                     x$f_l$training_input,
                                     left_test,
                                     x$f_l$hyperparameters$beta_prior_cov,
                                     x$f_l$hyperparameters$beta_prior_mean,
                                     x$f_l$hyperparameters$scale_factor,
                                     x$f_l$hyperparameters$length_scale,
                                     x$f_l$hyperparameters$outcome_sigma,
                                     x$f_l$prior_covariance,
                                     x$f_l$prior_cov_inverse,
                                     x$f_l$beta_bar)
        right_preds <- .gprd_predict(x$f_r$training_outcomes,
                                     x$f_r$training_input,
                                     right_test,
                                     x$f_r$hyperparameters$beta_prior_cov,
                                     x$f_r$hyperparameters$beta_prior_mean,
                                     x$f_r$hyperparameters$scale_factor,
                                     x$f_r$hyperparameters$length_scale,
                                     x$f_r$hyperparameters$outcome_sigma,
                                     x$f_r$prior_covariance,
                                     x$f_r$prior_cov_inverse,
                                     x$f_r$beta_bar)
    } else {
        left_test   <- c(x$f_l$training_input[,2], x$cutoff)
        right_test  <- c(x$cutoff, x$f_r$trainint_input[,2])
        left_preds  <- x$f_l
        right_preds <- x$f_r
        left_preds$predictive_mean <- c(left_preds$training_outcomes,
                                        left_preds$predictive_mean)
        left_preds$predictive_var  <- c(diag(left_pred$prior_covariance),
                                        left_preds$predictive_var)
        right_preds$predictive_mean <- c(right_preds$predictive_mean,
                                         right_preds$training_outcomes)
        right_preds$predictive_var  <- c(right_preds$predictive_var,
                                         diag(right_preds$prior_covariance))
    }
    ## Setup vectors to plot
    q_high     <- 1 - ((1 - ci_width) / 2)
    q_low      <- 1 - q_high
    left_mean  <- left_preds$predictive_mean
    sds        <- sqrt(left_preds$predictive_var)
    left_low   <- qnorm(q_low,  mean = left_mean, sd = sds)
    left_high  <- qnorm(q_high, mean = left_mean, sd = sds)
    right_mean <- right_preds$predictive_mean
    sds        <- sqrt(right_preds$predictive_var)
    right_low  <- qnorm(q_low,  mean = right_mean, sd = sds)
    right_high <- qnorm(q_high, mean = right_mean, sd = sds)
    ## Set up plot
    margins <- "if"(main_title == "", c(3, 3, 1, 1), c(3, 3, 3, 1)) + 0.1
    opar <- par(mar = margins)
    on.exit(par(opar))
    lim <- range(c(left_low, right_low, left_high, right_high,
                   x$f_l$training_outcomes, x$f_r$training_outcomes),
                 na.rm = TRUE)
    plot(x = c(left_test, right_test), y = c(left_mean, right_mean),
         ylim = lim, type = "n", ...,
         main = main_title, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(side = 1, tick = FALSE, line = -0.75)
    axis(side = 2, tick = FALSE, line = -0.75)
    mtext(side = 1, text = xlab, line =  1.5)
    mtext(side = 2, text = ylab, line =  1.5)
    points(x$f_l$training_input[ , 2], x$f_l$training_outcomes,
           pch = 19, col = data_color)
    points(x$f_r$training_input[ , 2], x$f_r$training_outcomes,
           pch = 19, col = data_color)
    polygon(x = c(left_test, rev(left_test)),
            y = c(left_low, rev(left_high)),
            border = NA, col = ci_color)
    polygon(x = c(right_test, rev(right_test)),
            y = c(right_low, rev(right_high)),
            border = NA, col = ci_color)
    lines(left_test, left_mean, col = line_color)
    lines(right_test, right_mean, col = line_color)
    if ( plot_cutoff ) {
        abline(v = x$cutoff, lty = 2)
    }
    return(invisible(NULL))
}
