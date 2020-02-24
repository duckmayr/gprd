#' Optimize hyperparameters for GPRD
#'
#' @param x A numeric vector; the forcing variable
#' @param y A numeric vector; the outcome variable
#' @param cutoff A numeric vector of length one; the treatment cutoff
#'     (default is 0)
#' @param estimator A character vector of length one; should be one of "global"
#'     or "piecewise". The default is "global".
#' @param b A numeric vector giving the prior mean for the coefficients for the
#'     mean function. If \code{NULL} (the default), all coefficients have a
#'     prior mean of 0.
#' @param B A numeric matrix giving the prior covariance for the coefficients
#'     for the mean function. If \code{NULL} (the default), the prior
#'     covariance for the coefficients is a diagonal matrix, with the diagonal
#'     elements set equal to 10.
#' @param method A character vector of length one giving the optimization
#'     routine to use; default is "CG", the conjugate gradient method.
#'     See \code{\link[stats]{optim}} for other options.
#' @param control A list of control parameters for optim
#' @param ... Other arguments passed to \code{\link[stats]{optim}}
#'
#' @importFrom stats qnorm
#' @importFrom stats optim
#' @seealso gprd
#'
#' @export
optimize_hyperparameters <- function(x, y, cutoff = 0, estimator = "piecewise",
                                     b = NULL, B = NULL, method = "CG",
                                     control = list(type = 2), ...) {
    if ( length(x) != length(y) ) {
        stop("x and y have a different number of observations.")
    }
    idx <- which(is.na(x) | is.na(y))
    if ( length(idx) > 0 ) {
        x <- x[-idx]
        y <- y[-idx]
    }
    ## We used to call the piecewise estimator the local estimator.
    ## For now, I'll fix it if someone's code calls it the local estimator.
    if ( estimator == "local" ) {
        estimator <- "piecewise"
        warning("Note the previously named 'local' estimator is now called ",
                "the 'piecewise' estimator. Using the 'local' terminology ",
                "is deprecated and will soon result in an error. ",
                "Correcting to 'piecewise' estimator.")
    }
    if ( estimator == "global" ) {
        x <- cbind(x, x > cutoff)
    } else if ( estimator == "piecewise" ) {
        x <- matrix(x)
    } else {
        stop("estimator should be one of global or piecewise.")
    }
    if ( is.null(b) ) {
        b <- rep(0, ncol(x) + 1)
    }
    if ( is.null(B) ) {
        B <- diag(10, nrow = ncol(x) + 1)
    }
    con <- list(fnscale = -1)
    con[names(control)] <- control
    hypers <- c(1, 1, rep(1, ncol(x)))
    X  <- cbind(1, x)
    ## Note slight difference in K0 to doc (no ell^2 division)
    if ( estimator == "piecewise" ) {
        res <- vector("list", 2)
        left_cases <- which(x[,1] < cutoff)
        xtmp <- x[left_cases, , drop = FALSE]
        Xtmp <- X[left_cases, , drop = FALSE]
        K0 <- .squared_distance(xtmp);
        Q0 <- Xtmp %*% B %*% t(Xtmp)
        M  <- Xtmp %*% b - y[left_cases]
        o  <- optim(hypers, .log_marginal_likelihood, .gradient, K0 = K0,
                    Q0 = Q0, M = M, method = method, control = con, ...)
        params <- o$par
        res[[1]] <- list(b = b, B = B, sigma_y = params[1], sigma_f = params[2],
                         ell = params[3:length(params)])
        right_cases <- which(x[,1] > cutoff)
        xtmp <- x[right_cases, , drop = FALSE]
        Xtmp <- X[right_cases, , drop = FALSE]
        K0 <- .squared_distance(xtmp);
        Q0 <- Xtmp %*% B %*% t(Xtmp)
        M  <- Xtmp %*% b - y[right_cases]
        o  <- optim(hypers, .log_marginal_likelihood, .gradient, K0 = K0,
                    Q0 = Q0, M = M, method = method, control = con, ...)
        params <- o$par
        res[[2]] <- list(b = b, B = B, sigma_y = params[1], sigma_f = params[2],
                         ell = params[3:length(params)])
    } else {
        K0 <- .squared_distance(x);
        Q0 <- X %*% B %*% t(X)
        M  <- X %*% b - y
        o  <- optim(hypers, .log_marginal_likelihood, .gradient, K0 = K0,
                    Q0 = Q0, M = M, method = method, control = con, ...)
        params <- o$par
        res <- list(b = b, B = B, sigma_y = params[1], sigma_f = params[2],
                    ell = params[3:length(params)])
    }
    return(res)
}
