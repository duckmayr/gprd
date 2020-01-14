#' Summary of gprd objects
#'
#' Print the treatment effect estimate and C.I.
#'
#' @param x A object of class \code{\link{gprd}}
#' @param digits An integer vector of length one giving the number of digits
#'     to print; default is 3
#'
#' @return Returns \code{NULL} invisibly
#'
#' @seealso gprd
#'
#' @export
summary.gprd <- function(x, digits = 3) {
    ## Sanity checks
    is_gprd <- inherits(x, "gprd")
    has_both_fs <- all(c("f_r", "f_l") %in% names(x))
    both_fs_are_gpr_objects <- all(sapply(x[c("f_r", "f_l")], inherits, "gpr"))
    if ( !is_gprd | !has_both_fs | !both_fs_are_gpr_objects ) {
        stop("x should be a gprd object returned from gprd().")
    }
    ## Prepare elements to print
    estimate <- sprintf(paste0("%0.", digits, "f"), x$tau_mean)
    ci <- sprintf(paste0("%0.", digits, "f"), x$tau_ci)
    ci_width <- as.numeric(names(x$tau_ci)[2]) - as.numeric(names(x$tau_ci)[1])
    ci_width <- round(100 * ci_width)
    header <- paste0("tau estimate        ", "[ ", ci_width, "% C.I. ]")
    n <- nchar(header)
    header <- paste0(strrep("-", n), "\n", header, "\n", strrep("-", n), "\n")
    ci <- paste0("[", ci[1], ", ", ci[2], "]")
    m <- n - nchar(ci) - nchar(estimate)
    body <- paste0(estimate, strrep(" ", m), ci, "\n")
    ## Print to console
    cat(header, body, sep = "")
    return(invisible(NULL))
}
