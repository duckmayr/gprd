#' Summary of gprd objects
#'
#' Print the treatment effect estimate and C.I.
#'
#' @param object A object of class \code{\link{gprd}}
#' @param digits An integer vector of length one giving the number of digits
#'     to print; default is 3
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns \code{NULL} invisibly
#'
#' @seealso gprd
#'
#' @export
summary.gprd <- function(object, digits = 3, ...) {
    ## Sanity checks
    is_gprd <- inherits(object, "gprd")
    all_fs_are_gpr_objects <- all(sapply(object[["f"]], inherits, "gpr"))
    if ( !is_gprd | !all_fs_are_gpr_objects ) {
        stop("object should be a gprd object returned from gprd().")
    }
    ## Prepare elements to print
    estimate <- sprintf(paste0("%0.", digits, "f"), object$tau_mean)
    ci <- sprintf(paste0("%0.", digits, "f"), object$tau_ci)
    ci_width <- round(100 * diff(as.numeric(rev(names(object$tau_ci)))))
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
