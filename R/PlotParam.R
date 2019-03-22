#' Set graphical parameters
#'
#' \code{PlotParam} can be used to set graphical parameters for a pdf file
#' or the default screen device.
#'
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}
#'
#' @param pdf.name A name of a pdf file.
#' @param pdf.param A list of parameters that define pdf graphics device. See \code{\link{pdf}}.
#' @param par.param A list of parameters that define graphcial parameters. See \code{\link{par}}.
#'
#' @note The \code{PlotParam} function is designed to restore graphical parameters \code{par}
#' to the default values after execution.
#'
#' @importFrom grDevices pdf dev.cur dev.off
#' @importFrom graphics par
#' @importFrom knitr all_labels
#'
#' @examples
#' pdf.param = list(height=5)
#' par.param = list(mar=c(4, 4, 3, 2))
#' pdf.name = NULL
#' PlotParam(pdf.name, pdf.param, par.param)
#' hist(rnorm(1000), col="blue")
#' PlotParam()
#'
#' @export

PlotParam <- function(pdf.name, pdf.param, par.param) {
    vec <- c(missing(pdf.name), missing(pdf.param), missing(par.param))
    if (all(vec)) {
        # Load old par on exit and shut down the graphical device.
        #if (names(dev.cur()) != "RStudioGD" & length(all_labels()) == 0) {

        if (names(dev.cur()) == "quartz_off_screen" & length(all_labels()) == 0) {
            par(old.par)
            invisible(dev.off())
        } else {
            # Load old par on exit. Does not shutdown the screen device.
            par(old.par)
        }
    } else {
        # Fix the issue of plotting multiple graphs in one graphical device.
        # It does not fix the issue in markdown though.
        if (names(dev.cur()) == "RStudioGD") {
            invisible(dev.off())
        }

        # Create/Modify pdf.param.
        if (!is.null(pdf.name)) {
            if (is.null(pdf.param)) {
                pdf.param <- list()
            }
            pdf.param$file <- pdf.name
            do.call(pdf, pdf.param)
        }
        # Save par to global variable.
        old.par <- NULL
        old.par <<- par(no.readonly=TRUE)

        # Create/Modify par.param.
        if (is.null(par.param)) {
            par.param <- list()
        }
        invisible(do.call(par, par.param))
    }
}
