#' Construct Matrix from Formula and Names.
#' 
#' Function is a helper function for adding tables as X-axis labels
#' representing the structure of crossed (nested) factor-levels.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @usage getMatrix(form, labels, split="\\\\.")
#' @param form		(formula) object from which the matrix should be constructed
#' @param labels	(character) the 'names" element of the list returned by e.g.
#'                  \code{\link{BoxPlot}} or \code{\link{boxplot}}
#' @param split		(character) symbol or string specifying the split symbol or
#'                  string separating factor-levels in 'labels'
#' 
#' @examples
#' bp  <- boxplot(mpg~cyl:gear, mtcars, plot=FALSE)
#' mat <- getMatrix(mpg~cyl:gear, bp$names)
#' mat
#' 
#' @export

getMatrix <- function(form, labels, split="\\.")
{
    stopifnot(is.character(labels))
    stopifnot(class(form) == "formula")
    form <- terms(form)
    fac  <- rownames(attr(form, "factors"))
    fac  <- if(attr(form, "response")) fac[-1]
    mat  <- matrix(unlist(strsplit(labels, split)), ncol=length(labels))	
    if(nrow(mat) > length(fac))
        warning("Factor-level combinations include the split-character ", paste("'", split, "'", sep="")," causing this problem!")
    stopifnot(nrow(mat) == length(fac))	
    rownames(mat) <- fac
    return(mat)
}