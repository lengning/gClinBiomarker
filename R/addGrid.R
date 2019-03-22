#' Add a Grid to an Existing Plot.
#' 
#' It is possible to use automatically determined grid lines (\code{x=NULL, y=NULL}) or specifying the number 
#' of cells \code{x=3, y=4} as done by \code{grid}. Additionally, x- and y-locations of grid-lines can be specified,
#' e.g. \code{x=1:10, y=seq(0,10,2)}.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @param x (integer, numeric) single integer specifies number of cells, numeric vector specifies vertical grid-lines
#' @param y (integer, numeric) single integer specifies number of cells, numeric vector specifies horizontal grid-lines
#' @param col (character) color of grid-lines
#' @param lwd (integer) line width of grid-lines
#' @param lty (integer) line type of grid-lines
#' 
#' @export

addGrid <- function(x=NULL, y=NULL, col="lightgray", lwd=1L, lty=3L)
{
    if(all(is.null(c(x,y))) || all(length(c(x,y))<2))               # call grid function
        grid(nx=x, ny=y, col=col, lwd=lwd, lty=lty)
    else
    {
        if(length(x) == 0)                                          # NULL
            xticks <- axTicks(side=1)
        else if(length(x) == 1)
        {
            U <- par("usr")
            xticks <- seq.int(U[1L], U[2L], length.out = x + 1)
        }
        else
            xticks <- x
        
        if(length(y) == 0)                                          # NULL
            yticks <- axTicks(side=2)
        else if(length(y) == 1)
        {
            U <- par("usr")
            yticks <- seq.int(U[3L], U[4L], length.out = y + 1)
        }
        else
            yticks <- y
        
        abline(v=xticks, col=col, lwd=lwd, lty=lty)
        abline(h=yticks, col=col, lwd=lwd, lty=lty)
    }
}