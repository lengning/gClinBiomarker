#' Choose an R standard color or get the RGB-code of an arbitrary color.
#' 
#' Function lets the user choose a color via a tcltk widget in case it is called with default value of
#' paramater 'continuous'. Otherwise, it plots a grid of differently colored tiles, which can be used to select as many
#' colors as desired by simply clicking at the corresponding tiles.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @param continuous (logical) TRUE = one or multiple colors can be chosen from the RGB color-space via a tcltk-widget
#'                            FALSE = a grid is drawn where R standard colors can be selected by clicking on colored-tiles
#' 
#' @examples 
#' ## get R standard color
#'  get.colors(FALSE)
#' 
#' ## get RGB-code (hexadecimal) via tcltk-widget
#'  \dontrun{
#'      get.colors()
#'  }
#'  
#' @export

get.colors <- function(continuous=TRUE)
{
    if(continuous)
    {
        tt <- tktoplevel()
        tkwm.title(tt,"Color Selection")
        color <- "blue"
        n <- 1
        canvas <- tkcanvas(tt,width="80",height="25",bg=color)
        ChangeColor <- function()
        {
            color <<- tclvalue(.Tcl(paste("tk_chooseColor",.Tcl.args(initialcolor=color,title="Choose a color"))))
            if (nchar(color)>0)
                tkconfigure(canvas, bg=color)
        }
        AcceptColor <- function()
        {
            tkmessageBox(message=paste("You selected:", color), icon="info", type="ok")
            cat("\n", paste(n,".", paste(rep(" ", 2-nchar(n)), collapse=""), sep=""),"color:", color)
            n <<- n+1
        }
        Exit <- function()
        {
            tkdestroy(tt)
        }
        ChangeColor.button <- tkbutton(tt,text="Choose Color",command=ChangeColor)
        AcceptButton <- tkbutton(tt,text="Accept Color", command=AcceptColor)
        ExitButton <- tkbutton(tt, text="        Exit        ", command=Exit)
        tkgrid(tklabel(tt, text="                                   "))
        tkgrid(canvas,ChangeColor.button, tklabel(tt, text="      "), AcceptButton, tklabel(tt, text="      "))
        tkgrid(tklabel(tt, text="                                   "))   
        tkgrid(tklabel(tt, text=" "), tklabel(tt, text=" "), tklabel(tt, text=" "), ExitButton)   
        tkgrid(tklabel(tt, text="                                   "))   
        tkfocus(tt)   
    }
    else
    {
        plot(0:26,0:26, xlab="", ylab="", main="Choose Colors")
        Colors <- colors()
        rows <- numeric()
        cols <- numeric()
        mat <- matrix(nrow=26,ncol=26)
        for(i in 1:26)                              # rows
        {
            for(j in 1:26)                            # columns
            {
                rect(i-1, j-1, i, j, col=Colors[(i-1)*26+j])
                mat[i,j] <- Colors[(i-1)*26+j]
            }
        }
        while(TRUE)
        {
            loc <- locator(1)
            if(is.null(loc))
                break
            x <- ceiling(loc$x)
            y <- ceiling(loc$y)
            cat("\nx =",x,"y =",y,"color:",mat[x,y])
            flush.console()
        } 
    }    
}
