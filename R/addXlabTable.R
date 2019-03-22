#' Create Tabular Environment from Matrix and Add it to Bottom Margin.
#' 
#' This is a helper function intended to add tabular environments as X-axis
#' labels representing specific structure, e.g. formulas in boxplots.
#' Function takes a matrix and constructs a tabular environment from this.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @param mat			(matrix) representing the table which will be added to lower margin;
#'                      rows represent factors, columns represent factor-levels, rownames will
#'                      be re-used as rownames of the table in the plot
#' @param merge			(logial) TRUE = neighboring cells of 'mat' with the same content will
#'                      be merged to a wider cell
#' @param bmargin		(numeric) value representing the height of the (b)ottom margin as proportion
#'                      of the total height of the figure region. Note, that this will determine the height
#'                      of the table in the bottom margin also depending on the amount of space available there.
#' @param Label			(list) specifying all parameters applicable in function 'text', x- and y-values
#'                      will be set automatically for factor-labels (rows of the table)
#' @param Text			(list) specifying all parameters applicable in function 'text', x- and y-values will
#'                      be set automatically for character strings appearing in the cells of the table
#' 
#' @examples 
#' bp  <- BoxPlot(mpg~cyl:gear, mtcars, Xaxis=NULL)
#' mat <- getMatrix(mpg~cyl:gear, bp$names)
#' addXlabTable(mat)
#' 
#' # changing appearance for labels and cell-content
#' 
#' bp  <- BoxPlot(mpg~cyl:gear, mtcars, Xaxis=NULL)
#' addXlabTable(mat, Label=list(font=2, col="red"), Text=list(font=3, col="blue"))
#' 
#' @export

addXlabTable <- function(mat, merge=TRUE, bmargin=.025, Label=list(), Text=list())
{
    stopifnot(is.logical(merge))
    stopifnot(is.matrix(mat))
    stopifnot(is.numeric(bmargin))
    stopifnot(0 <= bmargin && bmargin <= 1)
    
    old.par <- par(xpd=TRUE)								# stop clipping to plotting region
    
    usr <- par("usr")										# user-coordinates of the plotting region ("plt")
    fig <- par("fig")										# figure-coordinates
    plt <- par("plt")										# plotting-region coordinates on the scale of "usr"
    
    dxp <- plt[2] - plt[1]									# delta_x in plt-coordinates
    dyp <- plt[4] - plt[3]									# delta_y in plt-coordinates
    dxu <- usr[2] - usr[1]									# delta_x in usr-coordinates
    dyu <- usr[4] - usr[3]									# delta_y in usr-coordinates
    
    fuc    <- numeric(4)									# figure-region in user coordinates
    fuc[1] <- usr[2] - dxu * plt[2] / dxp  
    fuc[2] <- fuc[1] + (usr[2] - fuc[1]) / plt[2]
    fuc[3] <- usr[4] - dyu * plt[4] / dyp  
    fuc[4] <- fuc[3] + (usr[4] - fuc[3]) / plt[4]
    
    l  <- usr[1]											# left border
    r  <- usr[2]											# right border
    u  <- usr[3]											# upper border
    b  <- fuc[3] + bmargin * diff(fuc[3:4])					# bottom border
    w  <- abs(diff(c(l,r)))									# table width
    h  <- abs(diff(c(u, b)))								# table height
    rh <- h/nrow(mat)										# row height
    cw <- w/ncol(mat)										# cell width
    
    tmp.Text <- list()
    tmp.Text[names(Text)] <- Text
    Text <- tmp.Text
    
    # draw table
    
    tmp.l <- NA
    
    for( i in 1:(nrow(mat)+1) )
    {
        lines(c(l,r), rep(u-(i-1)*rh, 2))							# horizontal lines
        
        for( j in 1:(ncol(mat)+1) )									# vertical lines
        {
            if(i == nrow(mat) + 1)
                next
            
            if(j %in%  c(1, ncol(mat)+1) )							# outer lines
            {
                lines(rep(l+(j-1)*cw, 2), c(u-(i-1)*rh, u-i*rh))
            }
            else													# inner lines
            {
                if(mat[i,j-1] != mat[i,j] || !merge)				# draws lines between cells with different content
                {
                    swit  <- TRUE
                    lines(rep(l+(j-1)*cw, 2), c(u-(i-1)*rh, u-i*rh))
                }
            }
            
            if(j < ncol(mat) + 1)									# add text to cells
            {
                if( j == 1 )
                    tmp.l <- l
                
                if( merge )
                {
                    if( j < ncol(mat) )
                    {
                        if(mat[i, j] != mat[i, j+1])
                        {
                            tmp.r <- l + j * cw	
                            Text$x <- mean(c(tmp.l, tmp.r))
                            Text$y <- mean(c(u-(i-1)*rh, u-i*rh))
                            Text$labels <- mat[i,j]
                            Text$adj <- c(.5, .5)
                            do.call("text", Text)
                            tmp.l <- tmp.r
                        }						
                    }
                    else
                    {
                        tmp.r <- l + j * cw	
                        
                        Text$x <- mean(c(tmp.l, tmp.r))
                        Text$y <- mean(c(u-(i-1)*rh, u-i*rh))
                        Text$labels <- mat[i,j]
                        Text$adj <- c(.5, .5)
                        do.call("text", Text)
                    }
                }
                else
                {
                    Text$x <- l+(j-1)*cw+cw/2
                    Text$y <- mean(c(u-(i-1)*rh, u-i*rh))
                    Text$labels <- mat[i,j]
                    Text$adj <- c(.5, .5)
                    
                    do.call("text", Text)
                }				
            }
            
        }
    }
    
    # add rownames
    tmp.Label <- list(x = l - .01 * (fuc[2]-fuc[1]),
                      y = seq(u-rh/2, b+rh/2, by=ifelse(u-rh/2 > b+rh/2, -rh,rh)),
                      labels=rownames(mat),
                      adj=c(1, .5))
    
    tmp.Label[names(Label)] <- Label
    Label <- tmp.Label	
    do.call("text", Label)
    
    par(old.par)
}
