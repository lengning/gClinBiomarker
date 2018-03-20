#' Create Tabular Environment from Matrix and Add it to Margin.
#' 
#' This is a generalization of function \code{\link{addXlabTable}}, intended to be a helper 
#' function adding tables to plots instead of usual axis annotations.
#' These tables could represent specific structures, e.g. formulas in boxplots.
#' This function takes a matrix and constructs a tabular environment from this. For
#' all options of argument 'margin' the table will always be ordered from inside to outside,
#' e.g. in boxplots factor levels of the nested factor will appear near to the plot and the
#' outer factor will appear more distant (see examples).
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @param mat			(matrix) representing the table which will be added to margin;
#'                      rows represent factors, columns represent factor-levels, rownames will
#'                      be re-used as rownames of the table in the plot
#' @param margin		(character) string specifying the margin to which the table should be added,
#'                      partial matching is supported
#' @param merge			(logial) TRUE = neighboring cells of 'mat' with the same content will
#'                      be merged to a wider cell. If provided as vector of logicals, each elements will be
#'                      applied to rows of 'mat' separately, allowing to merge some variables and leaving others untouched.
#' @param propMargin	(numeric) value representing the height (bottom, top) or width (left, right) of the margin as proportion
#'                      of the total height of the figure region. Note, that this will determine the height
#'                      of the table in the bottom margin also depending on the amount of space available there.
#' @param Label			(list) specifying all parameters applicable in function 'text', x- and y-values
#'                      will be set automatically for factor-labels (rows of the table)
#' @param Text			(list) specifying all parameters applicable in function 'text', x- and y-values will
#'                      be set automatically for character strings appearing in the cells of the table
#' @param reorder		(logical) TRUE = 'mat' will be reorder to match the way horizontal boxplots are drawn,
#'                      i.e. the column order will be inverted, FALSE = 'mat' will be used in the given order
#' 
#' @seealso \code{\link{addXlabTable}}
#' 
#' @examples 
#' # these examples should work when run by the user
#' \dontrun{
#' old.par <- par(mar=c(8,4,4,2), yaxs="r")
#' bp  <- boxplot(mpg~cyl:gear, mtcars)
#' bxp(bp, axes=FALSE)
#' box()
#' axis(2)
#' mat <- getMatrix(mpg~cyl:gear, bp$names)
#'
#' # default is table below the plot (margin="bottom") 
#' addTableToMargin(mat)
#' par(old.par)
#'
#' # place table in left margin
#' old.par <- par(mar=c(4,8,4,3), xaxs="r")
#' bxp(boxplot(mpg~cyl:gear, mtcars, horizontal=TRUE ), axes=FALSE, 
#' 	   horizontal=TRUE)
#' box()
#' axis(1)
#' addTableToMargin(mat, margin="left", reorder=TRUE)
#' par(old.par)
#' 
#' # place table in upper margin
#' old.par <- par(mar=c(4,3,8,2), yaxs="r")
#' bxp(boxplot(mpg~cyl:gear, mtcars), axes=FALSE)
#' box()
#' axis(2)
#' addTableToMargin(mat, margin="top") 
#' par(old.par)
#'
#' # place table in right margin
#' old.par <- par(mar=c(4,3,3,8))
#' bxp(boxplot(mpg~cyl:gear, mtcars, horizontal=TRUE), axes=FALSE, 
#' 	   horizontal=TRUE)
#' box()
#' axis(1)
#' addTableToMargin(mat, margin="right", reorder=TRUE) 
#' par(old.par)
#' 
#' # changing appearance for labels and cell-content
#' # Note: table is added as is, without adapting to the way the boxplot is drawn
#' old.par <- par(mar=c(4,3,3,8))
#' bxp(boxplot(mpg~cyl:gear, mtcars, horizontal=TRUE), axes=FALSE, horizontal=TRUE)
#' box()
#' axis(1)
#' addTableToMargin(mat, margin="right", Label=list(font=2, col="red"), 
#'					Text=list(font=3, col="blue")) 
#' par(old.par)
#' }
#' 
#' @export

addTableToMargin <- function(mat, margin=c("bottom", "left", "top", "right"),
                             merge=TRUE, propMargin=.025, Label=list(), Text=list(),
                             reorder=FALSE)
{
    stopifnot(is.logical(merge))
    stopifnot(is.matrix(mat))
    stopifnot(is.numeric(propMargin))
    stopifnot(0 <= propMargin && propMargin <= 1)
    
    if(margin %in% c("left", "right") && reorder)
    {
        tmp.mat <- matrix(ncol=ncol(mat), nrow=nrow(mat))	#reorder matrix to be consistent with horizontal
        for(i in ncol(mat):1)
            tmp.mat[,(ncol(mat)+1)-i] <- mat[,i]
        rownames(tmp.mat) <- rownames(mat)
        mat <- tmp.mat
    }
    
    if(length(merge) < nrow(mat))
        merge <- rep(merge, ceiling(nrow(mat) / length(merge)))		# replicate if necessary
    
    margin <- match.arg(margin)
    
    old.par <- par(xpd=TRUE)								# stop clipping to plotting region
    
    usr <- par("usr")										# user-coordinates of the plotting region ("plt")
    fig <- par("fig")										# figure-coordinates
    plt <- par("plt")										# plotting-region coordinates on the scale of "usr"
    
    dxp <- plt[2] - plt[1]									# delta_x in plt-coordinates
    dyp <- plt[4] - plt[3]									# delta_y in plt-coordinates
    dxu <- usr[2] - usr[1]									# delta_x in usr-coordinates
    dyu <- usr[4] - usr[3]									# delta_y in usr-coordinates
    
    fuc    <- numeric(4)									# figure-region in user coordinates
    fuc[1] <- usr[2] - dxu * plt[2] / dxp  					# x1, x2, y1, y2
    fuc[2] <- fuc[1] + (usr[2] - fuc[1]) / plt[2]
    fuc[3] <- usr[4] - dyu * plt[4] / dyp  
    fuc[4] <- fuc[3] + (usr[4] - fuc[3]) / plt[4]
    
    if(margin %in% c("top", "bottom"))
    {
        l  <- usr[1] 										# left border
        r  <- usr[2]										# right border
        
        u  <- if(margin == "top")
            fuc[4] - propMargin * diff(fuc[3:4])		# upper border above upper figure region y
        else
            usr[3]									# bottom -> lower figure region y
        
        b  <- if(margin == "top")
            usr[4] 									# bottom border
        else
            fuc[3] + propMargin * diff(fuc[3:4])		
    }
    else
    {
        l  <- if(margin == "left")
            fuc[1] + propMargin * diff(fuc[1:2])		# left border near display region
        else
            usr[2]									# right -> right figure region
        
        r  <- if(margin == "left")
            usr[1] 									# right border
        else
            fuc[2] - propMargin * diff(fuc[1:2])									# left -> left figure region
        
        u  <- usr[4]										# upper border
        b  <- usr[3] 										# bottom border
    }
    
    w  <- abs(diff(c(l, r)))								# table width
    h  <- abs(diff(c(u, b)))								# table height
    
    if(margin %in% c("top", "bottom"))
    {
        rh <- h/nrow(mat)									# row height
        cw <- w/ncol(mat)									# cell width
    }
    else
    {
        cw <- w/nrow(mat)									# row height
        rh <- h/ncol(mat)									# cell width
    }
    
    tmp.Text <- list()
    tmp.Text[names(Text)] <- Text
    Text <- tmp.Text
    
    if(length(Text) > 0)
    {
        for(i in 1:length(Text))							# replicate arguments to a length equal to ncol(mat)
        {
            if(length(Text[[i]]) < ncol(mat))
                Text[[i]] <- rep(Text[[i]], ceiling(ncol(mat) / length(Text[[i]])))
        }
    }
    
    TextSpec <- Text										# keep current state of Text as (Spec)ified
    
    # draw table
    
    tmp.l <- NA
    
    for( i in 1:(nrow(mat)+1) )										# over factor-variables
    {
        if(margin %in% c("top", "bottom"))
            lines(c(l,r), rep(u-(i-1)*rh, 2))						# horizontal lines (top, bottom)
        else
            lines(rep(l+(i-1)*cw, 2), c(b,u))						# vertical line (left, right)																	
        
        for( j in 1:(ncol(mat)+1) )									# over factor-levels, vertical (t,b), horizontal (l,r) lines
        {
            if( i == (nrow(mat) + 1) )
                next
            
            if(length(TextSpec) > 0 && j < ncol(mat) + 1)			# anything specified by the user?
            {
                tmpText <- vector("list", length(TextSpec))			# extract j-th element of argument-vectors in 'TextSpec'
                names(tmpText) <- names(TextSpec)
                
                for(k in 1:length(TextSpec))
                    tmpText[[k]] <- TextSpec[[k]][j]
                
                Text <- tmpText
            }
            
            if(j %in%  c(1, ncol(mat)+1) )							# outer lines
            {
                if(margin %in% c("top", "bottom"))
                    lines(rep(l+(j-1)*cw, 2), c(u-(i-1)*rh, u-i*rh))
                else
                    lines(c(r-(i-1)*cw, r-i*cw), rep(b+(j-1)*rh, 2))
            }
            else																		# inner lines
            {
                if(any(mat[i:nrow(mat),j-1] != mat[i:nrow(mat),j]) || !merge[i])		# draws lines between cells with different content
                {
                    if(margin == "top")
                        lines(rep(l+(j-1)*cw, 2), c(b+(i-1)*rh, b+i*rh))
                    else if(margin == "bottom")
                        lines(rep(l+(j-1)*cw, 2), c(u-(i-1)*rh, u-i*rh))
                    else if(margin == "left")
                        lines(c(r-(i-1)*cw, r-i*cw), rep(b+(j-1)*rh, 2))
                    else
                        lines(c(l+(i-1)*cw, l+i*cw), rep(b+(j-1)*rh, 2))
                }
            }
            
            if(j < ncol(mat) + 1)									# add text to cells
            {
                if( j == 1 )
                {
                    if(margin %in% c("top", "bottom"))
                        tmp.l <- l
                    else
                        tmp.u <- u
                }					
                
                if( merge[i] )
                {
                    if( j < ncol(mat) )
                    {
                        if(any(mat[i:nrow(mat), j] != mat[i:nrow(mat), j+1]))		# cell different from neighboring cell
                        {
                            if(margin %in% c("top", "bottom"))
                            {
                                tmp.r <- l + j * cw	
                                Text$x <- mean(c(tmp.l, tmp.r))
                                if(margin == "bottom")
                                    Text$y <- mean(c(u-(i-1)*rh, u-i*rh))
                                else
                                    Text$y <- mean(c(b+(i-1)*rh, b+i*rh))
                                Text$labels <- mat[i,j]
                                Text$adj <- c(.5, .5)
                            }
                            else
                            {
                                tmp.b <- u - j * rh
                                if(margin == "left")
                                    Text$x <- mean(c(r-(i-1)*cw, r-i*cw))
                                else
                                    Text$x <- mean(c(l+(i-1)*cw, l+i*cw))
                                Text$y <- mean(c(tmp.u, tmp.b))
                                Text$labels <- mat[i,j]
                                Text$adj <- c(.5, .5)
                            }						
                            
                            do.call("text", Text)
                            
                            if(margin %in% c("top", "bottom"))
                                tmp.l <- tmp.r
                            else
                                tmp.u <- tmp.b
                        }						
                    }
                    else
                    {
                        if(margin %in% c("top", "bottom"))
                        {
                            tmp.r <- l + j * cw	
                            
                            Text$x <- mean(c(tmp.l, tmp.r))
                            if(margin == "bottom")
                                Text$y <- mean(c(u-(i-1)*rh, u-i*rh))
                            else
                                Text$y <- mean(c(b+(i-1)*rh, b+i*rh))
                            Text$labels <- mat[i,j]
                            Text$adj <- c(.5, .5)
                        }
                        else
                        {
                            tmp.b <- u - j * rh
                            
                            if(margin == "left")
                                Text$x <- mean(c(r-(i-1)*cw, r-i*cw))
                            else
                                Text$x <- mean(c(l+(i-1)*cw, l+i*cw))
                            Text$y <- mean(c(tmp.u, tmp.b))
                            Text$labels <- mat[i,j]
                            Text$adj <- c(.5, .5)
                        }
                        
                        do.call("text", Text)
                    }
                }
                else
                {
                    if(margin %in% c("top", "bottom"))
                    {
                        Text$x <- l+(j-1)*cw+cw/2
                        if(margin == "bottom")
                            Text$y <- mean(c(u-(i-1)*rh, u-i*rh))
                        else
                            Text$y <- mean(c(b+(i-1)*rh, b+i*rh))
                        Text$labels <- mat[i,j]
                        Text$adj <- c(.5, .5)
                    }
                    else
                    {
                        if(margin == "left")
                            Text$x <- r-(i-1)*cw-cw/2
                        else
                            Text$x <- l+(i-1)*cw+cw/2
                        Text$y <- mean(c(u-(j-1)*rh, u-j*rh))
                        Text$labels <- mat[i,j]
                        Text$adj <- c(.5, .5)
                    }
                    
                    do.call("text", Text)
                }				
            }
            
        }
    }
    
    # add rownames
    
    if(margin %in% c("top", "bottom"))
    {
        if(margin == "bottom")
        {
            tmp.Label <- list(	x = l - .01 * (fuc[2]-fuc[1]),
                               y = seq(u-rh/2, b+rh/2, by=-rh),
                               labels=rownames(mat),
                               adj=c(1, .5))
        }
        else
        {
            tmp.Label <- list(	x = l - .01 * (fuc[2]-fuc[1]),
                               y = seq(b+rh/2, u-rh/2, by=rh),
                               labels=rownames(mat),
                               adj=c(1, .5))
        }
        
        tmp.Label[names(Label)] <- Label
        Label <- tmp.Label	
    }
    else
    {
        if(margin == "left")
        {
            tmp.Label <- list(	x = seq(r-cw/2, l+cw/2, by=-cw),								
                               y = u + .01 * (fuc[4]-fuc[3]),
                               labels=rownames(mat),
                               adj=c(.5, 0))
        }
        else
        {
            tmp.Label <- list(	x = seq(l+cw/2, r-cw/2, by=cw),								
                               y = u + .01 * (fuc[4]-fuc[3]),
                               labels=rownames(mat),
                               adj=c(.5, 0))
        }
        tmp.Label[names(Label)] <- Label
        Label <- tmp.Label
    }
    
    do.call("text", Label)
    
    par(old.par)
}