#' Enhanced Box-Percentile Plots.
#' 
#' Please use function \code{BoxPlot} with \code{box.type="bp"} for box-percentile plots, since this function here
#' was designed as helper function for function \code{BoxPlot}.
#' It uses the R-code of function \code{bpplot} which generates box-percentile plots as implemented in 
#' package \code{Hmisc}. This R-code is changed and augmented to implement an formula-interface and to be able to 
#' add box-percentile plots to an existing plot. It is assumed that \code{form} is a simple formula object, 
#' e.g. \code{y~grp:time} with just simple interactions or even simpler with just a single grouping factor, e.g. \code{y~grp}.
#' 
#' @author Andre Schuetzenmeister (using source-code of function \code{bpplot} from package \code{Hmisc})
#' 
#' @param obj (data.frame, matrix) corresponding to the dataset to be used
#' @param form (formula) such as \code{y~grp}, where \code{y} is a numeric vector of data values to be split into
#'              groups according to the grouping variable \code{grp}.
#' @param var (character) vector specifying the columns in 'obj' to be used for plotting. The order of elements is retained in the boxplot.
#' @param main (character) string giving the main title of the plot, left out if \code{add=TRUE}
#' @param add (logical) TRUE=add box-percentile plot(s) to an existing plot
#' @param col (character) vector specifying the colors of percentile-boxes, defaults to no color
#' @param labels (character) string(s) for group labels, which are drawn under each boxplot. Note: The order must correspond to the order of group-levels
#'               if the formula interface was used, i.e. check the order of sort(unique(obj$grp)).
#' @param ylab (character) character string specifying the optional label of the Y-axis (vertically centered)
#' @param ylab.line (numeric) specifying the line 'ylab' is put (useful for custom figure margins). 
#' @param ylab.font (integer) specifying the font to be used for the Y-axis label
#' @param swag (numeric) values within ]0,1] specifying the swaging of percentile boxes, for better separating boxes from each other
#' @param line.col (character) color of boundry lines of the percentile boxes
#' @param line.lwd (integer) line width of the boundry lines of the percentile boxes
#' @param line.lty (integer) line type of the boundry lines of the percentile boxes
#' @param add.xlab (logical) TRUE = automically determined group-labels are plotted below each box/percentile-box
#' @param ... additional graphical parameters passed on
#' 
#' @seealso \link{BoxPlot}, function \code{bpplot} in package \code{Hmisc}
#' 
#' @export

BPplot <- function (obj, form=NULL, var=NULL, main = "Box-Percentile Plot", add=FALSE, col=NULL, labels=NULL,
                    ylab=NULL, ylab.line=2.5, ylab.font=1L, swag=.9, line.col="black", line.lwd=1L,
                    line.lty=1L, add.xlab=TRUE, ...) 
{
    stopifnot(swag > 0 && swag <= 1)
    
    bpx <- function (y, offset, swag)                         # helper function for Hmisc-function 'bpplot'
    {
        y <- y[!is.na(y)]
        n <- length(y)
        delta <- 1/(n + 1)
        prob <- seq(delta, 1 - delta, delta)
        quan <- sort(y)
        med <- median(y)
        q1 <- median(y[y < med])
        q3 <- median(y[y > med])
        first.half.p <- prob[quan <= med]
        second.half.p <- 1 - prob[quan > med]
        plotx <- c(first.half.p, second.half.p)*swag
        qx <- approx(quan, plotx, xout = q1)$y
        q1.x <- c(-qx, qx) + offset
        qx <- approx(quan, plotx, xout = q3)$y
        q3.x <- c(-qx, qx) + offset
        q1.y <- c(q1, q1)
        q3.y <- c(q3, q3)
        med.x <- c(-max(first.half.p)*swag, max(first.half.p)*swag) + offset
        med.y <- c(med, med)
        return(list(x1 = (-plotx) + offset, y1 = quan, x2 = plotx + 
                        offset, y2 = quan, q1.y = q1.y, q1.x = q1.x, q3.y = q3.y, 
                    q3.x = q3.x, med.y = med.y, med.x = med.x))
    }
    
    if(!is.null(form) && class(form)=="formula")
    {
        char <- as.character(form)
        form <- formula(paste(paste(char[c(2,1,3)], collapse=""), "-1", sep="")) 
        mf <- model.frame(form, data=obj, na.action=na.pass)                    # NAs will be removed if not calling model.frame with "na.pass"
        mm <- model.matrix(form, mf)                                            
        mm <- apply(mm, 1:2, function(x) ifelse(x==0, NA, 1))                   # substitute NAs for 0s
        xnames <- boxplot(form, obj, plot=FALSE)$names 
        dat <- apply(mm, 2, function(x) obj[,char[2]] * x)                      # assign observations to classes
        all.x <- list()
        for(i in 1:ncol(dat))
            all.x[[i]] <- na.omit(dat[,i])
    }
    else
    {
        if(is.null(var))
        {
            warning("Boxplot cannot be drawn because there is neither formula 'form' nor variable names 'var' provided!")
            return(1)
        }
        xnames <- var
        all.x <- list()
        for(i in 1:length(var))
            all.x[[i]] <- na.omit(obj[,var[i]])
    }
    
    n <- length(all.x)
    centers <- seq(from = 1, by = 1, length = n)
    ymax <- max(sapply(all.x, function(x){
        if(all(is.na(x)))
            return(NA) 
        else 
            return(max(x, na.rm = TRUE))
    }), na.rm=TRUE)
    ymin <- min(sapply(all.x, function(x){
        if(all(is.na(x)))
            return(NA)
        else
            return(min(x, na.rm = TRUE))
    }), na.rm=TRUE)
    xmax <- max(centers) + 0.5
    xmin <- 0.5
    
    if(!add)
    {
        plot(c(xmin, xmax), c(ymin, ymax), type = "n", main = main, 
             xaxt = "n", xlab=NA, ylab=NA, ...)
        if(!is.null(ylab))
            mtext(ylab, side=2, line=ylab.line, font=ylab.font)
        if(is.null(labels) && add.xlab)
            axis(1, at=1:length(all.x), labels=xnames)
    }
    for (i in 1:n) {
        if(length(is.na(all.x[[i]])) < 2)           # function 'bxp' needs at least 2 non-NA values
            next
        plot.values <- bpx(all.x[[i]], centers[i], swag=swag)
        if(!is.null(col))
        {
            if(length(col) != length(all.x))
                col <- rep(col, ceiling(length(all.x)/length(col)))[1:length(all.x)] 
            if(length(line.col) != length(all.x))
                line.col <- rep(line.col, ceiling(length(all.x)/length(line.col)))[1:length(all.x)]
            
            polygon(c(plot.values$x1, rev(plot.values$x2)), 
                    c(plot.values$y1, rev(plot.values$y2)), col=col[i], border=line.col[i])
        }
        lines(plot.values$x1, plot.values$y1, col=line.col[i], lwd=line.lwd, lty=line.lty)
        lines(plot.values$x2, plot.values$y2, col=line.col[i], lwd=line.lwd, lty=line.lty)
        lines(plot.values$q1.x, plot.values$q1.y, col=line.col[i], lwd=line.lwd, lty=line.lty)
        lines(plot.values$q3.x, plot.values$q3.y, col=line.col[i], lwd=line.lwd, lty=line.lty)
        lines(plot.values$med.x, plot.values$med.y, col=line.col[i], lwd=line.lwd, lty=line.lty)
    }    
}