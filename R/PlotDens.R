#' Plot Density
#' 
#' Generate density plot for a continuous covariate.
#' 
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}, and previous team members (see DESCRIPTION)
#' 
#' @param data Input data frame. Rows are patients and columns are variables (e.g. demographics variables, time to event variables, 
#' biomarker variables, treatment indicator, etc.). One patient per row. 
#' @param var Name of a variable (e.g. a biomarker variable). Should be in colnames of \code{data}.  
#' @param log2 If TRUE, computes binary (i.e. base 2) logarithm. Default is FALSE.
#' @param col The color of the line segments. Default is blue.
#' @param text.font Legend text font size. Default is 3.
#' @param main The main title. Default is \code{paste("Distribution of", var)}. A variable name will be added at the end of the title.
#' @param add.num The constant to add to all values. Helps to avoid applying log transformation on 0 or negative values. Default is 0.
#' @param xlab X axis label. Default is \code{paste(var, ifelse(log2 == TRUE, "(log2 scale)", ""))}. A variable name will be added at the beginning of the label.
#' @param pdf.name Name of output pdf file. If it's NULL (default), the plots will be displayed but not saved as pdf.
#' @param pdf.param A list of parameters that define pdf graphics device. See \code{\link{pdf}}. Default is \code{list(width=6, height=4.5)}. 
#' @param par.param A list of parameters that define graphcial parameters. See \code{\link{par}}. Default is \code{list(mar=c(4,4,3,2))}.
#' 
#' @return Histogram and density will be shown. Summary statistics will be also shown on the plots.
#' 
#' @examples
#' data(input)
#' PlotDens(data=input, var="KRAS.exprs", log2=TRUE)
#' 
#' @export

PlotDens <- function(data,
                     var,
                     log2=FALSE,
                     col="blue",
                     text.font=3,
                     main=paste("Distribution of", var),
                     add.num=0,
                     xlab=paste(var, ifelse(log2 == TRUE, "(log2 scale)", "")),
                     pdf.name=NULL,
                     pdf.param=list(width=6, height=4.5),
                     par.param=list(mar=c(4,4,3,2))) {
    
    if (!(is.data.frame(data))) {
        stop("An input is not a data frame!")
    }
    
    if (missing(var) | length(var) == 0) {
        stop("Please specify at least one covariate!")
    }
    
    sel <- which(!var %in% colnames(data))
    if (length(sel) > 0) {
        stop(paste("A covariate", var[sel], "is not in the data frame! "))
    }  
    
    PlotParam(pdf.name, pdf.param, par.param)  
    
    V <- data[, var] + add.num
    
    if (log2 == TRUE) {
        if (any(data[, var] <= 0, na.rm=T)) {
            stop(paste(var, " contains values less than or equal 0! 
                       No log transformation possible! You may set add.num to add a constant to all values."))
        }
        V <- log2(V)
    }
    
    sd.v <- round(sd(V, na.rm=TRUE), 2)
    s.v <- summary(V)
    leg.text <- paste(c("Mean", "SD", "Median", "Range"),
                      c(round(s.v, 2)["Mean"], 
                        sd.v, 
                        round(s.v, 2)["Median"],
                        paste("(", round(s.v["Min."], 2), " - ", round(s.v["Max."], 2), ")", sep="")),
                      sep=": ")
    
    y2 <- max(density(V, na.rm=T)$y, na.rm=TRUE)
    hist(V, prob=T, main=main, xlab=xlab, col="grey", ylim=c(0, y2))
    lines(density(V, na.rm=T), col=col, lwd=2)
    legend("topright", leg.text, bty="n", text.font=text.font)
    mtext(side=3, line=0, paste("N=", sum(!is.na(V))))
    box()
    
    PlotParam()
}