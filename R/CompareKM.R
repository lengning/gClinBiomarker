#' Compare KM curve of full population vs. BEP
#'
#' This function provides K-M curves to compare full population vs. BEP.
#'
#' @author Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param tte column name that indicates the time to event variable
#' @param cens column name that indicates the censoring variable associated to tte. 1 indicates event and 0 indicates censoring
#' @param col.itt,col.bep,col.ci Color for itt curve, bep curve and confidence interval curve (CI)
#' @param surv.conf.type confidence interval type. Default is "plain". see conf.type in survfit
#' @param shaded.ci Whether add background shade to disply CI
#' @param pdf.name name of output pdf file. If it is NULL (default), the plots will be displayed but not saved as pdf
#' @param pdf.param A list of parameters that define pdf graphics device. See \code{\link{pdf}}.
#' for example, pdf.param=list(width=10, height=5)
#' @param par.param A list of parameters that define graphcial parameters. See \code{\link{par}}.
#' for example, par.param=list(mfrow=c(3,1), mar=c(4,4,3,2))
#' @param xlim,xat,ylab,xlab,main see \code{\link{plot}}
#' @param ... additional parameters for \code{\link{plot}}
#'
#' @note This function generates KM curves to compare full population vs BEP, within each treatment arm.
#'
#' @importFrom graphics plot axis mtext grid box polygon lines legend
#' @importFrom survival survfit
#'
#' @inheritParams SummarySingle 
#' @inheritParams PlotTabForestBiomarker
#' 
#'
#' @examples
#' data(input)
#' sample.data <- input
#' CompareKM(data=sample.data, tte="OS",cens="OS.event", main="OS ITT", bep="BEP")
#' @export
#' @export

CompareKM <- function(data, tte, cens, trt=NULL, bep, bep.indicator=1,
		      bep.name="Biomarker Evaluable", itt.name="All",
		      col.itt="palegreen4", col.bep="lightpink2", col.ci="lightcyan",shaded.ci=TRUE,
		      xlim=NULL, xat=NULL, ylab=paste(tte,"Survival Probability"), xlab="Time", main="",
		      surv.conf.type="plain",
                      pdf.name = NULL, pdf.param=list(height=5), par.param=list(mar=c(4,4,3,2)),...){

  if(is.null(xlim)){ xlim <- c(0, max(data[,tte], na.rm=TRUE)) + c(0, 0.1*max(data[,tte], na.rm=TRUE)) }

  lev <- NULL
  nlev <- 1
  if(!is.null(trt)){
  lev <- levels(factor(data[,trt]))
  nlev <- length(lev)
  }


  if(is.null(par.param$mfrow))par.param$mfrow <- c(1,nlev) # 1 by nlev subplots
  if(is.null(pdf.param$height))pdf.param$height <- 5
  if(is.null(pdf.param$width))pdf.param$width <- nlev * pdf.param$height # width = nlev * height to display the subplots

  PlotParam(pdf.name, pdf.param, par.param)
  for(i in 1:nlev){

    if(!is.null(trt))tmp <- data[which(data[,trt] == lev[i]),]
    if(is.null(trt))tmp <- data

    sf <- survfit(as.formula(paste("Surv(",tte,",",cens,")~1")), data=tmp, conf.type=surv.conf.type)
    plot(sf, col=col.itt, ylab=ylab, xlab=xlab, xlim=xlim, axes=FALSE, mark.time=FALSE,...)
    axis(2)
    if(!is.null(xat)){axis(1, at=xat, labels=xat)}
    if(is.null(xat)){axis(1)}
    mtext(side=3, paste(main, "\n", lev[i]), col="black", line=1.2, cex=1.5, font=2)
    grid()
    box()

    if(shaded.ci == TRUE){
      # if there are NAs in the CI (may happen at the end of the curve -- remove)
      idx <- which(!is.na(sf[["upper"]]) & !is.na(sf[["lower"]]))

      # shaded CI region
      x <- c(0,sf[["time"]][idx]) #add at time 0
      xp <- x
      up <- sf[["upper"]][idx]
      lo <- sf[["lower"]][idx]

      # however, add a timepoint at the end with missing CI to continue until then (first missing)
      # this might not work in all situations!!
      idx2 <- which(is.na(sf[["upper"]]) & is.na(sf[["lower"]]))
      if(length(idx2) > 0){
        idx2 <- idx2[1]
        xp[length(xp)+1] <- sf[["time"]][idx2]
        up[length(up)+1] <- up[length(up)]
        lo[length(lo)+1] <- lo[length(lo)]
      }

      yp <- c(1,up,rev(lo),1)

      #since it is a step function need to repeat some points to get the steps in the CIs (!!!)
      xp <- c(xp[1],rep(xp[2:length(xp)],each=2), xp[length(xp)])
      xp <- c(xp,rev(xp))
      yp <- rep(yp, each=2)
      polygon(xp,yp,density=NULL,col=col.ci, lwd=0.5, border=NA)
    }


    lines(sf, mark.time=FALSE, col=col.itt, lwd=3)

    sfflag <- survfit(as.formula(paste("Surv(",tte,",",cens,")~1")), data=tmp[which(tmp[,bep]==bep.indicator),], conf.type="log")
    lines(sfflag, mark.time=FALSE, col=col.bep, lwd=3)

    legend("topright", lty=1, lwd=3, col=c(col.itt, col.bep), legend=c(itt.name,bep.name), bg="white")

  }
  PlotParam()
}
