#' Generate bar plot to summarize response outcome by subgroups
#'
#' This function generates bar plots to compare response outcome summarization (e.g. response rate) across different groups
#' (e.g. treatment arm, biomarker variable, demographics variable, full population vs bep, etc.)
#'
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param outcome.var column name of the response variable (e.g. best overall response).
#' entries with missing rsp value will be ignored in analysis.
#' @param binary whether summarize the response rate by responder/nonresponder (if binary=TRUE),
#' or by its original category (e.g. PD/SD/PR/CR). If binary is TRUE, responder categories
#' and nonresponder categories should be specified in rsp.response and rsp.nonresponse (all values in the outcome column
#' should be included in rsp.response and rsp.nonresponse),
#' at the same time rsp.levels will be ignored. If binary is FALSE, rsp.levels should be
#' specified to provide order of the categories. At the same time rsp.response and
#' rsp.nonresponse will be ignored.
#' @param rsp.response categories that should be considered as responder.
#' @param rsp.nonresponse categories that should be considered as non responder.
#' @param rsp.levels vector that indicates how to sort the response categories. This parameter will be ignored if binary is TRUE.
#' @param rsp.name Display name for the responders. Default is "rsp".
#' @param col color for different categories
#' @param plot.count default is FALSE. By default percentages will be shown.
#' If it is TRUE, will show counts instead
#' @param trt column name that indicates treatment variable. If this is not NULL,
#' the program will summarize the response outcome within each treatment arm.
#' @param trt.name treatment name to display
#' @param show.combine.trt if this is TRUE, will show summarization of combined data across treatment arms as well
#' @param compare.bep.itt whether want to generate two groups of bars to compare the summary statistics
#' in full population (e.g. ITT) vs in BEP. If this is TRUE, paramemeter bep should be specified. If this is FALSE,
#' parameters bep, bep.name, itt.name, bep.indicator will be ignored
#' @param compare.var whether want to generate multiple groups of bars to compare the summary statistics
#' in subgroups defined by var (e.g. categorical demographics or biomarker variable)
#' . If this is TRUE, paramemeter var should be specified. If this is FALSE,
#' parameters var, var.name, show.combine.var will be ignored
#' @param show.combine.var if this is TRUE, will show summarization of combined data across var levels as well
#' @param main main title of the barplot. Default is "Association of biomarker expression with response rate".
#' @param sub footnote under the bar plot. Default is NULL.
#' @param cex cex
#' @param digits see \code{\link{pdf}}. Number of digits to display for the response rate.
#' @param horiz whether show barplot horizontally
#' @param pdf.name name of output pdf file. If it's NULL, the plots will be displayed but not saved as pdf. Default is "rsp.barplot.pdf".
#' @param pdf.param list of parameters that define pdf graphics device. See \code{\link{pdf}}. Default is \code{list(width=6, height=4.5)}.
#' @param par.param list of parameters that define graphcial parameters. See \code{\link{par}}. Default is \code{list(mar=c(4,4,3,2))}.
#' @inheritParams SummaryVars
#' @inheritParams PlotTabForestBiomarker
#'
#' @return a bar plot that summarizes response outcome
#'
#' @examples
#' data(input)
#' PlotRspBar (input, outcome.var="Response",
#'            binary=TRUE,
#'            rsp.response=c("CR","PR"),
#'            rsp.nonresponse=c("SD", "PD","NON CR/PD","NE",NA),
#'            trt="Arm",
#'            compare.var=TRUE,
#'            var="KRAS.mutant")
#'
#' @export

PlotRspBar <- function(data, outcome.var,
                       binary=FALSE,
                       rsp.response = c("CR","PR"),
                       rsp.nonresponse = c("SD", "PD","NON CR/PD","NE"),
                       rsp.levels=c("CR", "PR","SD","NON CR/PD", "PD","NE"),
                       rsp.name = "Rsp",
                       col=NULL,
                       plot.count=FALSE,digits=1,
                       trt=NULL, trt.name =NULL, show.combine.trt=TRUE,
                       compare.bep.itt=FALSE, bep = NULL, bep.name = NULL, itt.name="Full",bep.indicator=1,
                       compare.var=FALSE, var=NULL, var.name=NULL, show.combine.var=TRUE, horiz=FALSE,
                       main="Association of response rate",
                       sub=NULL,cex=1,
                       pdf.name=NULL,
                       pdf.param=NULL,
                       par.param=NULL) {


  if(is.null(col)){
      if(binary==FALSE) col <- colorRampPalette(c("tomato","deepskyblue"))(length(rsp.levels))
      if(binary==TRUE) col <- c("tomato","gray")
  }
  stopifnot(class(data) == "data.frame")
  if(!all(c(var, trt, bep,outcome.var) %in% colnames(data)))stop("outcome.var, var, trt and bep should have matched column names in the input data!")


  if(length(which(is.na(data[,outcome.var])))>0){
    data <- data[!is.na(data[,outcome.var]),]
    message(paste("entries with missing outcome.var are removed!", nrow(data),"entries left"))
  }

  if(compare.bep.itt & is.null(bep)){
    compare.bep.itt <- FALSE
    message("compare.bep.itt=TRUE but bep is not specified. Reset test.bep as FALSE")
  }
  if(compare.var & is.null(var)){
    compare.var <- FALSE
    message("compare.var=TRUE but bep is not specified. Reset test.bep as FALSE")
  }

  if(binary)if(!all(unique(data[,outcome.var])%in%c(rsp.response,rsp.nonresponse)))
    stop("all unique values in outcome.var column should be included in rsp.response or rsp.nonresponse!")
  if(!binary)if(!all(unique(data[,outcome.var])%in%c(rsp.levels)))
    stop("all unique values in outcome.var column should be included in rsp.levels!")

  # generate response
  if(binary)
    data$rspvar <- factor(ifelse(data[,outcome.var]%in%rsp.response,rsp.name,paste0("non-",rsp.name)),levels=c(rsp.name,paste0("non-",rsp.name)))
  if(!binary) data$rspvar <- factor(data[,outcome.var], levels=rsp.levels)

  data$sectionvar <- data$trtvar <- rep("All", length(data$rspvar))

  if(compare.bep.itt){
    if(is.null(bep.name))bep.name <- bep
    data.bep <- subset(data, data[,bep]%in%bep.indicator)
    data.0 <- data
    data <- rbind(data,data.bep)
    data$sectionvar <- c(rep(itt.name,nrow(data.0)), rep(bep.name, nrow(data.bep)))
    data$sectionvar <- factor(data$sectionvar, levels=c(itt.name, bep.name))
  }

  if(compare.var){
    if(is.null(var.name))var.name <- var

    var.levs <- levels(factor(data[,var]))
    data$sectionvar <- factor(paste0(var.name,"(",data[,var],")"), levels=paste0(var.name,"(",var.levs,")"))
  }

  if(!is.null(trt)){
    data$trtvar <- factor(data[,trt])
  }


  section.levels <- levels(factor(data$sectionvar))
  trt.levels <- levels(factor(data$trtvar))

  tab.list <- sapply(section.levels,function(i){
    sapply(trt.levels,function(j){
      table(subset(data, sectionvar==i & trtvar==j)$rspvar)
    }, simplify=TRUE)
  },simplify=FALSE)


  tab.list0 <- tab.list
  # If more than one trt level, add an 'All' row
  if(length(trt.levels)>1 & show.combine.trt) for(i in section.levels){
    tab.list[[i]] <- cbind(tab.list[[i]], rowSums(tab.list[[i]]))
    colnames(tab.list[[i]])[ncol(tab.list[[i]])] <- "All"
  }

  tab.list1 <- tab.list
  # Add a 'All' row for var
  if(compare.var){
    if(length(section.levels)>1 & show.combine.var){
      tab.list[[length(tab.list)+1]] <- matrix(
        sapply(1:ncol(tab.list[[1]]),function(i)rowSums(sapply(tab.list,function(j)j[,i,drop=FALSE]))),
        nrow=nrow(tab.list[[1]]))
      rownames(tab.list[[length(tab.list)]]) <- rownames(tab.list[[1]])
      colnames(tab.list[[length(tab.list)]]) <- colnames(tab.list[[1]])
      names(tab.list) [length(tab.list)]<- "All"
    }
  }


  n.in.section <- ncol(tab.list[[1]])
  n.section <- length(tab.list)

  if(!(n.in.section==1&n.section==1))for(i in 1:n.section) colnames(tab.list[[i]]) <- paste0(names(tab.list)[i],":", colnames(tab.list[[i]]))

  tab.table <- do.call(cbind, tab.list)
  if(is.null(dim(tab.table)))tab.table <- matrix(tab.table, ncol=1, dimnames=list(rownames(tab.list[[1]]),"All")
  )

    plottab <- t(t(tab.table)/colSums(tab.table))
    space <- intersect((1:(n.section-1))*n.in.section +1, 1:ncol(plottab))
    myspace <- rep(0.2, ncol(plottab))
    myspace[space] <- 0.8
    if(n.in.section==1&n.section==1)myspace <- .2

    if(binary)bi.txt <- paste(paste("N=",colSums(tab.table,na.rm=TRUE),sep=""),
                    paste("; ",rsp.name,"=",round(plottab[rsp.name,],digits+2)*100,"%",sep=""), sep="")
    else bi.txt <- paste("N=",colSums(tab.table,na.rm=TRUE),sep="")
    max.bi.txt <- max(nchar(bi.txt))
    max.name.txt <- max(nchar(colnames(plottab)))

    # Parameters for the pdf file if it's created.
   if(!is.null(pdf.name)){
     if (is.null(pdf.param$height)) {
        pdf.param$height <- 8
    }

    if (is.null(pdf.param$width)) {
        pdf.param$width <- 8
    }
   }
    # Define par parameters
    if(is.null(par.param$mar)) {
        p1 <- max(5,.5*max.name.txt)
        p2 <- max(5.5,.55*max.bi.txt)
        if(horiz)par.param$mar <- c(5,p1,4,p2)
        if(!horiz) par.param$mar <- c(p1,4,p2,5)
    }

    if(is.null(par.param$font.axis)) {
        par.param$font.axis <- 3
    }

    PlotParam(pdf.name, pdf.param, par.param)

    num.axis <- ifelse(horiz,4,3)
    label.axis <- ifelse(horiz,3,4)
    perc.axis <- ifelse(horiz,1,2)

    if(!plot.count){
      xlab <- ifelse(horiz,"(%)","")
      ylab <- ifelse(horiz,"","(%)")
      bb <- barplot(plottab*100, horiz=horiz, col=col, space=myspace,
                  names=colnames(plottab),las=2, axes=FALSE, xlab=xlab,ylab=ylab,
                  cex.axis=0.7*cex)
    }
    if(plot.count){
      xlab <- ifelse(horiz,"count","")
      ylab <- ifelse(horiz,"","count")
      bb <- barplot(tab.table, horiz=horiz, col=col, space=myspace,
                                names=colnames(plottab),las=2, axes=FALSE,  xlab=xlab,
                    ylab=ylab,cex.axis=0.7*cex)
    }
    title(main,line=max(max.bi.txt*.45,4.5))

    axis(perc.axis, las=1)

    if (binary) {
        axis(num.axis, at=bb, bi.txt,
             las=2, cex.axis=0.8*cex, line=-0.5, tick=F)
    } else {
        axis(num.axis, at=bb, bi.txt, las=2, cex.axis=0.8*cex, line=-0.5, tick=F)
    }

    box()

    myoutcome.levels <- levels(data$rspvar)
    n.levels <- nlevels(data$rspvar)
    xmax0 <- 100
    if(plot.count) xmax0 <- max(colSums(tab.table))
    xmax <- xmax0
    if(binary)xmax <- xmax0/3
    nn <- nchar(myoutcome.levels)
    nn.unit <- xmax/sum(nn)
    tmp <- 0
    for (kk in 1:n.levels) {
      mtext(myoutcome.levels[kk], col=col[kk], side=label.axis, outer=FALSE, at=tmp, adj=0, las=3)
      tmp <- tmp +  nn.unit*nn[kk]
      }

   # abline(h=bb[space]/2+bb[space-1]/2)

    if (!is.null(sub)) {
        mtext(sub, side=1, outer=FALSE, cex=0.8*cex, line=4, adj=0)
    }

    PlotParam()

out <- list(count=t(tab.table), perc=t(plottab))
}
