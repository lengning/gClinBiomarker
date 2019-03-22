#' Generate KM curve(s) for full population or subgroups defined by single factor or multiple factors
#'
#' The function generates KM curves for full population or subgroups. The subgroups may be defined as a single factor or multiple factors.
#'
#' @author Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param bep name of the column which indicates subpopulation (e.g. biomarker evaluable population)
#' If parameter bep is not defined, the KM curve(s) will be draw using all samples.
#' If bep is defined, the KM curve(s) will be draw using only samples in BEP.
#' @param varlist name (or names) of the column which indicates the subgroups (e.g. treatment group). It is supposed to be a vector.
#' This is an alternative option for specifying variable of interest (instead of specifying trt and var).
#' trt and var will be ignored if varlist is not NULL.
#' Compare to specifying trt and var, user can input any numbers of variables to varlist (a vector of column names).
#' Any specified column is expected to be categorical. If one column is in character class and var.levels is not specified,
#' it will be converted to a factor by factor() function. If varlist.levels is defined, the column will be converted to
#' a factor following the level order in varlist.levels.
#' In the legend, the subgroups will be ordered based on the order of factor levels.
#' The parameter varlist can also be a vector of multiple column names.
#' @param varlist.levels levels in the subgroups. It should be a vector if the parameter varlist is a single column name.
#' It should be a list if more than one columns are specified in the prarameter varlist.
#' The elements in the list should match the columns defined in parameter varlist.
#' Each element of the list should contain a vector, elements in the vector defines levels of the corresponding column.
#' @param varlist.labels preferred labels for the varlist.
#' varlist.levels should be provided if subgroupd.labels is specified. The order in varlist.labels should match varlist.levels.
#' It should be a vector if the parameter varlist is a single column name.
#' It should be a list if more than one columns are specified in the prarameter varlist.
#' The elements in the list should match the columns defined in parameter varlist.
#' Each element of the list should contain a vector, elements in the vector defines labels of the corresponding column.
#' @param plot.nrisk whether show number of patients at risk at the below the graph. If it is specified as TRUE, number of patients
#' at risk will be summarized by subgroup.
#' @param nrisk.interval interval to summarize number of patients at risk . Default is to summarize every 2 (months)
#' @param cex.nrisk font size for the number of patients at risk.
#' @param plot.grid whether show horizontal grids
#' @param grids horizontal grids
#' @param plot.legend whether show legend
#' @param legend.loc,legend.x,legend.y legend location. a single keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
#' @param plot.median whether show median TTE of each subgroup. (won't show if median TTE is NA)
#' @param median.cex font size of marked median. This parameter will be ignored if plot.median=FALSE
#' @param xlim,ylab,xlab,main,col,lty,lwd,sub,ylim see \code{\link{plot}}
#' @param y.percentage whether show percentage in y axis (0-100) or probability (0-1). Default is probability
#' @param return.data if it is TRUE, input data frame will be returned. If var is cont., an additional column called
#' var_group will be added to the data form, which stores the dichotomized values
#' @param surv.conf.type type of confidence interval. Default is "plain". See survfit
#' @param var.levels,var.labels parameter for old versions, please dont use
#' @param  ... additional parameters for \code{\link{plot}}
#'
#' @note This function generates KM curve(s) for full population (when parameter var is not defined)
#' or var (when parameter var is defined).
#'
#' @importFrom graphics plot axis mtext grid box polygon lines legend
#' @importFrom survival survfit
#'
#' @inheritParams CompareKM
#' @inheritParams SummarySingle
#' @inheritParams PlotTabForestBiomarker
#'
#' @examples
#' data(input)
#' sample.data <- input
#' PlotKM(data=sample.data, tte="OS",cens="OS.event", main="OS ITT by treatment", var="Arm")
#' @export

PlotKM <- function(data, tte, cens,
                   trt=NULL, var=NULL,
                   var.class=NULL,
                   var.name=NULL,
                   percentile.cutoff=0.5,quantile.type=2,cutoff.digits=2,equal.in.high = TRUE,
                   numerical.cutoff=NULL,
                   varlist=NULL, varlist.levels=NULL, varlist.labels=NULL,
                   bep=NULL, bep.indicator=1,
                   plot.nrisk=TRUE, nrisk.interval=2, cex.nrisk=.8,
                   plot.grid=TRUE, grids=seq(0,1,0.1), plot.legend=TRUE,legend.loc="topright", legend.x=NULL, legend.y=NULL,
                   col=NULL, lty=NULL, lwd=3,surv.conf.type="plain",
                   xlab="Months To Event Or Censoring", ylim=c(0,1), xlim=NULL,  ylab="Survival Probability",
                   main="",sub="", plot.median=FALSE,median.cex=.8,digits=2,y.percentage=FALSE,
                   pdf.name=NULL, pdf.param=list(height=5), par.param=list(mar=c(12,9,3,2)), return.data=FALSE,
                   var.levels=NULL, var.labels=NULL){


    if(length(var)>1){
        message("more than one elements in 'var', trt parameter will be ignored")
        varlist <- var
        varlist.levels <- var.levels
        varlist.labels <- var.labels
    }
    if(!is.null(varlist)){
        var <- trt <- NULL
        message("'varlist' is specified, trt and var parameters will be ignored")
    }

    stopifnot(class(data) == "data.frame")
    if(!is.null(bep))if(! bep %in% colnames(data))stop("bep should in column names in the input data!")
    if(!is.null(var))if(! all(c(trt,var,varlist) %in% colnames(data)))stop("names in 'var','trt','var.list' should be in column names in the input data!")
    if(!is.null(bep))data <- data[which(data[,bep]==bep.indicator),]

    var.store <- var
    percentile.footnote <- NULL
    if(!is.null(numerical.cutoff)) {
      percentile.cutoff <- NULL
      message("numerical cutoff specified")
    }
    
    ll <- c(trt,var,varlist)
    if(length(ll)>1){
        whichNA <- sapply(data[ll],function(i)which(is.na(i)))
        whichNA.v <- unique(unlist(whichNA))
        if(length(whichNA.v)>0){
            data <- data[-whichNA.v,]
            message("entries who have NA in trt, var, or varlist are removed")
        }
    }

    # in the cases when trt and var are specified - creat var.list
    if(!is.null(trt)) varlist <- trt
    if(!is.null(var)){
        possible.class <-c("categorical","numeric")
        if(is.null(var.class)||!all(var.class%in%possible.class)){
            if(class(data[,var])%in%c("numeric","integer"))var.class <- "numeric"
            if(class(data[,var])%in%c("logical"))class(data[,var]) <- "character"
            if(class(data[,var])%in%c("character","factor"))var.class <- "categorical"
        }
        data$bm.tmp <- rep(NA, length(data[[1]]))
        
        
        if(var.class=="numeric"){
            if(!is.null(percentile.cutoff)){
                percentile.cutoff <- sort(unique(c(0,1,percentile.cutoff)))
                for(i in 2:length(percentile.cutoff)){
                    qt1 <- round(quantile(data[[var]], percentile.cutoff[i-1], type=quantile.type),cutoff.digits)
                    qt2 <- round(quantile(data[[var]], percentile.cutoff[i], type=quantile.type),cutoff.digits)
                    if(equal.in.high){
                        if(percentile.cutoff[i]!=1){
                            data$bm.tmp[which(data[[var]]>=qt1 & data[[var]]< qt2)]<- paste0(var.name,"[",percentile.cutoff[i-1]*100," - ",percentile.cutoff[i]*100,"%)")
                            if(i==2)percentile.footnote <- paste0(percentile.footnote, percentile.cutoff[i-1]*100,"%: ",qt1,". ")
                            percentile.footnote <- paste0(percentile.footnote,percentile.cutoff[i]*100,"%: ",qt2,". ")
                        }
                        if(percentile.cutoff[i]==1){
                            data$bm.tmp[which(data[[var]]>=qt1 & data[[var]]<= qt2)] <- paste0(var.name,"[",percentile.cutoff[i-1]*100," - ",percentile.cutoff[i]*100,"%]")
                            if(i==2)percentile.footnote <- paste0(percentile.footnote,percentile.cutoff[i-1]*100,"%: ",qt1,". ")
                            percentile.footnote <- paste0(percentile.footnote, percentile.cutoff[i]*100,"%: ",qt2,". ")
                        }
                    }
                    if(!equal.in.high){
                        if(percentile.cutoff[i]!=0){
                            data$bm.tmp[which(data[[var]]>qt1 & data[[var]]<= qt2)] <- paste0(var.name,"(",percentile.cutoff[i-1]*100," - ",percentile.cutoff[i]*100,"%]")
                            if(i==2)percentile.footnote <- paste0(percentile.footnote, percentile.cutoff[i-1]*100,"%: ",qt1,". ")
                            percentile.footnote <- paste0(percentile.footnote,percentile.cutoff[i]*100,"%: ",qt2,". ")
                            }
                        if(percentile.cutoff[i]==0){
                            data$bm.tmp[which(data[[var]]>=qt1 & data[[var]]<= qt2)] <- paste0(var.name,"[",percentile.cutoff[i-1]*100," - ",percentile.cutoff[i]*100,"%]")
                            if(i==2)percentile.footnote <- paste0(percentile.footnote,percentile.cutoff[i-1]*100,"%: ",qt1,". ")
                            percentile.footnote <- paste0(percentile.footnote, percentile.cutoff[i]*100,"%: ",qt2,". ")
                            }
                    }

                }}

            if(!is.null(numerical.cutoff)){
                numerical.cutoff <- sort(unique(c(min(data[[var]]),max(data[[var]]),numerical.cutoff)))
                for(i in 2:length(numerical.cutoff)){
                    qt1 <- numerical.cutoff[i-1]
                    qt2 <- numerical.cutoff[i]
                    if(i==2)qt1 <- qt1 - 10^(-cutoff.digits)
                    if(i==length(numerical.cutoff)) qt2 <- qt2 + 10^(-cutoff.digits)
                    if(equal.in.high){
                        if(numerical.cutoff[i]!=max(data[[var]]))
                            data$bm.tmp[which(data[[var]]>=qt1 & data[[var]]< qt2)]  <- paste0(var.name,"[",numerical.cutoff[i-1]," - ",numerical.cutoff[i],")")
                        if(numerical.cutoff[i]==max(data[[var]]))
                            data$bm.tmp[which(data[[var]]>=qt1 & data[[var]]<= qt2)] <- paste0(var.name,"[",numerical.cutoff[i-1]," - ",numerical.cutoff[i],"]")
                    }
                    if(!equal.in.high){
                        if(numerical.cutoff[i]!=min(data[[var]]))
                            data$bm.tmp[which(data[[var]]>qt1 & data[[var]]<= qt2)] <- paste0(var.name,"(",numerical.cutoff[i-1]," - ",numerical.cutoff[i],"]")
                        if(numerical.cutoff[i]==min(data[[var]]))
                            data$bm.tmp[which(data[[var]]>=qt1 & data[[var]]<= qt2)] <- paste0(var.name,"[",numerical.cutoff[i-1]," - ",numerical.cutoff[i],"]")
                    }
                }}
            var <- paste0(var,"_groups")
            data[,var] <- data$bm.tmp
        }
        varlist <- c(varlist,var)
    }



    if(!is.null(varlist.labels) & is.null(varlist.levels)) stop("varlist.levels should be provided if varlist.labels is specified!")

    if(length(varlist)==1){
        if(!is.null(varlist.levels))if(nlevels(factor(data[[varlist]]))!=length(varlist.levels))
            stop(paste("number of elements in varlist.levels should match number of unique values in",varlist ))

        if(!is.null(varlist.labels))if(nlevels(factor(data[[varlist]]))!=length(varlist.labels))
            stop(paste("number of elements in varlist.labels should match number of unique values in",varlist ))

    }

    if(length(varlist)>1){
        if(!is.null(varlist.levels))if(length(varlist)!=length(varlist.levels))
            stop(paste("number of elements in varlist.levels should match number of column names in parameter 'varlist'"))

        if(!is.null(varlist.labels))if(length(varlist)!=length(varlist.labels))
            stop(paste("number of elements in varlist.labels should match number of column names in parameter 'varlist'"))


        for(i in 1:length(varlist)){

            if(!is.null(varlist.levels))if(nlevels(factor(data[[varlist[i]]]))!=length(varlist.levels[[i]]))
                stop(paste("number of elements in varlist.levels should match number of unique values in",varlist[i] ))

            if(!is.null(varlist.labels))if(nlevels(factor(data[[varlist[i]]]))!=length(varlist.labels[[i]]))
                stop(paste("number of elements in varlist.labels should match number of unique values in",varlist[i] ))

        }}

    var.ori <- varlist
    var.levels <- varlist.levels
    var.labels <- varlist.labels
    var <- "tmp.subgroup"
    n.subs <- length(var.ori)
    if(n.subs==0) data$tmp.subgroup <- ""
    if(n.subs==1) data$tmp.subgroup <- data[[var.ori]]
    if(n.subs>1) data$tmp.subgroup <-  apply(data[,var.ori],1,function(i)paste0(i,collapse=","))

    if(is.null(var.levels)){
        tmp.levels <- sapply(data[,var.ori],function(i)levels(factor(i)), simplify=FALSE)
        if(n.subs>1)data[,var] <- factor(data[,var], levels = apply(expand.grid(tmp.levels[n.subs:1])[,n.subs:1],1,function(i)paste0(i,collapse=",")))
        if(n.subs<=1)data[,var] <- factor(data[,var])
    }
    if(!is.null(var.levels)){
        if(length(var.ori)>1) var.levels <- apply(expand.grid(var.levels[n.subs:1])[,n.subs:1],1,function(i)paste0(i,collapse=","))
        data[,var] <- factor(data[,var], levels=var.levels)
        if(!is.null(var.labels)) {
            if(length(var.ori)>1) var.labels <- apply(expand.grid(var.labels[n.subs:1])[,n.subs:1],1,function(i)paste0(i,collapse=","))
            levels(data[,var]) <- var.labels
        }}

    if(is.null(par.param$mar))par.param$mar <- c(12,9,3,2)

    col.v <- c("blue","red","darkgreen","brown","darkgrey","skyblue","purple","cyan","pink","orange")
    strat.vec <- data[,var]
    nlev <- nlevels(strat.vec)
    if(n.subs<=1){
        if(is.null(col)) {
            col <-  col.v[1:nlev]
            if(!is.null(var.store)) col <- rep(1,nlev) # when only bm is specified, use different line types
        }
        if(is.null(lty)){
            lty <- 1
            if(!is.null(var.store)) lty <- 1:nlev # when only bm is specified, use different line types
        }
    }
    # if more than one factors, use color to distinguish first several factors and use lty to distingush the last factor
    if(n.subs>1) {
        nfirst <- length(unique(apply(data[,var.ori[-length(var.ori)], drop=FALSE],1,function(i)paste0(i,collapse=","))))
        nlast <- length(unique(data[,var.ori[length(var.ori)]]))
        if(is.null(col))col <- col.v[rep(1:nfirst, each=nlast)]
        if(is.null(lty))lty <- rep(1:nlast, nfirst)
    }

    if(is.null(var.labels))var.labels <- levels(strat.vec)

    fit <- survfit(as.formula(paste("Surv(",tte,",",cens,") ~ ", var)),conf.type=surv.conf.type, data=data)


    # xlim
    if(is.null(xlim)){
        xlim2 <- max(data[,tte],na.rm=TRUE)*1.05
        if(plot.nrisk)
            xlim1 <- -0.5
        else
            xlim1 <- 0
        xlim <- c(xlim1, xlim2)
    }


    if(plot.nrisk) {
        fit2 <- fit
        time.pt <- seq(0,xlim[2],nrisk.interval)
        ix = 0
        n.risk <- c()
        if(nlev==1) fit2$strata <- length(fit$n.risk)
        for (kk in 1:(length(fit2$strata)))
        {
            fit.n.risk = fit2$n.risk[(ix+1) : (ix+fit2$strata[kk])]
            fit.time = fit2$time[(ix+1) : (ix+fit2$strata[kk])]
            tmp = findInterval(time.pt, fit.time)
            n.risk <- rbind(n.risk, ifelse(tmp<length(fit.time), fit.n.risk[tmp+1], 0))
            ix = ix + fit2$strata[kk]
        }
        dimnames(n.risk)[[2]] = time.pt
    }

    if(plot.nrisk){
        if(par.param$mar[1] < 4+nlev) par.param$mar[1] <- 4+ nlev
    }



    if(nlev>1)meds <- summary(fit)$table[,"median"]
    if(nlev==1)meds <- summary(fit)$table["median"]

    PlotParam(pdf.name, pdf.param, par.param)

    plot(fit,col=col,lwd=lwd,xlab="", ylab=ylab,lty=lty,
         main=main, sub=sub, axes=FALSE, ylim=ylim, xlim=xlim, conf.int=F, mark.time=TRUE)
    box()

    mtext(xlab,side=1, line=2)


    axis(1,at=seq(0,xlim[2],nrisk.interval),seq(0,xlim[2],nrisk.interval))
    if(y.percentage==FALSE)axis(2,at=seq(ylim[1],ylim[2],0.1), seq(ylim[1],ylim[2],0.1),las=2); abline(h=0, col="gray")
    if(y.percentage==TRUE)axis(2,at=seq(ylim[1],ylim[2],0.1), seq(ylim[1],ylim[2],0.1)*100,las=2); abline(h=0, col="gray")


    if(plot.grid) abline(h=grids, col="gray",lty=3)
                                      
    if(plot.legend & nlev > 1){
        if(!is.null(legend.loc))legend(legend.loc,paste0(var.labels,", MST ", round(meds,digits)), lwd=2, col=col, lty=lty, bg="white")
        if(is.null(legend.loc))legend(x=legend.x, y=legend.y,paste0(var.labels,", MST ", round(meds,digits)), lwd=2, col=col, lty=lty, bg="white")
    }                                  
    
    i <- 1
    if(plot.nrisk){
        for(i in 1:nlev){
            mtext(side=1, at=xlim[1]-1.2, line=i+3,text=levels(strat.vec)[i],col=col[i],adj=1,cex=cex.nrisk*3/4)
            mtext(side=1, at=time.pt, line=i+3,text=n.risk[i,],col=col[i],cex=cex.nrisk)
        }
    }
    if(!is.null(percentile.footnote))mtext(side=1, at=xlim[1]-1.2, line=i+6,text=percentile.footnote,cex=0.8)

    if(plot.median){
        lines(c(0,max(meds, na.rm=T)),c(.5,.5), col="gray", lty=2)
        jj <- 0
        for(i in 1:nlev){
            if(!is.na(meds[i])){
                text(x=meds[i],y=0.05+(ylim[1]+(diff(ylim)/10)*jj),labels=paste0(var.labels[i],"\nmedian ",round(meds[i], digits)), col=col[i],cex=median.cex)
                lines(c(meds[i],meds[i]), c(0,.5),lty=3, lwd=1, col=col[i])
                jj <- jj+1
            }
            if(is.na(meds[i])){
                text(x=xlim[1]+(diff(xlim)/10),y=0.05+(ylim[1]+(diff(ylim)/10)*jj),labels=paste0(var.labels[i],"\nmedian NA"), col=col[i],cex=median.cex)
                jj <- jj+1
            }
        }
    }

    if (length(all_labels()) == 0) {
        PlotParam()
    }

    out <- ""
    if(return.data)out <- data
    out
}
