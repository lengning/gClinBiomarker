#' Generate forest plot and summarization table for 2-arm or within-arm comparison for multiple variables
#'
#' This function creates a forest plot along with table with summary statistics to infer effect of multiple clinical variables, within a single arm
#' or across two treatment arms.  The outcome could be survival, binary or continuous. This function can be used to summarize a group of
#' variables. This function may be used to compare effect of these variables in full population vs biomarker evaluable population. Or compare effect
#' of these variables in different biomarker subgroups
#'
#' @author  Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param compare.subgroup whether want to compare across multiple subgroups, 
#' in each subgroup defined by parameter 'subgroup'. Alternative is to compare BEP vs ITT
#' @param subgroup The column which defines subgroups. If compare.subgroup is TRUE, the program will generate forest plot of the vars within each subgroup
#' @param compare.bep.itt whether want to compare BEP vs ITT
#'
#' @inheritParams PlotTabForestBiomarker
#' @export
#' @examples
#' data(input)
#' PlotTabForestMulti(data=input,
#'                       outcome.class=c("survival"),
#'                       outcome.var=c("PFS","PFS.event"),
#'                       trt="Arm",
#'                       var=c("Sex","Age"),
#'                       var.class="categorical",bep="BEP")



PlotTabForestMulti <- function(data,
                                  outcome.class=c("survival", "binary","continuous"),
                                  outcome.var, #c(OS,OS.event)
                                  trt=NULL,
                                  var, #KRAS...
                                  var.class=NULL, var.name=NULL,
                                  percentile.cutoff=0.5,
                                  greater=TRUE, less=TRUE,
                                  within.bin=FALSE,compare.bep.itt=TRUE, compare.subgroup=FALSE,
                                  show.itt=FALSE, show.bep=FALSE,
                                  subgroup=NULL,
                                  bep = NULL, bep.name = "BEP", itt.name="All",bep.indicator=1,
                                  covariate=NULL, #Sex
                                  strata=NULL, #Age
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  quantile.type=2,
                                  placebo.code=NULL,
                                  active.code=NULL,
                                  var.code=NULL,
                                  tabforest=FALSE,
                                  alpha=0.05,
				  surv.conf.type="plain",
				  ties="efron",
                                  main=NULL,
                                  sub=NULL,
                                  clip=NULL,
                                  xlab=NULL,
                                  cex.headings=1.1,
                                  cex.note=1,
                                  cols="darkgreen",
                                  only.stat=FALSE,
                                  pdf.name=NULL,
                                  pdf.param=list(width=6, height=4.5),
                                  par.param=list(cex=1.2, cex.main=1.5, cex.sub=1, cex.axis=1)) {

  if(compare.bep.itt & compare.subgroup) stop("compare.bep.itt & compare.subgroup cannot both be true!")
  if(compare.bep.itt)if(is.null(bep))stop("compare.bep.itt is TRUE, bep needs to be specified!")

  if(!is.null(subgroup)){
    stopifnot(subgroup%in%colnames(data))
    data[[subgroup]] <- factor(data[[subgroup]])
    groups.level <- levels(data[[subgroup]])
    ngroups <- nlevels(data[[subgroup]])
    if(any(greater, less)) {
      greater<- TRUE
      less <- TRUE
    }
  }
  data.list <- list(ITT=data)
  names(data.list)[1] <- itt.name
  if(compare.bep.itt){
    data.list <- list(ITT=data, BEP=data[which(data[[bep]]%in%bep.indicator),])
    names(data.list)[2] <- bep.name
    }
  if(compare.subgroup){
    data.list <- sapply(groups.level, function(i)data[which(data[[subgroup]]==i),], simplify=F)
    if(show.bep){
      if(is.null(bep)){
        if(any(is.na(data[[subgroup]])))message("show.bep is TRUE but bep is not specified, will define the non NA entries in subgroup column as BEP")
        data$BEPnew <- ifelse(is.na(data[[subgroup]]),0,bep.indicator)
        bep <- "BEPnew"
      }
      data.list <- c(list(BEP=data[which(data[[bep]]%in%bep.indicator),]), data.list)
    }
    if(show.itt) data.list <- c(list(ITT=data), data.list)
    }

  if(length(unique(data[,trt]))==1) {
    trt <- NULL
  }
  if(is.null(trt)) nArms <- 1
  if (!is.null(trt)) { # multi-arm study
    Treatment <- data[,trt]
    Arms <- levels(factor(Treatment))
    if (!is.null(placebo.code)) { # reference arm specified
      if(length(placebo.code)!=1)stop("placebo.code should have length 1")
      if(!placebo.code %in% Arms)stop("placebo.code should be an element in treatment column")
      Arms <- c(placebo.code, setdiff(Arms, placebo.code)) # first placebo then others
      if (!is.null(active.code)) {
        if(length(active.code)!=1)stop("active.code should have length 1")
        if(!all(active.code %in% Arms))stop("active.code should be elements in treatment column")
        if(length(intersect(placebo.code, active.code))>1)
          stop("code cannot be in both active.code and placebo.code!")
        Arms <- c(Arms[1], active.code) #order by specified input
      }
    }
    nArms <- length(Arms)
    if(!nArms%in%c(1,2))stop("only 1 or 2 arms are allowed")
    data[,trt] <- factor(data[,trt],levels=Arms)
    if(is.null(placebo.code))placebo.code <- Arms[1]
    if(is.null(active.code))active.code <- Arms[-1]
  }



  res.list <- sapply(data.list, function(dd){
    sapply(1:length(var), function(v){
    var.class.tmp <- var.class
    if(!is.null(var.class))var.class.tmp <- var.class[v]
    var.name.tmp <- var.name
    if(!is.null(var.name))var.name.tmp <- var.name[v]
    PlotTabForestBiomarker (data=dd,
                                       outcome.class=outcome.class,
                                       outcome.var=outcome.var,
                                       trt=trt,
                                       var=var[v],
                                       var.class=var.class.tmp,
                                       var.name=var.name.tmp,
                                       percentile.cutoff=percentile.cutoff,
                                       numerical.cutoff=NULL,
                                       greater=greater, less=less,
                                       within.bin=within.bin,
                                       show.itt=FALSE, show.bep=FALSE,
                                       bep = NULL, bep.name = "BEP", itt.name="ITT",bep.indicator=1,
                                       rsp.cat = rsp.cat,
                                       rsp.response = rsp.response,
                                       rsp.nonresponse = rsp.nonresponse,
                                       covariate=covariate, 
                                       strata=strata,
				       surv.conf.type=surv.conf.type, ties=ties,
                                       tabforest=tabforest,
                                       quantile.type=quantile.type,
                                       placebo.code=placebo.code,
                                       active.code=active.code,
                                       alpha=alpha,
                                       only.stat=TRUE)},simplify=FALSE)
  },simplify=FALSE)

  res.list2 <- res.list
  for(i in 1:length(res.list)){
    for(j in 1:length(var)){
      nlev <- (nrow(res.list2[[i]][[j]])-1)/2
      res.list2[[i]][[j]][(1:nlev)*2,1] <- paste(names(res.list)[i], res.list[[i]][[j]][(1:nlev)*2,1])
  }}

  final.tab <- res.list[[1]][[1]][1,]
  hl <- NULL
  ct <- 0
  for(j in 1:length(var)){
    nlev <- (nrow(res.list2[[1]][[j]])-1)/2
    for(k in 1:nlev){
    for(i in 1:length(res.list)){
        ct <- ct+1
        final.tab <- rbind(final.tab,res.list2[[i]][[j]][c(2+(k-1)*2, 1+(k*2)),])
      }
      }
  hl <- c(hl, nrow(final.tab)-1)
    }
  tabletext <- final.tab

  if (is.null(main)) {
    main.text <- ifelse(nArms==1, "Within arm", "Across arm")
    if(compare.bep.itt)main.text <- paste0(main.text, ", Compare ",bep.name, " vs. " , itt.name )
    if(compare.subgroup) main.text <- paste0(main.text, ", Compare ",subgroup," subgroup")
    main.text <- paste0(main.text, "\n", outcome.var[1])

    } else {
    main.text <- main
  }

  if (is.null(sub)) {
    sub1.text <- NULL
    if(length(covariate) > 0)sub1.text <- paste("Results adjusted by ", paste(covariate, collapse=" , "), sep="")
    sub2.text <- NULL
    if(length(strata) > 0)
      sub2.text <-  paste("Results stratified by ", paste(strata, collapse=" , "), sep="")
    if (is.null(sub1.text)  & is.null(sub2.text) ) {
      sub.text <- "Unadjusted, unstratified analysis"
    } else {
      sub.text <- paste(sub1.text, sub2.text, sep=";")
    }
  } else {
    sub.text <- sub
  }
 if(!only.stat){
  PlotParam(pdf.name, pdf.param, par.param)

  if (is.null(clip) & outcome.class=="survival") {
    r1 <- as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][1]))
    r2 <- as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][2]))
    good1 <- !is.na(r1) & is.finite(r1) & r1!= 0
    good2 <- !is.na(r2) & is.finite(r2)
    xrange <- c(min(round(r1[good1], 2)), max(as.numeric(round(r2[good2], 2))))
    clip <- exp(c(-max(abs(log(xrange))), max(abs(log(xrange)))))
  }
  if (is.null(clip) & outcome.class=="binary") {
    r1 <- as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][1]))
    r2 <- as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][2]))
    good1 <- !is.na(r1) & is.finite(r1) & r1!= 0
    good2 <- !is.na(r2) & is.finite(r2)
    xrange <- c(min(round(r1[good1], 2)), max(as.numeric(round(r2[good2], 2))))
    mm <- max(abs(xrange))
    clip <- c(-mm,mm)
  }
  if (is.null(clip) & outcome.class=="continuous") {
    r1 <- as.numeric(sapply(tabletext[-1, 5], function(z)strsplit(z, " - ")[[1]][1]))
    r2 <- as.numeric(sapply(tabletext[-1, 5], function(z)strsplit(z, " - ")[[1]][2]))
    good1 <- !is.na(r1) & is.finite(r1) & r1!= 0
    good2 <- !is.na(r2) & is.finite(r2)
    xrange <- c(min(round(r1[good1], 2)), max(as.numeric(round(r2[good2], 2))))
    mm <- max(abs(xrange))
    clip <- c(-mm,mm)
  } 

  wid <- max(nchar(sapply(tabletext[,1], function(z)strsplit(z, "\n")[[1]][1])),na.rm=T)/6

  note <- ""
  if(length(cols)==nrow(tabletext)/2) cols <- rep(cols,each=2)

  if(tabforest){
      if(is.null(xlab)) {
          if(nArms==2)xlab <- c(paste(active.code, "better", sep=" "),
          paste(placebo.code, "better", sep=" "))
          if(nArms==1)xlab <- c("","")
      }
    
    xlog <- TRUE
    if(outcome.class=="binary") xlog <- FALSE
    num1 <- 5
    num2 <- 6
    wid2 <- c( wid,2, 1.5, 1, 1, 2, 1, 5)
    if(outcome.class=="continuous"){
      num1 <- 4
      num2 <- 5
      wid2 <- c( wid,2, 1.5,  1, 2, 1, 5)
    }
   
        
    PlotTabForest(label.text=tabletext[-c(1), ],
                mean=as.numeric(tabletext[-1, num1]),
                lower=as.numeric(sapply(tabletext[-1, num2], function(z)strsplit(z, " - ")[[1]][1])),
                upper=as.numeric(sapply(tabletext[-1, num2], function(z)strsplit(z, " - ")[[1]][2])),
                headings=c(tabletext[1, ], c("Forest plot")),
                cols=cols,
                xlog=xlog,
                xticks=NULL,
                box.size=rep(2.5, nrow(tabletext)-1),
                main=main.text,
                sub=sub.text,
                hline=1:nrow(tabletext),
                vline=c(1),
                group.hline=hl,
                note=note,clip=clip,
                widths=c( wid,2, 1.5, 1, 1, 2, 1, 5),
                sub.main=xlab,
                cex.headings=cex.headings,
                cex.note=cex.note,
                par.param=par.parm
  )
    }

  if(!tabforest){
    num1 <- 5
    num2 <- 6
    if(outcome.class=="continuous"){
      num1 <- 4
      num2 <- 5
    }
      hz <- vector("list",1)
      for(i in 1:length(hl)){
        if(hl[i] < nrow(tabletext)){
          hz[[i]] <- gpar(lwd=2, col="#99999999")
          names(hz)[i] <- hl[i]+2
        }
      }

      tabletext2 <- tabletext
      tabletext2[seq(1,nrow(tabletext2),2),num2] <- paste0("(",tabletext2[seq(1,nrow(tabletext2),2),num2],")")
    
    if(is.null(xlab)) {
        if(nArms==2)xlab <- paste("<-- ", active.code, "better [",res.list[[1]][[1]][1,num1],"] ",placebo.code, "better -->\n",note)
        if(nArms==1)xlab <- res.list[[1]][[1]][1,num1]
    }
      xlog <- TRUE
      if(outcome.class=="binary") {
          if(nArms==2)xlab <- paste("<-- ", placebo.code, "better [",res.list[[1]][[1]][1,num1],"] ",active.code, "better -->\n",note)
          xlog <- FALSE
      }



          forestplot(tabletext2,
                 mean=c(NA,as.numeric(tabletext[-1,num1])),
                 lower=c(NA,as.numeric(sapply(tabletext[-1, num2], function(z)strsplit(z, " - ")[[1]][1]))),
                 upper=c(NA,as.numeric(sapply(tabletext[-1, num2], function(z)strsplit(z, " - ")[[1]][2]))),
                 xlab=xlab,
                 hrzl_lines=hz,align="l",
                 lwd.xaxis=2, lwd.ci=2,col=fpColors(box=cols, line=cols),
                 xlog=xlog,
                 title=paste(main.text,"\n",sub.text),
                 #graphwidth=unit(100, 'mm'),
                 colgap=unit(cex.note*4,"mm"),
                 line.margin =unit(cex.note*2,"mm"),
                 txt_gp=fpTxtGp(label=gpar(cex=cex.note),
                                ticks=gpar(cex=cex.note),
                                xlab=gpar(cex = cex.note))
      )

      }
  PlotParam()
}
  out <- tabletext
}















