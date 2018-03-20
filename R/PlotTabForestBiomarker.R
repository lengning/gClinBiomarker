#' Generate forest plot and summarization table for 2-arm or within-arm comparison
#'
#' This function creates a forest plot along with table with summary statistics to infer biomarker effects, within a single arm
#' or across two treatment arms.  The outcome could be survival, binary or continuous. This function can be used to summarize a single
#' biomarker variable
#'
#' @author  Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param data input data frame. Rows are patients and columns are variables (e.g. demographics variables, time to event variables,
#' biomarker variables, treatment indicator, etc.). One patient per row.
#' @param outcome.class type of the outcome variable. Default is \code{c("survival", "binary", "continuous")}. Continuous is not available now
#' @param outcome.var name of the outcome varible. If the outcome.class is binary or coutinuous, only one value should be provided.
#' If the outcome.class is survival, two values should be provided - name of the 'time to event' variable and 'censorship' variable
#'  For the censoring variable, 1 indicates event and 0 indicates censoring. In all cases but when outcome.class=binary and rsp.cat=TRUE,
#'  patients with missing outcome variable (NA) will be excluded from BEP.
#' @param trt name of the treatment variable. If this is NULL, within-arm analysis will be performed
#' @param var name of the biomarker variable. only one variable should be specified.
#' @param var.class class of the variable. valid categories are "numeric", "categorical". If the class is continuous,
#' user needs to specify percentile.cutoff or numerical.cutoff to dichotomize the continuous measure into subgroups
#' @param var.name display name for the biomarker variable
#' @param percentile.cutoff percentile to dichotomize continuous biomarker measure. This could be a vector with multiple elements.
#' Values should be between 0 and 1
#' @param numerical.cutoff raw value to dichotomize continuous biomarker measure. numerical.cutoff and percentile.cutoff
#' cannot be both specified
#' @param greater whether calculate summary statistics within the subgroup whose biomarker value is greater than or equal to
#' cutoff value. If this is TRUE, in 2-arm study, across-arm HR within biomarker high group will be calculated.
#' In single arm study HR of biomarker high vs low will be calculated.
#' @param less whether calculate summary statistics within the subgroup whose biomarker value is less than the cutoff value.
#' greater and less can both be TRUE
#' @param greater.by.less whether show "greater" bin and "less" bin in consecutive rows. Default is FALSE. If it is TRUE,
#' parameters greater and less will be both set as TRUE
#' @param across.and.within whether show across- and within- arm results in the same figure. Default is FALSE. This parameter
#' will be ignored if number of arm is 1. If it is TRUE, within-arm analysis results will be shown below the across-arm results.
#' @param equal.in.high whether include equal in high group. Default is TRUE. If it is TRUE, ">=" and "<" will be
#' applied. Otherwise "<=" and ">" will be applied.
#' @param within.bin whether calculate summary statistics within bin (e.g. > cutoff1 and <= cutoff2). If within.bin is TRUE,
#' greater and less will be set as FALSE.
#' @param show.itt whether calculate summary statistics using all patients in full population (e.g. ITT). This will be ignored in 1arm case
#' @param show.bep whether calculate summary statistics using all patients in BEP (biomarker evaluable population). This will be ignored in 1arm case
#' @param bep name of the column which indicates biomarker evaluable population. If it is null, patients who have non NA records
#' in biomarker variable will be used as BEP.
#' @param bep.name preferred display name of the biomarker evaluable population.
#' If it is NULL, bep will be used.
#' @param itt.name preferred display name of the full population (e.g. ITT).
#' If it is NULL, "All" will be used.
#' @param bep.indicator In the subpopulation column, which value is used
#' to define the biomarker evaluable population.
#' @param covariate a vector specifying the covariate variables to be adjusted in the model. Default is set to NULL, meaning no adjustment.
#' @param strata name of the stratification variables. Default is set to NULL, meaning no stratification.
#' @param placebo.code name of the control arm of the treatment variable. If you want to specify placebo code using this parameter, both placebo.code and active.code need to be provided.
#' @param active.code of the treatment/experimental arm of the treatment variable. If you want to specify active code using this parameter, both placebo.code and active.code need to be provided.
#' @param var.code ordered levels of the biomarker variable. This will be ignored for continuous biomarker.
#' If the biomarker is categorical and this is NULL, biomarker subgroups will be ordered by the order from factor() function
#' @param rsp.cat whether the response outcome variable is coded as binary (1 as responder and 0 as non-responder),
#' If rsp.cat is TRUE, responder categories
#' and nonresponder categories should be specified in rsp.response and rsp.nonresponse (all values in the outcome column
#' should be included in rsp.response and rsp.nonresponse)
#' . If rsp.cat is FALSE, the response outcome variable should be coded as binary (0/1).
#' At the same time rsp.response and rsp.nonresponse will be ignored.
#' @param rsp.response categories that should be considered as responder.
#' @param rsp.nonresponse categories that should be considered as non responder.
#' @param tabforest Default is FALSE. If it is FALSE, forest plot will be generated using forestplot::forestplot() function.
#' If it is TRUE, a table will be generated with forest plots incorpriated
#' @param quantile.type an integer between 1 and 9 selecting one of the nine quantile algorithms. See \code{\link{quantile}}. Default is 2.
#' @param alpha type I error rate. Default is 0.05.
#' @param ties Default is "efron". To match internal sas results, use "exact". See parameter "ties" in coxph.
#' @param surv.conf.type confidence interval type. Default is "plain". see conf.type in survfit
#' @param cutoff.digits,digits cutoff.digits:number of digits for rounding when calculating cutoff. will only be used when percentile.cutoff is specified. digits: number of digits for the summary statistics display
#' @param main,main.prefix main title (prefix of title) of the forest plot. Default is "Association of biomarker effect within treatment arms".
#' @param sub sub title under the forest plot. Default is NULL.
#' @param clip range of the x-axis of the forest plot. Default is NULL.
#' @param xticks,xticks.digits x axis tick marks for the forest plot
#' @param xlab xlab for forest plot
#' @param cex.headings amount of magnification of headings of the forest plot relative to cex. Default is 1.1.
#' @param cex.note amount of magnification of the note. Default is 1.
#' @param cols Color of the 'effect size' displayed in the forest plot.
#' @param only.stat if it is TRUE, only summary statistics will be generated. No figure will be generated
#' @param pdf.name name of output pdf file. If it's NULL, the plots will be displayed but not saved as pdf. Default is "forestplot::forestplot.pdf".
#' @param pdf.param a list of parameters that define pdf graphics device. See \code{\link{pdf}}. Default is \code{list(width=6, height=4.5)}.
#' @param par.param a list of parameters that define graphcial parameters. See \code{\link{par}}. Default is \code{list(mar=c(4,4,3,2))}.
#'
#' @importFrom grid gpar
#' @importFrom forestplot forestplot fpTxtGp fpColors
#' @export
#'
#' @examples
#' data(input)
#' PlotTabForestBiomarker(data=input,
#'                       outcome.class=c("survival"),
#'                       outcome.var=c("PFS","PFS.event"),
#'                       trt="Arm",
#'                       var="KRAS.mutant",
#'                       var.class="categorical")



PlotTabForestBiomarker <- function(data,
outcome.class=c("survival", "binary"),
outcome.var, #c(OS,OS.event)
trt=NULL,
var=NULL, #KRAS...
var.class=NULL,
var.name=NULL,
percentile.cutoff=0.5,
numerical.cutoff=NULL,
greater=TRUE,
less=FALSE,
greater.by.less = FALSE,
across.and.within = FALSE,
equal.in.high = TRUE,
within.bin=FALSE,
show.itt=TRUE,
show.bep=TRUE,
bep = NULL,
bep.name = "BEP",
itt.name="All",
bep.indicator=1,
covariate=NULL, #Sex
strata=NULL, #Age
tabforest=FALSE,
quantile.type=2,
digits=2, cutoff.digits=2,
placebo.code=NULL,
active.code=NULL,
rsp.cat = TRUE,
rsp.response = c("CR","PR"),
rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
var.code=NULL,
surv.conf.type="plain",
ties="efron",
alpha=0.05,
main=NULL,
main.prefix=NULL,
sub=NULL,
clip=NULL, xticks=NULL, xticks.digits=1,
xlab=NULL,
cex.headings=1.1,
cex.note=.8,
cols=NULL,
only.stat=FALSE,
pdf.name=NULL,
pdf.param=list(width=12, height=4.5),
par.param=list(cex=1, cex.main=1, cex.sub=1, cex.axis=1)) {

    stopifnot(class(data) == "data.frame")
    percentile.footnote <- NULL


    outcome.class <- match.arg(outcome.class, c("survival", "binary","continuous"))
    if(outcome.class=="binary") {
        covariate <- strata <- NULL
        message("Covariate adjustment and stratification are not supported for binary outcome")
    }
    if(outcome.class=="continuous") {
        covariate <- strata <- NULL
        message("Stratification is not supported for continuous outcome")
    }

    if(outcome.class == "survival" | (outcome.class=="binary" & rsp.cat==F)) {
        tmp <- rowMeans(data[,outcome.var])
        if (anyNA(tmp)) {
            message("Some patients have missing outcome. Exclude these patients from ITT.")
            data <- data[which(!is.na(tmp)),]
        }
    }
    if(is.null(var))show.itt <- TRUE

    possible.class <-c("categorical","numeric")
    if(is.null(var.class)||!all(var.class%in%possible.class)){
        if(class(data[,var])%in%c("numeric","integer"))var.class <- "numeric"
        if(class(data[,var])%in%c("logical"))class(data[,var]) <- "character"
        if(class(data[,var])%in%c("character","factor"))var.class <- "categorical"
    }

    var.class <- match.arg(var.class,c("numeric","categorical"))
    stopifnot(var%in%colnames(data))
    if(is.null(bep)){
        if(any(is.na(data[[var]])))message("Some NAs in var column, will define the non NA entries as BEP")
        data$BEPnew <- ifelse(is.na(data[[var]]),0,bep.indicator)
        bep <- "BEPnew"
    }
    if(!is.null(percentile.cutoff) & !is.null(numerical.cutoff)) stop("Cannot specify both percentile.cutoff and numerical.cutoff")

    Outcome <- data[, outcome.var]
    Biomarker <- data[, var]



    if(var.class=="categorical"){
        if (!is.null(var.code)) {
            if (length(var.code) == length(levels(Biomarker))) {
                Biomarker <- factor(Biomarker, levels=var.code)
            } else {
                stop("Length of var.code needs to be the same of levels of Biomarker...\n")
            }
        } else {
            var.code <- levels(Biomarker)
        }}



    if(!is.null(trt))if(length(unique(data[,trt]))==1)  trt <- NULL

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

    if(nArms==1){
        if(show.itt){
            show.itt <- FALSE
            message("only 1 arm; show.itt is set to FALSE")
        }
        if(show.bep){
            show.bep <- FALSE
            message("only 1 arm; show.bep is set to FALSE")
        }
        if(within.bin){
            message("only 1 arm; within.bin is set to FALSE")
            within.bin <- FALSE
        }
    }

    if(within.bin & any(c(greater, less))){
        message("within.bin is TRUE, greater and less will be ignored")
        greater <- less <- FALSE
    }

    if(greater.by.less) {
        greater <- less <- TRUE
        within.bin <- FALSE
    }

    if(is.null(cols)){
        cols <- "mediumblue"
        if(within.bin) cols <- "darkorchid"
        if(nArms==1) cols <- "chocolate4"
    }

    if(is.null(var.name))var.name <- var
    if(is.null(bep.name))bep.name <- "BEP"
    if(is.null(itt.name))itt.name <- "All"
    ncut <- max(length(percentile.cutoff), length(numerical.cutoff))
    if(var.class=="numeric" & ncut==0)stop("numeric var but no cutoff was specified!")


    data.bep <- data[which(data[[bep]]%in%bep.indicator),]
    bm.list <- list()
    #if(show.itt)bm.list[[itt.name]] <- rep(T, length(data.bep[[1]]))
    if(show.bep&nArms==2) bm.list[[bep.name]] <- ifelse(data.bep[[bep]]%in%bep.indicator,T,F)
    if(!is.null(var)){
        if(var.class=="categorical"){
            b.class <- names(table(data.bep[,var]))
            for(i in b.class)bm.list[[paste0(var.name,"(",i,")")]] <- ifelse(data.bep[[var]]==i,T,F)
            ncut <- length(b.class)
        }
        if(var.class=="numeric"){
            if(any(is.na(data.bep[[var]])))stop(paste("in BEP patients," ,var,"contains NA"))

            if(!greater.by.less){
                if(greater)if(!is.null(percentile.cutoff)) for(i in percentile.cutoff){
                    qt <- round(quantile(data.bep[[var]], i, type=quantile.type),cutoff.digits)
                    if(equal.in.high)bm.list[[paste0(var.name,"(>=",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]>=qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(>",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]>qt,T,F)
                }
                if(less)if(!is.null(percentile.cutoff)) for(i in percentile.cutoff){
                    qt <- round(quantile(data.bep[[var]], i, type=quantile.type),cutoff.digits)
                    if(equal.in.high)bm.list[[paste0(var.name,"(<",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]<qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(<=",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]<=qt,T,F)
                }
                if(greater)if(!is.null(numerical.cutoff)) for(i in numerical.cutoff){
                    qt <- i
                    if(equal.in.high)bm.list[[paste0(var.name,"(>=",i,")")]] <- ifelse(data.bep[[var]]>=qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(>",i,")")]] <- ifelse(data.bep[[var]]>qt,T,F)
                }
                if(less)if(!is.null(numerical.cutoff)) for(i in numerical.cutoff){
                    qt <- i
                    if(equal.in.high)bm.list[[paste0(var.name,"(<",i,")")]] <- ifelse(data.bep[[var]]<qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(<=",i,")")]] <- ifelse(data.bep[[var]]<=qt,T,F)
                }}

            if(greater.by.less){
                if(!is.null(percentile.cutoff)) for(i in percentile.cutoff){
                    qt <- round(quantile(data.bep[[var]], i, type=quantile.type),cutoff.digits)
                    if(equal.in.high)bm.list[[paste0(var.name,"(>=",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]>=qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(>",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]>qt,T,F)
                    if(equal.in.high)bm.list[[paste0(var.name,"(<",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]<qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(<=",i*100,"%, ",qt,")")]] <- ifelse(data.bep[[var]]<=qt,T,F)
                }
                if(!is.null(numerical.cutoff)) for(i in numerical.cutoff){
                    qt <- i
                    if(equal.in.high)bm.list[[paste0(var.name,"(>=",i,")")]] <- ifelse(data.bep[[var]]>=qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(>",i,")")]] <- ifelse(data.bep[[var]]>qt,T,F)
                    if(equal.in.high)bm.list[[paste0(var.name,"(<",i,")")]] <- ifelse(data.bep[[var]]<qt,T,F)
                    if(!equal.in.high)bm.list[[paste0(var.name,"(<=",i,")")]] <- ifelse(data.bep[[var]]<=qt,T,F)
                }
            }

            if(within.bin)if(!is.null(percentile.cutoff)){
                percentile.cutoff <- sort(unique(c(0,1,percentile.cutoff)))
                for(i in 2:length(percentile.cutoff)){
                    qt1 <- round(quantile(data.bep[[var]], percentile.cutoff[i-1], type=quantile.type),cutoff.digits)
                    qt2 <- round(quantile(data.bep[[var]], percentile.cutoff[i], type=quantile.type),cutoff.digits)
                    if(i==2)qt1 <- qt1 - 10^(-cutoff.digits)
                    if(i==length(percentile.cutoff)) qt2 <- qt2 + 10^(-cutoff.digits)
                    if(equal.in.high){
                        if(percentile.cutoff[i]!=1){
                            bm.list[[paste0(var.name,"[",percentile.cutoff[i-1]*100,"-",percentile.cutoff[i]*100,"%, ",qt1,"-",qt2,")")]]  <-
                            ifelse(data.bep[[var]]>=qt1 & data.bep[[var]]< qt2,T,F)
                        }
                        if(percentile.cutoff[i]==1){
                            bm.list[[paste0(var.name,"[",percentile.cutoff[i-1]*100,"-",percentile.cutoff[i]*100,"%, ",qt1,"-",qt2,"]")]] <-
                            ifelse(data.bep[[var]]>=qt1 & data.bep[[var]]<= qt2,T,F)
                        }
                    }
                    if(!equal.in.high){
                        if(percentile.cutoff[i]!=0){
                            bm.list[[paste0(var.name,"(",percentile.cutoff[i-1]*100,"-",percentile.cutoff[i]*100,"%, ",qt1,"-",qt2,"]")]]  <-
                            ifelse(data.bep[[var]]>qt1 & data.bep[[var]]<= qt2,T,F)
                        }
                        if(percentile.cutoff[i]==0){
                            bm.list[[paste0(var.name,"[",percentile.cutoff[i-1]*100,"-",percentile.cutoff[i]*100,"%, ",qt1,"-",qt2,"]")]] <-
                            ifelse(data.bep[[var]]>=qt1 & data.bep[[var]]<= qt2,T,F)
                        }
                    }

                }}

            if(within.bin)if(!is.null(numerical.cutoff)){
                numerical.cutoff <- sort(unique(c(min(data.bep[[var]]),max(data.bep[[var]]),numerical.cutoff)))
                for(i in 2:length(numerical.cutoff)){
                    qt1 <- numerical.cutoff[i-1]
                    qt2 <- numerical.cutoff[i]
                    if(i==2)qt1 <- qt1 - 10^(-cutoff.digits)
                    if(i==length(numerical.cutoff)) qt2 <- qt2 + 10^(-cutoff.digits)
                    if(equal.in.high){
                        if(i!=length(numerical.cutoff))
                        bm.list[[paste0(var.name,"[",numerical.cutoff[i-1],"-",numerical.cutoff[i],")")]] <-
                        ifelse(data.bep[[var]]>=qt1 & data.bep[[var]]< qt2,T,F)
                        if(i==length(numerical.cutoff))
                        bm.list[[paste0(var.name,"[",numerical.cutoff[i-1],"-",numerical.cutoff[i],"]")]] <-
                        ifelse(data.bep[[var]]>=qt1 & data.bep[[var]]<= qt2,T,F)
                    }
                    if(!equal.in.high){
                        if(i!=2)
                        bm.list[[paste0(var.name,"(",numerical.cutoff[i-1],"-",numerical.cutoff[i],"]")]] <-
                        ifelse(data.bep[[var]]>qt1 & data.bep[[var]]<= qt2,T,F)
                        if(i==2)
                        bm.list[[paste0(var.name,"[",numerical.cutoff[i-1],"-",numerical.cutoff[i],"]")]] <-
                        ifelse(data.bep[[var]]>=qt1 & data.bep[[var]]<= qt2,T,F)
                    }
                }}
        }

    }


    #### stratification,
    # one for ITT and one for BEP

    if(is.null(covariate)) {
        Covariate <- Covariate.bep <- NULL
    } else {
        Covariate.bep <- data.bep[, covariate]
        Covariate <- data[,covariate]
    }

    if(is.null(strata)) {
        Strat.fac <- Strat.fac.bep <- NULL
    } else {
        Strat.fac.bep <- data.bep[, strata]
        Strat.fac <- data[,strata]
    }


    ################## survival #########################
    if(outcome.class=="survival"){
        res <- NULL
        if(nArms==2){
            if(!is.null(var)) res <- t(sapply(bm.list,function(jj)SummaryTwoGroups(outcome.var=data.bep[,outcome.var],
            subgroup.var=jj, treatment.var=data.bep[,trt],
            placebo.code=placebo.code, active.code=active.code, outcome.class="survival", alpha=alpha,surv.conf.type=surv.conf.type, ties=ties,
            covariate.var=Covariate.bep,
            strat.factor.var=Strat.fac.bep)))

            if(show.itt) {res <- rbind(
                SummaryTwoGroups(outcome.var=data[,outcome.var],
                subgroup.var=rep(T, length(data[[1]])), treatment.var=data[,trt],
                placebo.code=placebo.code, active.code=active.code, outcome.class="survival", alpha=alpha,
		surv.conf.type=surv.conf.type, ties=ties,
                covariate.var=Covariate,
                strat.factor.var=Strat.fac)
                ,res)
                rownames(res)[1] <- itt.name
            }

            ac <- Arms[2]
            # interaction p value: per arm? If originally cont., then use cont. in modeling
            # no stratification??
            fit1 <- coxph(Surv(data.bep[,outcome.var[1]], data.bep[,outcome.var[2]]) ~ data.bep[,trt] * data.bep[,var] )
            fit2 <- coxph(Surv(data.bep[,outcome.var[1]], data.bep[,outcome.var[2]]) ~ data.bep[,trt] + data.bep[,var])
            L1 <- summary(fit1)[[5]][2]
            n1 <- summary(fit1)[[9]][2]
            L2 <- summary(fit2)[[5]][2]
            n2 <- summary(fit2)[[9]][2]
            stat <- -2*L2 + 2*L1
            inter.p <- pchisq(stat, df=n1-n2, lower.tail=FALSE)
            if(max(c(length(percentile.cutoff), length(numerical.cutoff)))>1 ) inter.p<- NULL

            if(across.and.within){
                res.ori <- res
                data.bep.soc <- data.bep[which(data.bep[[trt]]==placebo.code),]
                bm.list.subonly <- bm.list[setdiff(names(bm.list),c(itt.name, bep.name))]
                bm.list.soc <- data.frame(bm.list)[which(data.bep[[trt]]==placebo.code),]
                bm.list.soc <- bm.list.soc[setdiff(names(bm.list.soc),c(itt.name, bep.name))]
                if(is.null(covariate))Covariate.bep.soc <- NULL
                else Covariate.bep.soc <- data.bep.soc[, covariate]
                if(is.null(strata))Strat.fac.bep.soc <- NULL
                else Strat.fac.bep.soc <- data.bep.soc[, strata]
                res.soc <- t(sapply(bm.list.soc,function(jj)SummaryTwoGroups(outcome.var=data.bep.soc[,outcome.var],
                subgroup.var=rep(T,length(data.bep.soc[[1]])), treatment.var=jj,
                placebo.code="FALSE", active.code="TRUE", outcome.class="survival", alpha=alpha,
		surv.conf.type=surv.conf.type, ties=ties,
                covariate.var=Covariate.bep.soc,
                strat.factor.var=Strat.fac.bep.soc)))
                rownames(res.soc) <- paste0(placebo.code,":", names(bm.list.subonly))

                data.bep.active <- data.bep[which(data.bep[[trt]]==active.code),]
                bm.list.active <- data.frame(bm.list)[which(data.bep[[trt]]==active.code),]
                bm.list.active <- bm.list.active[setdiff(names(bm.list.active),c(itt.name, bep.name))]
                if(is.null(covariate))Covariate.bep.active <- NULL
                else Covariate.bep.active <- data.bep.active[, covariate]
                if(is.null(strata))Strat.fac.bep.active <- NULL
                else Strat.fac.bep.active <- data.bep.active[, strata]
                res.active <- t(sapply(bm.list.active,function(jj)SummaryTwoGroups(outcome.var=data.bep.active[,outcome.var],
                subgroup.var=rep(T,length(data.bep.active[[1]])), treatment.var=jj,
                placebo.code="FALSE", active.code="TRUE", outcome.class="survival", alpha=alpha,
		surv.conf.type=surv.conf.type, ties=ties,
                covariate.var=Covariate.bep.active,
                strat.factor.var=Strat.fac.bep.active)))
                rownames(res.active) <- paste0(active.code,":", names(bm.list.subonly))

                res <- rbind(res.ori, res.soc, res.active)
            }


        }

        if(nArms==1){
            placebo.code <- ""
            active.code <- ""
            res <- t(sapply(bm.list,function(jj)SummaryTwoGroups(outcome.var=data.bep[,outcome.var],
            subgroup.var=rep(T,length(data.bep[[1]])), treatment.var=jj,
            placebo.code="FALSE", active.code="TRUE", outcome.class="survival", alpha=alpha,
	    surv.conf.type=surv.conf.type, ties=ties,
            covariate.var=Covariate.bep,
            strat.factor.var=Strat.fac.bep)))
            inter.p <- NULL
        }

        if(nArms==2 & !across.and.within)code.v <- rep(c(placebo.code, active.code),nrow(res))
        if(nArms==2 & across.and.within) code.v.ori <- rep(c(placebo.code, active.code),nrow(res.ori))
        if(nArms==1 | across.and.within) {
            if(nArms==1)nn <- names(bm.list)
            if(across.and.within) nn <- names(bm.list.subonly)
            vv <- c("less","greater")
            code.l <- sapply(nn, function(i){
                j1 <- grep("<",i)
                j2 <- grep(">",i)
                j <- c(j1, j2)
                if(length(j1)==1)out <- c("Greater","Less")
                if(length(j2)==1) out <- c("Less","Greater")
                if(length(j)==0)out <- c("No","Yes")
                out})
            code.v <- as.vector(code.l)
            if(across.and.within) code.v <- c(code.v.ori,code.v, code.v)
        }

        tabletext <- rbind(c( "Subgroup","Group", "Event/N", "MST", "HR", "CI", "raw P"),
        cbind(as.vector(sapply(rownames(res),function(z)c(z, ""))),
        code.v,
        as.vector(t(cbind(paste(res[, 1], "/", res[, 2], sep=" "),
        paste(res[, 4], "/", res[, 5], sep=" ")))),
        as.vector(t(round(res[, c(3, 6)], 2))),
        as.vector(t(cbind(rep("", nrow(res)), round(res[, 7], digits)))),
        as.vector(t(cbind(rep("", nrow(res)), paste(round(res[, 8], digits), round(res[, 9], digits), sep=" - ")))),
        as.vector(t(cbind(rep("", nrow(res)), round.signif(res[, 10], digits))))))


        if(!only.stat){
            if (is.null(main)) {
                main.text <- ifelse(nArms==1, "Within-arm Effect of Biomarker", "Across-arm Effect of Biomarker")
                if(across.and.within) main.text <- "Across-arm and Within-arm Effect of Biomarker"
                main.text <- paste0(main.prefix, " ",main.text, "\n", outcome.var[1],", ",var.name)
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

            PlotParam(pdf.name, pdf.param, par.param)

            if (is.null(clip)) {
                good1 <- !is.na(res[, 8]) & is.finite(res[, 8]) & res[, 8] != 0
                good2 <- !is.na(res[, 9]) & is.finite(res[, 9])
                xrange <- c(min(round(res[good1, 8], digits)), max(as.numeric(round(res[good2, 9], digits))))
                mm <- max(max(abs(log(xrange))),1)
                clip <- exp(c(-mm,mm))
            }
            if(is.null(xticks))xticks <- round(c(exp(seq(log(clip[1]),log(clip[2]),length.out=5)),1),digits)


            wid <- max(nchar(sapply(rownames(res), function(z)strsplit(z, "\n")[[1]][1])))/6

            if(within.bin) ncut <- ncut+1
            hl <- 0
            if(show.itt) hl <- c(hl,max(hl)+2)
            if(show.bep) hl <- c(hl, max(hl)+2)
            # if(within.bin) hl <- c(hl, max(hl)+length(bm.list)*2)
            if(!greater.by.less){
                if(greater) hl <- c(hl, max(hl)+ncut*2)
                if(less) hl <- c(hl, max(hl)+ncut*2)
            }
            if(greater.by.less) hl <- c(hl, max(hl)+ncut*4)
            if(across.and.within & !greater.by.less) hl <- c(hl, max(hl)+ncut*2, max(hl)+ncut*4)
            if(across.and.within & greater.by.less) hl <- c(hl, max(hl)+ncut*4, max(hl)+ncut*8)


            note <- ""
            if(length(cols)==nrow(tabletext)/2) cols <- rep(cols,each=2)
            if(!is.null(inter.p)) note <- paste0("unadj P = ", round.signif(inter.p, 2), "(interaction)")
            if(tabforest){
                if(is.null(xlab)) {
                    if(nArms==2)xlab <- c(paste(active.code, "better", sep=" "),
                    paste(placebo.code, "better", sep=" "))
                    if(nArms==1)xlab <- c("","")
                    if(across.and.within) xlab <- c("","")
                }
                PlotTabForest(label.text=tabletext[-c(1), ],
                mean=as.numeric(tabletext[-1, 5]),
                lower=as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][1])),
                upper=as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][2])),
                headings=c(tabletext[1, ], c("Forest plot")),
                cols=cols,
                xlog=TRUE,
                xticks=xticks,
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

                hz <- vector("list",1)
                for(i in 1:length(hl)){
                    if(hl[i] < nrow(tabletext)){
                        hz[[i]] <- grid::gpar(lwd=2, col="#444444")
                        names(hz)[i] <- hl[i]+2
                    }
                }

                tabletext2 <- tabletext
                tabletext2[seq(1,nrow(tabletext2),2),6] <- paste0("(",tabletext2[seq(1,nrow(tabletext2),2),6],")")
                if(is.null(xlab)) {
                    if(nArms==2)xlab <- paste("<-- ", active.code, "better [HR] ",placebo.code, "better -->\n",note)
                    if(nArms==1 | across.and.within)xlab <- "HR"
                }
                if(across.and.within) {
                    cols <- rep(cols, nrow(res))
                    cols[(length(bm.list)+1): length(cols)] <- "chocolate4"
                }
                forestplot::forestplot(tabletext2,
                mean=c(NA,as.numeric(tabletext[-1,5])),
                lower=c(NA,as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][1]))),
                upper=c(NA,as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][2]))),
                xlab=xlab,
                hrzl_lines=hz,align="l",
                lwd.xaxis=2, lwd.ci=2,col=forestplot::fpColors(box=cols, line=cols),
                lwd.zero=3,
                xlog=TRUE, clip=clip, xticks=xticks,
                title=paste(main.text,"\n",sub.text),
                #graphwidth=unit(100, 'mm'),
                colgap=unit(cex.note*4,"mm"),
                line.margin =unit(cex.note*2,"mm"),
                txt_gp=forestplot::fpTxtGp(label=grid::gpar(cex=cex.note),
                ticks=grid::gpar(cex=cex.note),
                xlab=grid::gpar(cex = cex.note))
                )

            }

            PlotParam()
        }
        out <- tabletext
    }


    ################## binary #########################
    if(outcome.class=="binary"){

        if(!rsp.cat)if(!all(data[,outcome.var]%in%c(0,1))){
            stop("rsp.cat is FALSE, all elements in outcome.var should be 0 or 1!")
        }
        if(rsp.cat)if(!all(unique(data[,outcome.var])%in%c(rsp.response,rsp.nonresponse)))
        stop("all unique values in outcome.var column should be included in rsp.response or rsp.nonresponse!")

        outcome.var.ori <- outcome.var
        # generate response
        if(rsp.cat){
            data$rspvar <- ifelse(data[,outcome.var]%in%rsp.response,1,0)
            data.bep$rspvar <- ifelse(data.bep[,outcome.var]%in%rsp.response,1,0)
            outcome.var <- 'rspvar'
        }


        res <- NULL
        if(nArms==2){
            if(!is.null(var)) res <- t(sapply(bm.list,function(jj)SummaryTwoGroups(outcome.var=data.bep[,outcome.var],
            subgroup.var=jj, treatment.var=data.bep[,trt],
            placebo.code=placebo.code, active.code=active.code, outcome.class="binary", alpha=alpha,
            covariate.var=Covariate.bep,
            strat.factor.var=Strat.fac.bep)))

            if(show.itt) {res <- rbind(
                SummaryTwoGroups(outcome.var=data[,outcome.var],
                subgroup.var=rep(T, length(data[[1]])), treatment.var=data[,trt],
                placebo.code=placebo.code, active.code=active.code, outcome.class="binary", alpha=alpha,
                covariate.var=Covariate,
                strat.factor.var=Strat.fac)
                ,res)
                rownames(res)[1] <- itt.name
            }

            ac <- Arms[2]
            # interaction p value: per arm? If originally cont., then use cont. in modeling
            # no stratification??
            fit1 <- glm(data[,outcome.var]~as.character(Treatment)*Biomarker, subset=as.character(Treatment) %in% c(placebo.code, ac), family=binomial)
            fit2 <- glm(data[,outcome.var]~as.character(Treatment)+Biomarker, subset=as.character(Treatment) %in% c(placebo.code, ac), family=binomial)
            L1 <- summary(fit1)$deviance
            n1 <- summary(fit1)[[7]]
            L2 <- summary(fit2)$deviance
            n2 <- summary(fit2)[[7]]
            stat <- L2-L1
            inter.p <- pchisq(stat, df=n2-n1, lower.tail=FALSE)
            if(max(c(length(percentile.cutoff), length(numerical.cutoff)))>1 ) inter.p<- NULL

            if(across.and.within){
                res.ori <- res
                data.bep.soc <- data.bep[which(data.bep[[trt]]==placebo.code),]
                bm.list.subonly <- bm.list[setdiff(names(bm.list),c(itt.name, bep.name))]
                bm.list.soc <- data.frame(bm.list)[which(data.bep[[trt]]==placebo.code),]
                bm.list.soc <- bm.list.soc[setdiff(names(bm.list.soc),c(itt.name, bep.name))]
                if(is.null(covariate))Covariate.bep.soc <- NULL
                else Covariate.bep.soc <- data.bep.soc[, covariate]
                if(is.null(strata))Strat.fac.bep.soc <- NULL
                else Strat.fac.bep.soc <- data.bep.soc[, strata]
                res.soc <- t(sapply(bm.list.soc,function(jj)SummaryTwoGroups(outcome.var=data.bep.soc[,outcome.var],
                subgroup.var=rep(T,length(data.bep.soc[[1]])), treatment.var=jj,
                placebo.code="FALSE", active.code="TRUE", outcome.class="binary", alpha=alpha,
                covariate.var=Covariate.bep.soc,
                strat.factor.var=Strat.fac.bep.soc)))
                rownames(res.soc) <- paste0(placebo.code,":", names(bm.list.subonly))

                data.bep.active <- data.bep[which(data.bep[[trt]]==active.code),]
                bm.list.active <- data.frame(bm.list)[which(data.bep[[trt]]==active.code),]
                bm.list.active <- bm.list.active[setdiff(names(bm.list.active),c(itt.name, bep.name))]
                if(is.null(covariate))Covariate.bep.active <- NULL
                else Covariate.bep.active <- data.bep.active[, covariate]
                if(is.null(strata))Strat.fac.bep.active <- NULL
                else Strat.fac.bep.active <- data.bep.active[, strata]
                res.active <- t(sapply(bm.list.active,function(jj)SummaryTwoGroups(outcome.var=data.bep.active[,outcome.var],
                subgroup.var=rep(T,length(data.bep.active[[1]])), treatment.var=jj,
                placebo.code="FALSE", active.code="TRUE", outcome.class="binary", alpha=alpha,
                covariate.var=Covariate.bep.active,
                strat.factor.var=Strat.fac.bep.active)))
                rownames(res.active) <- paste0(active.code,":", names(bm.list.subonly))

                res <- rbind(res.ori, res.soc, res.active)
            }


        }

        if(nArms==1){
            placebo.code <- ""
            active.code <- ""
            res <- t(sapply(bm.list,function(jj)SummaryTwoGroups(outcome.var=data.bep[,outcome.var],
            subgroup.var=rep(T,length(data.bep[[1]])), treatment.var=jj,
            placebo.code="FALSE", active.code="TRUE", outcome.class="binary", alpha=alpha,
            covariate.var=Covariate.bep,
            strat.factor.var=Strat.fac.bep)))
            inter.p <- NULL
        }

        if(nArms==2 & !across.and.within)code.v <- rep(c(placebo.code, active.code),nrow(res))
        if(nArms==2 & across.and.within) code.v.ori <- rep(c(placebo.code, active.code),nrow(res.ori))
        if(nArms==1 | across.and.within) {
            if(nArms==1)nn <- names(bm.list)
            if(across.and.within) nn <- names(bm.list.subonly)
            vv <- c("less","greater")
            code.l <- sapply(nn, function(i){
                j1 <- grep("<",i)
                j2 <- grep(">",i)
                j <- c(j1, j2)
                if(length(j1)==1)out <- c("Greater","Less")
                if(length(j2)==1) out <- c("Less","Greater")
                if(length(j)==0)out <- c("No","Yes")
                out})
            code.v <- as.vector(code.l)
            if(across.and.within) code.v <- c(code.v.ori,code.v, code.v)
        }

        tabletext <- rbind(c( "Subgroup","Group", "nRsp/N", "Rsp Rate", "deltaRR", "CI", "raw P"),
        cbind(as.vector(sapply(rownames(res),function(z)c(z, ""))),
        code.v,
        as.vector(t(cbind(paste(res[, 9], "/", res[, 7], sep=" "),
        paste(res[, 10], "/", res[, 8], sep=" ")))),
        as.vector(t(round(res[,c(5, 6)], 2))),
        as.vector(t(cbind(rep("", nrow(res)), round(res[, 1], 2)))),
        as.vector(t(cbind(rep("", nrow(res)), paste(round(res[, 2], 2), round(res[, 3], 2), sep=" - ")))),
        as.vector(t(cbind(rep("", nrow(res)), round.signif(res[, 4], 2))))))

        if(!only.stat){
            if (is.null(main)) {
                main.text <- ifelse(nArms==1, "Within-arm Effect of Biomarker", "Across-arm Effect of Biomarker")
                if(across.and.within) main.text <- "Across-arm and Within-arm Effect of Biomarker"
                main.text <- paste0(main.prefix, " ",main.text, "\n", outcome.var.ori[1],", ",var.name)
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

            PlotParam(pdf.name, pdf.param, par.param)

            if (is.null(clip)) {
                xrange <- c(min(round(res[, 2], digits)), max(as.numeric(round(res[, 3], digits))))
                mm <- max(abs(xrange))
                clip <- c(-mm,mm)
            }
            if(is.null(xticks))xticks <- round(seq(clip[1],clip[2],length.out=5),digits)


            wid <- max(nchar(sapply(rownames(res), function(z)strsplit(z, "\n")[[1]][1])))/6

            if(within.bin) ncut <- ncut+1
            hl <- 0
            if(show.itt) hl <- c(hl,max(hl)+2)
            if(show.bep) hl <- c(hl, max(hl)+2)
            # if(within.bin) hl <- c(hl, max(hl)+length(bm.list)*2)
            if(!greater.by.less){
                if(greater) hl <- c(hl, max(hl)+ncut*2)
                if(less) hl <- c(hl, max(hl)+ncut*2)
            }
            if(greater.by.less) hl <- c(hl, max(hl)+ncut*4)
            if(across.and.within & !greater.by.less) hl <- c(hl, max(hl)+ncut*2, max(hl)+ncut*4)
            if(across.and.within & greater.by.less) hl <- c(hl, max(hl)+ncut*4, max(hl)+ncut*8)


            note <- ""
            if(length(cols)==nrow(tabletext)/2) cols <- rep(cols,each=2)
            if(!is.null(inter.p)) note <- paste0("* Unadj P = ", paste(round.signif(inter.p, 2),"(interaction)", collapse=" ; "))
            if(tabforest){
                if(is.null(xlab)) {
                    if(nArms==2)xlab <- c(paste(placebo.code, "better", sep=" "),
                    paste(active.code, "better", sep=" "))
                    if(nArms==1)xlab <- c("","")
                    if(across.and.within) xlab <- c("","")
                }
                PlotTabForest(label.text=tabletext[-c(1), ],
                mean=as.numeric(tabletext[-1, 5]),
                lower=as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][1])),
                upper=as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][2])),
                headings=c(tabletext[1, ], c("Forest plot")),
                cols=cols,
                xlog=FALSE,
                xticks=xticks,
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

                hz <- vector("list",1)
                for(i in 1:length(hl)){
                    if(hl[i] < nrow(tabletext)){
                        hz[[i]] <- grid::gpar(lwd=2, col="#444444")
                        names(hz)[i] <- hl[i]+2
                    }
                }

                tabletext2 <- tabletext
                tabletext2[seq(1,nrow(tabletext2),2),6] <- paste0("(",tabletext2[seq(1,nrow(tabletext2),2),6],")")
                if(is.null(xlab)) {
                    if(nArms==2)xlab <- paste("<-- ", placebo.code, "better [deltaRR] ",active.code, "better -->\n",note)
                    if(nArms==1 | across.and.within)xlab <- "deltaRR"
                }
                if(across.and.within) {
                    cols <- rep(cols, nrow(res))
                    cols[(length(bm.list)+1): length(cols)] <- "chocolate4"
                }
                forestplot::forestplot(tabletext2,
                mean=c(NA,as.numeric(tabletext[-1,5])),
                lower=c(NA,as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][1]))),
                upper=c(NA,as.numeric(sapply(tabletext[-1, 6], function(z)strsplit(z, " - ")[[1]][2]))),
                xlab=xlab,
                hrzl_lines=hz,align="l",
                lwd.xaxis=2, lwd.ci=2,col=forestplot::fpColors(box=cols, line=cols),
                lwd.zero=3,
                xlog=FALSE, clip=clip, xticks=xticks,
                title=paste(main.text,"\n",sub.text),
                #graphwidth=unit(100, 'mm'),
                colgap=unit(cex.note*4,"mm"),
                line.margin =unit(cex.note*2,"mm"),
                txt_gp=forestplot::fpTxtGp(label=grid::gpar(cex=cex.note),
                ticks=grid::gpar(cex=cex.note),
                xlab=grid::gpar(cex = cex.note))
                )

            }

            PlotParam()
        }
        out <- tabletext
    }

    if(outcome.class=="continuous")  {
        res <- NULL
        if(nArms==2){
            if(!is.null(var)) res <- t(sapply(bm.list,function(jj)SummaryTwoGroups(outcome.var=data.bep[,outcome.var],
            subgroup.var=jj, treatment.var=data.bep[,trt],
            placebo.code=placebo.code, active.code=active.code,
            outcome.class="continuous", alpha=alpha,
            covariate.var=Covariate.bep,
            strat.factor.var=Strat.fac.bep)))

            if(show.itt) {res <- rbind(
                SummaryTwoGroups(outcome.var=data[,outcome.var],
                subgroup.var=rep(T, length(data[[1]])), treatment.var=data[,trt],
                placebo.code=placebo.code, active.code=active.code, outcome.class="continuous", alpha=alpha,
                covariate.var=Covariate,
                strat.factor.var=Strat.fac)
                ,res)
                rownames(res)[1] <- itt.name
            }

            ac <- Arms[2]
            # no stratification??
            fit1 <- lm(data[,outcome.var]~as.character(Treatment)*Biomarker, subset=as.character(Treatment) %in% c(placebo.code, ac))
            inter.p <- coef(summary(fit1))[4,4]
            if(max(c(length(percentile.cutoff), length(numerical.cutoff)))>1 ) inter.p<- NULL

            if(across.and.within){
                res.ori <- res
                data.bep.soc <- data.bep[which(data.bep[[trt]]==placebo.code),]
                bm.list.subonly <- bm.list[setdiff(names(bm.list),c(itt.name, bep.name))]
                bm.list.soc <- data.frame(bm.list)[which(data.bep[[trt]]==placebo.code),]
                bm.list.soc <- bm.list.soc[setdiff(names(bm.list.soc),c(itt.name, bep.name))]
                if(is.null(covariate))Covariate.bep.soc <- NULL
                else Covariate.bep.soc <- data.bep.soc[, covariate]
                if(is.null(strata))Strat.fac.bep.soc <- NULL
                else Strat.fac.bep.soc <- data.bep.soc[, strata]
                res.soc <- t(sapply(bm.list.soc,function(jj)SummaryTwoGroups(outcome.var=data.bep.soc[,outcome.var],
                subgroup.var=rep(T,length(data.bep.soc[[1]])), treatment.var=jj,
                placebo.code="FALSE", active.code="TRUE", outcome.class="continuous", alpha=alpha,
                covariate.var=Covariate.bep.soc,
                strat.factor.var=Strat.fac.bep.soc)))
                rownames(res.soc) <- paste0(placebo.code,":", names(bm.list.subonly))

                data.bep.active <- data.bep[which(data.bep[[trt]]==active.code),]
                bm.list.active <- data.frame(bm.list)[which(data.bep[[trt]]==active.code),]
                bm.list.active <- bm.list.active[setdiff(names(bm.list.active),c(itt.name, bep.name))]
                if(is.null(covariate))Covariate.bep.active <- NULL
                else Covariate.bep.active <- data.bep.active[, covariate]
                if(is.null(strata))Strat.fac.bep.active <- NULL
                else Strat.fac.bep.active <- data.bep.active[, strata]
                res.active <- t(sapply(bm.list.active,function(jj)SummaryTwoGroups(outcome.var=data.bep.active[,outcome.var],
                subgroup.var=rep(T,length(data.bep.active[[1]])), treatment.var=jj,
                placebo.code="FALSE", active.code="TRUE", outcome.class="continuous", alpha=alpha,
                covariate.var=Covariate.bep.active,
                strat.factor.var=Strat.fac.bep.active)))
                rownames(res.active) <- paste0(active.code,":", names(bm.list.subonly))

                res <- rbind(res.ori, res.soc, res.active)
            }


        }

        if(nArms==1){
            placebo.code <- ""
            active.code <- ""
            res <- t(sapply(bm.list,function(jj)SummaryTwoGroups(outcome.var=data.bep[,outcome.var],
            subgroup.var=rep(T,length(data.bep[[1]])), treatment.var=jj,
            placebo.code="FALSE", active.code="TRUE", outcome.class="continuous", alpha=alpha,
            covariate.var=Covariate.bep,
            strat.factor.var=Strat.fac.bep)))
            inter.p <- NULL
        }

        if(nArms==2 & !across.and.within)code.v <- rep(c(placebo.code, active.code),nrow(res))
        if(nArms==2 & across.and.within) code.v.ori <- rep(c(placebo.code, active.code),nrow(res.ori))
        if(nArms==1 | across.and.within) {
            if(nArms==1)nn <- names(bm.list)
            if(across.and.within) nn <- names(bm.list.subonly)
            vv <- c("less","greater")
            code.l <- sapply(nn, function(i){
                j1 <- grep("<",i)
                j2 <- grep(">",i)
                j <- c(j1, j2)
                if(length(j1)==1)out <- c("Greater","Less")
                if(length(j2)==1) out <- c("Less","Greater")
                if(length(j)==0)out <- c("No","Yes")
                out})
            code.v <- as.vector(code.l)
            if(across.and.within) code.v <- c(code.v.ori,code.v, code.v)
        }

        tabletext <- rbind(c( "Subgroup","Group", "Mean", "delta", "CI", "raw P"),
        cbind(as.vector(sapply(rownames(res),function(z)c(z, ""))),
        code.v,
        as.vector(t(round(res[,c(5, 6)], 2))),
        as.vector(t(cbind(rep("", nrow(res)), round(res[, 1], 2)))),
        as.vector(t(cbind(rep("", nrow(res)), paste(round(res[, 2], 2), round(res[, 3], 2), sep=" - ")))),
        as.vector(t(cbind(rep("", nrow(res)), round.signif(res[, 4], 2))))))

        if(!only.stat){
            if (is.null(main)) {
                main.text <- ifelse(nArms==1, "Within-arm Effect of Biomarker", "Across-arm Effect of Biomarker")
                if(across.and.within) main.text <- "Across-arm and Within-arm Effect of Biomarker"
                main.text <- paste0(main.prefix, " ",main.text, "\n", outcome.var[1],", ",var.name)
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

            PlotParam(pdf.name, pdf.param, par.param)

            if (is.null(clip)) {
                xrange <- c(min(round(res[, 2], digits)), max(as.numeric(round(res[, 3], digits))))
                mm <- max(abs(xrange))
                clip <- c(-mm,mm)
            }
            if(is.null(xticks))xticks <- round(seq(clip[1],clip[2],length.out=5),digits)


            wid <- max(nchar(sapply(rownames(res), function(z)strsplit(z, "\n")[[1]][1])))/6

            if(within.bin) ncut <- ncut+1
            hl <- 0
            if(show.itt) hl <- c(hl,max(hl)+2)
            if(show.bep) hl <- c(hl, max(hl)+2)
            # if(within.bin) hl <- c(hl, max(hl)+length(bm.list)*2)
            if(!greater.by.less){
                if(greater) hl <- c(hl, max(hl)+ncut*2)
                if(less) hl <- c(hl, max(hl)+ncut*2)
            }
            if(greater.by.less) hl <- c(hl, max(hl)+ncut*4)
            if(across.and.within & !greater.by.less) hl <- c(hl, max(hl)+ncut*2, max(hl)+ncut*4)
            if(across.and.within & greater.by.less) hl <- c(hl, max(hl)+ncut*4, max(hl)+ncut*8)


            note <- ""
            if(length(cols)==nrow(tabletext)/2) cols <- rep(cols,each=2)
            if(!is.null(inter.p)) note <- paste0("* Unadj P = ", paste(round.signif(inter.p, 2), "interaction",collapse=" ; "))
            if(tabforest){
                if(is.null(xlab)) {
                    if(nArms==2)xlab <- c(paste(placebo.code, "better", sep=" "),
                    paste(active.code, "better", sep=" "))
                    if(nArms==1)xlab <- c("","")
                    if(across.and.within) xlab <- c("","")
                }
                PlotTabForest(label.text=tabletext[-c(1), ],
                mean=as.numeric(tabletext[-1, 4]),
                lower=as.numeric(sapply(tabletext[-1, 5], function(z)strsplit(z, " - ")[[1]][1])),
                upper=as.numeric(sapply(tabletext[-1, 5], function(z)strsplit(z, " - ")[[1]][2])),
                headings=c(tabletext[1, ], c("Forest plot")),
                cols=cols,
                xlog=FALSE,
                xticks=xticks,
                box.size=rep(2.5, nrow(tabletext)-1),
                main=main.text,
                sub=sub.text,
                hline=1:nrow(tabletext),
                vline=c(1),
                group.hline=hl,
                note=note,clip=clip,
                widths=c( wid,2, 1.5, 1, 2, 1, 5),
                sub.main=xlab,
                cex.headings=cex.headings,
                cex.note=cex.note,
                par.param=par.parm
                )
            }
            if(!tabforest){

                hz <- vector("list",1)
                for(i in 1:length(hl)){
                    if(hl[i] < nrow(tabletext)){
                        hz[[i]] <- grid::gpar(lwd=2, col="#444444")
                        names(hz)[i] <- hl[i]+2
                    }
                }

                tabletext2 <- tabletext
                tabletext2[seq(1,nrow(tabletext2),2),5] <- paste0("(",tabletext2[seq(1,nrow(tabletext2),2),5],")")
                if(is.null(xlab)) {
                    if(nArms==2)xlab <- paste("<-- ", placebo.code, "better [delta] ",active.code, "better -->\n",note)
                    if(nArms==1 | across.and.within)xlab <- "delta"
                }
                if(across.and.within) {
                    cols <- rep(cols, nrow(res))
                    cols[(length(bm.list)+1): length(cols)] <- "chocolate4"
                }
                forestplot::forestplot(tabletext2,
                mean=c(NA,as.numeric(tabletext[-1,4 ])),
                lower=c(NA,as.numeric(sapply(tabletext[-1, 5], function(z)strsplit(z, " - ")[[1]][1]))),
                upper=c(NA,as.numeric(sapply(tabletext[-1, 5], function(z)strsplit(z, " - ")[[1]][2]))),
                xlab=xlab,
                hrzl_lines=hz,align="l",
                lwd.xaxis=2, lwd.ci=2,col=forestplot::fpColors(box=cols, line=cols),
                lwd.zero=3,
                xlog=FALSE, clip=clip, xticks=xticks,
                title=paste(main.text,"\n",sub.text),
                #graphwidth=unit(100, 'mm'),
                colgap=unit(cex.note*4,"mm"),
                line.margin =unit(cex.note*2,"mm"),
                txt_gp=forestplot::fpTxtGp(label=grid::gpar(cex=cex.note),
                ticks=grid::gpar(cex=cex.note),
                xlab=grid::gpar(cex = cex.note))
                )

            }

            PlotParam()
        }
        out <- tabletext}
    out
}
