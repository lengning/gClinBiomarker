#' Summary outputs from cox proportional model
#'
#' This function summarizes outputs from cox proportional model of one or multiple variables.
#'
#' @author Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param var a vector that indicates variables of interest.
#' Hazard ratio, CI and pvalue of these variables will be calculated. All elements in var should be also in input data's column names.
#' It will be ignored if fit or form is not null.
#' @param strata a vector that indicates stratification factors. All elements should be in input data's column names
#' strata will be ignored if var is NULL. It will be ignored if form or fit is not NULL
#' @param form a formula that includes variables of interest. It will be ignored if fit is not null.
#' For example, 'Age+Sex' or "Age+Sex*Trt+Country"
#' @param fit a fitted object from coxph() function.
#' If fit is not NULL, all other parameters will be ignored
#' @param additive Whether an additive model should be used. If additive=FALSE, separate cox PH models will be fitted to each elements in var
#' (stratification factor will be incorpriated for each model).
#' If additive = FALSE, User can only use 'var' to speficy variables of interest (cannot specify the model via parameter 'form')

#' @note The function generates a table that contains hazard ratio, CI, wald p value, and number of patients in each sub category.
#' User may input column names (via var and strata), a formula (via form), or a fitted coxph object (via fit)
#' If user chooses to input column names, a additive model will be formed if additive=TRUE. If additive=FALSE, separate cox PH models will be formed for each var.
#' More coplex models may be specified using form (e.g. model with interactions).
#' To calculate log rank test p value across different subgroups, see LogRankTab()
#'
#' @importFrom stats as.formula complete.cases fisher.test kruskal.test sd
#' @importFrom coin cmh_test pvalue
#'
#' @inheritParams SummaryVars
#' @inheritParams PlotTabForestBiomarker
#' @inheritParams PlotKM
#' 
#'
#' @examples
#' data(input)
#' sample.data <- input
#' CoxTab(data=sample.data,tte="OS", cens="OS.event",  var=c('Sex',"Country","Age"))
#' CoxTab(data=sample.data,tte="OS", cens="OS.event",  form="Age+Sex*Arm+Country")
#'
#' @export

CoxTab <- function(data=NULL, tte=NULL, cens=NULL, var=NULL, var.class=NULL, ordered.factor.levels.list=NULL, strata=NULL, form=NULL, fit=NULL,
		   additive=TRUE, digits=2,
		     bep = NULL, bep.indicator=1){
  fit0 <- fit
  if(is.null(fit) & is.null(form) & is.null(var)) stop("st least one of var, form, fit need to be not null")
  if(additive==FALSE & !is.null(form))stop("form cannot be specified if additive is FALSE! Please specify var")

  if(is.null(fit)){ # If fit is not NULL, all other variables will be ignored
  stopifnot(class(data) == "data.frame")
  if(!is.null(bep))data <- data[which(data[[bep]]==bep.indicator),]
  if(is.null(form)) {
	form <- paste0(var, collapse="+")
    	if(!is.null(strata)) form <- paste(form, "+", paste(paste0("strata(",strata,")"), collapse="+"))
    }
    if(is.null(var)) var <- all.vars(as.formula(paste("~",form)))
    form <- as.formula(paste("Surv(",tte,",",cens,")~", form)) # if form not NULL, var will be ignored

  if(!all(c(var, tte, cens) %in% colnames(data)))stop("bep, var, tte and cens should have matched column names in the input data!")
 ## variable class
  possible.class <-c("categorical","numeric","ordered.factor")

  if(!all(var.class%in%possible.class))stop(paste('var.class should be in', paste(possible.class,collapse=",")))

  if(length(var)== length(var.class)) names(var.class) <- var
  if(length(var)!= length(var.class)) if(!all(names(var.class)%in%var))stop("length of var.class doesn't match length of var -
			In this case the program request that names of the vector var.class is a subset of the var vector")
  if(is.null(var.class)) var.class <- rep("", length(var))
  if(!all(var.class%in%possible.class)){
    for(vv in 1:length(var)){
	if(class(data[,var[vv]])%in%c("numeric","integer"))var.class[vv] <- "numeric"
    	if(class(data[,var[vv]])%in%c("logical"))data[,var[vv]] <- "character"
    	if(class(data[,var[vv]])%in%c("character","factor"))var.class[vv] <- "categorical"
  }}
  if(!all(var.class%in%possible.class))stop(paste('var.class should be in', paste(possible.class,collapse=",")))

  for(vv in 1:length(var)){
	 if(var.class[vv]=="ordered.factor"){
	          if(is.null(ordered.factor.levels.list[[var[vv]]])) stop(paste(var[vv],": if class is ordered.factor,
				  ordered.factor.levels need to be specified."))
	 if(!is.null(ordered.factor.levels.list[[var[vv]]]))
		 ordered.factor.levels<- as.character(ordered.factor.levels.list[[var[vv]]])
		if(!identical(sort(ordered.factor.levels),sort(levels(factor(data[,var[vv]])))))stop(paste(var[vv],": ordered factor levels should match unique elements in the variable!"))
		  data[,var[vv]] <- factor(data[,var[vv]],levels=ordered.factor.levels, ordered=TRUE)
  }
  	if(var.class[vv]=="categorical")data[,var[vv]] <- factor(data[,var[vv]])
  }

## Generate 'fit' when fit is null
     if(additive){
	 fit <- coxph(form,data=data)

	  model <- summary(fit)
	  hr <- round(model$coefficients[,"exp(coef)",drop=F],digits)
	  hrci.l <- round(model$conf.int[,"lower .95",drop=F],digits)
	  hrci.h <- round(model$conf.int[,"upper .95",drop=F],digits)
	  res <- cbind(hr, hrci.l, hrci.h,signif(model$coefficients[,"Pr(>|z|)",drop=F],digits))
	  colnames(res) <- c("HR","CI.low","CI.high","p-value")
}

  if(!additive){
     	fit.list <- sapply(var, function(vv){
		  form <- vv
		  if(!is.null(strata))form <- paste(form, "+", paste(paste0("strata(",strata,")"), collapse="+"))
                  form <- as.formula(paste("Surv(",tte,",",cens,")~", form))
	     	  fit <- coxph(form,data=data)
		  model <- summary(fit)
		  hr <- round(model$coefficients[,"exp(coef)",drop=F],digits)
		  hrci.l <- round(model$conf.int[,"lower .95",drop=F],digits)
		  hrci.h <- round(model$conf.int[,"upper .95",drop=F],digits)
		  res <- cbind(hr, hrci.l, hrci.h,signif(model$coefficients[,"Pr(>|z|)",drop=F],digits))
		  colnames(res) <- c("HR","CI.low","CI.high","p-value")
    		  list(res=res, xlevels=fit$xlevels)
},simplify=F)
  res <- do.call(rbind, sapply(fit.list, function(jj)jj$res))
  fit <- vector("list",1)
  fit$xlevels <- sapply(fit.list, function(jj)jj$xlevels[[1]])
  }
}

if(!is.null(fit0)){

   	model <- summary(fit)
	  hr <- round(model$coefficients[,"exp(coef)",drop=F],digits)
	  hrci.l <- round(model$conf.int[,"lower .95",drop=F],digits)
	  hrci.h <- round(model$conf.int[,"upper .95",drop=F],digits)
	  res <- cbind(hr, hrci.l, hrci.h,signif(model$coefficients[,"Pr(>|z|)",drop=F],digits))
	  colnames(res) <- c("HR","CI.low","CI.high","p-value")
}

  nfactors = 0
  for(i in 1:length(var) ){
    nfactors = nfactors + as.numeric( is.factor(data[,var][[i]]))
  }

  if( nfactors > 0 ) {

    vars <- fit$xlevels
    ninter <- grep(":", rownames(res))
    res0 <- res
    res <- cbind(res, "", "")
    colnames(res)[ncol(res)-1:0] <- c("n.trt","n.ref")
    for(i in 1:length(vars)){
      if(is.null(vars[[i]]))next
      nn <- names(vars)[i]
      whichi <- setdiff(grep(nn, rownames(res)), ninter)
      rownames(res)[whichi] <- paste(nn," (", vars[[i]][-1], "/", vars[[i]][1], ")",sep="")
      res[whichi, c("n.trt")] <- sapply(vars[[i]][-1],function(j)length(which(data[[names(vars)[i]]]==j)))
      res[whichi, c("n.ref")] <- sapply(vars[[i]][1],function(j)length(which(data[[names(vars)[i]]]==j)))
    }
    if(length(ninter)>0){
      for(whichi in ninter){
        nn <- rownames(res)[whichi]
        which2 <-nn
    	#which2 <- c()
        #for(i in names(vars)){
        #  if(length(grep(i, nn))>0)
        #    which2 <- c(which2,i)
        #}
        rownames(res)[whichi] <- paste("Interaction (", paste(which2, collapse=":"),")",sep="")
      }
    }

  }

  res
}
