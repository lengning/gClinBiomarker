#' Summary multiple covariates, test across treatment and/or population 
#' 
#' This function performs demographics imbalance checking of a single covariate across multiple groups.
#' 
#' @author Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
#' 
#' @param var a vector of covariate names - the clinical covariate to test
#' @param var.name preferred display names of the clinical covariates 
#' If it is NULL, var will be used.
#' @param ordered.factor.levels.list a list indicates ordered levels for ordered.factor. Each ordered.factor
#' should have a corresponding element in this list. 
#' @param var.class a vector that indicates class of the variables. possible categories are "numeric", "categorical" and
#' "ordered.factor".  "ordered.factor" can be used to categorical variable with
#' ordered levels - e.g. IC score 0/1/2/3. If class is ordered.factor ,
#' ordered.factor.levels need to be specified.
#' If the user doesn't specify class of all variables (the length of the var.class is less than length of var), 
#' The program will try to use the class of the column.
#' "numeric","integer" will be treated as "numeric"
#' "logical""character","factor" will be treated as "categorical".
#' In this case the program request that names of the vector var.class is a subset of the var vector. 
#' 
#' @return output object is a matrix with summary statistics. It can be passed to knitr::kable(). 
#' 
#' @note trt allows for more than 2 levels. However, only 2 levels are allowed for subgroup.
#' For more general use, a user can specify trt to get summary statistics for any
#' sub-group defination (and leave subgroup as NULL).
#' @note This function provides summary statistics of a vector of clinical covariates. Using default parameters,
#' the function provides a table to compare summary statistics in All population vs. in BEP (biomarker evaluable population),
#' within treatment arm.
#' 
#' @inheritParams SummarySingle
#' 
#' @export

SummaryVars <- function (data, var, var.name = NULL, 
			trt = NULL, trt.name = NULL, 
      subgroup = NULL, subgroup.name = NULL, subgroup.indicator=1, compare.subgroup=FALSE,itt.name="All",
			var.class=NULL, ordered.factor.levels.list=NULL,
			cont.show = c("N" ,"Mean","Median", "Min-Max","NA's"),
			digits = 2, trt.order = NULL, test.subgroup=FALSE, 
				 na.action = "error") 
{
  stopifnot(na.action%in% c("na.omit", "error"))
  stopifnot(class(data) == "data.frame")

  if(!all(c(var, trt, subgroup) %in% colnames(data)))stop("var, trt and subgroup should have matched column names in the input data!")
  if(!is.null(subgroup)) if(nlevels(as.factor(data[,subgroup]))<2)stop("subpopulation column has only one unique value!")
  possible.show <- c("N" ,"Mean","SEM", "SD","Median",
                     "Min","Max" ,"Min-Max","1st Qrtl.","3rd Qrtl.",
                     "IQR" ,"NA's")
  if(length(setdiff(cont.show, possible.show))>0) 
    stop(paste("possible cont.show elements are:", possible.show ))
  if(test.subgroup & is.null(subgroup)){
    test.subgroup <- FALSE
    message("test.subgroup=TRUE but subgroup is not specified. Reset test.subgroup as FALSE")
  }

  if(test.subgroup & compare.subgroup==FALSE){
       compare.subgroup <- T
       message("test.subgroup is TRUE but compare.subgroup is FALSE. Set compare.subgroup to TRUE")
        }

  possible.class <-c("categorical","numeric","ordered.factor")
  if(!all(var.class%in%possible.class))stop(paste('var.class should be in', paste(possible.class,collapse=",")))

  if(length(var)== length(var.class)) names(var.class) <- var
  if(length(var)!= length(var.class)) if(!all(names(var.class)%in%var))stop("length of var.class doesn't match length of var - 
  In this case the program request that names of the vector var.class is a subset of the var vector")

  
  k.var <- length(var)
  for(i in 1:k.var){
   tmp <- try(SummarySingle(data=data, var=var[i], 
			trt = trt, trt.name = trt.name, 
                        subgroup = subgroup, subgroup.name = subgroup.name, subgroup.indicator=subgroup.indicator, 
			compare.subgroup=compare.subgroup, itt.name=itt.name,
			var.class=var.class[i], ordered.factor.levels=ordered.factor.levels.list[[var[i]]],
			cont.show = cont.show, 
			digits = digits, trt.order = trt.order, test.subgroup=test.subgroup, 
				 na.action = na.action))
   if(class(tmp)=="try-error") stop(paste("error in", var[i] ))
   if(i==1)kcol <- ncol(tmp)
   line1 <- rep("", kcol) 
   tmp2 <- rbind(line1, tmp)
   rownames(tmp2)[1] <- ifelse(is.null(var.name[i]),var[i], var.name[i]) 
   if(i==1)mat <- tmp2 
   if(i>1)mat <- rbind(mat, tmp2)
  }
 
 mat 

}
