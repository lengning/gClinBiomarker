#' Summary a single covariate, test across treatment and/or population 
#' 
#' This function performs demographics imbalance checking of a single covariate across multiple groups.
#' 
#' @author Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
#' 
#' @param data Input data frame. Rows are patients and columns are variables (e.g. demographics variables, time to event variables, 
#' biomarker variables, treatment indicator, etc.). One patient per row. 
#' @param var name of the clinical covariate to test. name should be in the column names of data. entries with empty value (nchar()==0) will be imputed as NA
#' @param trt name of the treatment column. If trt is specified, the analysis will be performed within treatment arm.
#' if it is NULL, the comparison will be performed using all samples.
#' @param trt.name preferred display name of the treatment variable
#' If it is NULL, trt will be used.
#' @param subgroup name of the column which indicates subpopulation (e.g. biomarker evaluable population)
#' @param subgroup.name preferred display name of the subpopulation (e.g. biomarker evaluable population).
#' If it is NULL, subgroup will be used.
#' @param itt.name preferred display name of the full population (e.g. ITT).
#' If it is NULL, "All" will be used.
#' @param subgroup.indicator In the subpopulation column, which value is used
#' to define the subpopulation (e.g. biomarker evaluable population). Default is 1. It can also be character or logical.
#' Default is 1. The non-subpopulation enrties is not allowed to be specified as NA.
#' @param var.class class of the variable. possible categories are "numeric", "categorical" and
#' "ordered.factor".  "ordered.factor" can be used to categorical variable with
#' ordered levels - e.g. IC score 0/1/2/3. If class is ordered.factor ,
#' ordered.factor.levels need to be specified.
#' If it is not specified, will try to use the class of the column.
#' "numeric","integer" will be treated as "numeric"
#' "logical""character","factor" will be treated as "categorical".
#' @param ordered.factor.levels ordered levels of the ordered factor. 
#' @param cont.show what summary statistics to show for a continuous covariate.
#' Default is c("N" ,"Mean","Median", "Min-Max","NA's").
#' Possible options are "N" ,"Mean","SEM", "SD","Median",
#' "Min","Max" ,"Min-Max","1st Qrtl.","3rd Qrtl.","IQR" ,"NA's"
#' @param digits digits for rounding
#' @param trt.order If the user wants to display the treatments in a certain
#' order, it can be defined here. All elements in trt.order should be the same
#' unique values in the treatment column.
#' @param na.action defaults to "na.omit". Possible options are "na.omit", "error"
#' When it is specified as "na.omit", entries with missing trt or subgroup
#' will be automatically removed before calculation.
#' @param compare.subgroup If it is TRUE,
#' the output will show summary statistics of subgroup and others. Default is FALSE. If it is FALSE,
#' will show summary statistics of subgroup vs. All patients
#' @param test.subgroup whether test across subpopulations within treatment arm. If class is numeric,
#' kruskal wallis rank sum test will be performed. If class is categorical, fisher's exact test will be performed.
#' If class is ordered.factor, cmh test will be performed. The test is always performed between subgroup vs others.
#' P value columns will be included in the output table if it is specified as TRUE.
#' Testing is not recommendated if either subgroup of non-subgroup has small sample size.
#' 
#' @note This function provides summary statistics of a single clinical covariate. Using default parameters,
#' the function provides a table to compare summary statistics in All patients vs. in BEP (biomarker evaluable population),
#' within treatment arm
#' @note trt allows for more than 2 levels. However, only 2 levels are allowed for subgroup.
#' For more general use, a user can specify trt to get summary statistics for any
#' sub-group defination (and leave subgroup as NULL).
#' 
#' @importFrom stats as.formula complete.cases fisher.test kruskal.test sd
#' @importFrom coin cmh_test pvalue
#' 
#' @export

SummarySingle <- function (data, var, 
			trt = NULL, trt.name = NULL, 
      subgroup = NULL, subgroup.name = NULL, subgroup.indicator=1, compare.subgroup=FALSE,itt.name="All",
			var.class=NULL, ordered.factor.levels=NULL,
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
  data[[var]][which(nchar(as.character(data[[var]]))==0)] <- NA
  possible.class <-c("categorical","numeric","ordered.factor")
  if(is.null(var.class)||!all(var.class%in%possible.class)){
  if(class(data[,var])%in%c("numeric","integer"))var.class <- "numeric"
  if(class(data[,var])%in%c("logical"))class(data[,var]) <- "character"
  if(class(data[,var])%in%c("character","factor"))var.class <- "categorical"
  }
  if(is.null(var.class)||!all(var.class%in%possible.class))stop(paste('var.class should be in', paste(possible.class,collapse=",")))

  if(var.class=="ordered.factor"){
    if(is.null(ordered.factor.levels)) stop("If class is ordered.factor,
		ordered.factor.levels need to be specified.")
    if(!is.null(ordered.factor.levels)){
      ordered.factor.levels<- as.character(ordered.factor.levels)
      if(!identical(sort(ordered.factor.levels),sort(levels(factor(data[,var])))))stop("ordered factor levels should match unique elements in the variable!")
      data[,var] <- factor(data[,var],levels=ordered.factor.levels, ordered=TRUE)
      Lev <- ordered.factor.levels
    }
  }
  
  if( (!is.null(subgroup))) subgroup.name <- ifelse(is.null(subgroup.name), subgroup, subgroup.name)
  if( (!is.null(trt))) trt.name <- ifelse(is.null(trt.name), trt, trt.name)
  

# TRT missing value
  if (!is.null(trt)) {
    if (any(is.na(data[, trt]))) {
      if (na.action == "na.omit") {
        data <- data[complete.cases(data[, trt]), ]
        warning("There were observations without 'trt' information! These were omitted (na.action=\"na.omit\")!")
      }
      else stop("There were observations without 'trt' information!")
    }
  }

# subgroup missing value
  if (!is.null(subgroup)) {
    if (length(subgroup) == 1 && subgroup %in% colnames(data)) {
      if (any(is.na(data[, subgroup]))) {
        if (na.action == "na.omit") {
          data <- data[complete.cases(data[, subgroup]), ]
          warning("There were observations without 'subgroup' information! These were omitted (na.action=\"na.omit\")!")
        }
        else stop("There were observations without 'subgroup' information!")
      }
    }
  }

 


  # if trt is NULL, create a column with all 1s as indicator
  if (is.null(trt)) {
    data$trt <- rep(1, nrow(data))
    trt.lev <- unique(data$trt)
  } else {
    if (!is.null(trt.order)) {
      tmp <- levels(factor(data[, trt]))
      if(!identical(sort(trt.order),sort(tmp))) stop("trt.order should match levels in trt column!")
      data$trt <- factor(data[, trt], levels = trt.order)
    }
    else data$trt <- factor(data[, trt])
    trt.lev <- levels(data$trt)
  }

  # if subgroup is NULL, create a column with all 1s as indicator
   data$ITT <- rep(1, nrow(data))
   compare.itt <- !compare.subgroup 
    if (!is.null(subgroup)) {
      if(compare.itt){
        data[,subgroup] <- ifelse(data[,subgroup]%in%subgroup.indicator,1,0)
        subgroup.name <- c(itt.name,subgroup.name)
        subgroup <- c("ITT",subgroup)
      }
      if(!compare.itt){
        subgroup.l <- paste0(subgroup.name,"_", c(subgroup.indicator,paste0("not_",subgroup.indicator)))
        
        tmp <- unique(data[which(!data[,subgroup]%in%subgroup.indicator),subgroup])
        if(length(tmp==1)) subgroup.l[2] <- paste0(subgroup.name,"_",  tmp)
        data[,subgroup.l[1]] <- ifelse(data[,subgroup]%in%subgroup.indicator,1,0)
        data[,subgroup.l[2]] <- ifelse(data[,subgroup]%in%subgroup.indicator,0,1)
        
        subgroup <- subgroup.l
        subgroup.name <- subgroup.l
      }
  }
  if(is.null(subgroup)) {
	  subgroup <- "ITT"
	  subgroup.name <-  " "
	  }
  # output matrix (trt by population)
  res <- matrix(ncol = sum(1, length(trt.lev) * length(subgroup)), 
                nrow = 0)

  if (var.class== "categorical") {
    data[,var] <- factor(data[,var])
    Lev <- levels(data[, var])
  }
  result <- vector("list", length(trt.lev))
  names(result) <- rep("",length(trt.lev))
  
  if (!is.null(trt.name)) 
      names(result) <- trt.lev
  

  res.ind <- 1
  for (i in 1:length(trt.lev)) {
    trt.mat <- data[which(data$trt == trt.lev[i]), ]
    for (j in 1:length(subgroup.name)) {
      subgroup.mat <- trt.mat[which(trt.mat[, subgroup[j]]==1), ]
      Nna <- length(which(is.na(subgroup.mat[, var])))
      Na <- nrow(subgroup.mat) - Nna

      if (var.class%in%c("categorical","ordered.factor") ){
        res.subgroup <- data.frame(1)
          res.subgroup$"Total (non-NA)" <- Na
          res.subgroup$"NA's" <- Nna
          for (k in Lev) {
            res.subgroup[, k] <- length(which(subgroup.mat[, var] ==  k))
            res.subgroup[, k] <- paste(res.subgroup[, k], " ", 
                                    "(", 100 * round(res.subgroup[, k]/Na, digits + 
                                                       2), "%)", sep = "")
          }
        #}
        res.subgroup <- t(res.subgroup[, -1])
      } else {
        smry <- summary(subgroup.mat[, var], digits = 16)
        res.subgroup <- data.frame(1)
          res.subgroup$N <- Na
          res.subgroup$Mean <- round(ifelse(is.nan(smry["Mean"]), 
                                              NA, smry["Mean"]),  digits=digits)
          res.subgroup$SEM <- round(sd(subgroup.mat[, var], 
                                         na.rm = TRUE)/sqrt(length(which(!is.na(subgroup.mat[, 
                                                            var])))),  digits=digits)
          res.subgroup$SD <- round(sd(subgroup.mat[, var], 
                                        na.rm = TRUE),  digits=digits)
          res.subgroup$Median <- round(smry["Median"], 
                                         digits=digits)
          res.subgroup$Min <- round(smry["Min."],  digits=digits)
          res.subgroup$Max <- round(smry["Max."],  digits=digits)
          if (any(is.na(smry[c("Min.", "Max.")]))) 
            res.subgroup$"Min-Max" <- NA
          else res.subgroup$"Min-Max" <- paste(round(smry["Min."], 
                                                        digits=digits), 
					round(smry["Max."],  digits=digits), 
                                          sep = "...")
        
          res.subgroup$"1st Qrtl." <- round(smry["1st Qu."],  digits=digits)
          res.subgroup$"3rd Qrtl." <- round(smry["3rd Qu."],   digits=digits)
          res.subgroup$IQR <- round(smry["3rd Qu."] -   smry["1st Qu."],  digits=digits)
          res.subgroup$"NA's" <- Nna
          res.subgroup <- t(res.subgroup[cont.show])

      }
      if (j == 1) 
        p.res <- res.subgroup
      else p.res <- cbind(p.res, res.subgroup)
    }
    colnames(p.res) <- subgroup.name
    # If all ITT are in subgroup, no test can be performed
    if(test.subgroup){
      if(length(which(trt.mat[,subgroup[2]]!=1))==0){
        test.subgroup  <-  FALSE
        message("All ITT patients are in subgroup. subgroup.test is set to FALSE")
      }  
    }
    if(test.subgroup){
      if(var.class=="numeric")
        tt <- kruskal.test(list(trt.mat[which(trt.mat[,subgroup[2]]!=1),var],
                     trt.mat[which(trt.mat[,subgroup[2]]==1),var]))$p.value
    if(var.class=="categorical")
      tt <- fisher.test(table(trt.mat[,c(var,subgroup[2])]))$p.value
    if(var.class=="ordered.factor"){
      #require(coin)
      tmp.mat <- trt.mat
      tmp.mat[[subgroup[2]]] <- as.factor(tmp.mat[[subgroup[2]]])
      tt <- pvalue(cmh_test(as.formula(paste0(var,' ~ ',subgroup[2])),data=tmp.mat))
    }
    p.res <- cbind(p.res, pvalue=c(round(tt,digits),rep("",nrow(p.res)-1)))
  }
    if(length(result)>1)colnames(p.res) <- paste0(colnames(p.res),"(",names(result)[i],")")  
    result[[i]] <- p.res
  }
  
  colns <- unlist(sapply(result,colnames,simplify=F))
  rowns <- unique(unlist(sapply(result,rownames,simplify=F)))
  outmat <- matrix("",nrow=length(rowns),ncol=length(colns),dimnames=list(rowns, colns))
  for(i in 1:length(result)) outmat[rownames(result[[i]]),colnames(result[[i]])] <- result[[i]]
  outmat
}
