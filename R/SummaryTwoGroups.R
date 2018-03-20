#' Summary statistics of two-group comparison
#'
#' This function returns a summary of two group comparison in terms of effect size, lower and upper limits of CI, and p-value.
#' Three classes of "outcome.var" variable can be analyzed using this function:
#' 1) Cox proportional hazards model for survival outcome.var using coxph()
#' 2) T-test for continuous outcome.var using lm()
#' 3) Z-test for binary outcome.var using prop.test().
#'
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param outcome.var a vector specifying the outcome variable. For 'binary' outcome.var, it should be a vector of 1 or 0. In case of a 'survival' variable, this will be a matrix of two columns: 1) time to event 2) censorship.
#' @param subgroup.var a vector of row index specifying the subgroup to be included for the analysis. If NULL (default), all data will be used.
#' @param surv.conf.type confidence interval type. See conf.type in survfit. Default is "plain"
#' @param ties Default is "efron". To match internal sas results, use "exact". See parameter "ties" in coxph.
#' @param treatment.var the name of the treatment variable.
#' @param placebo.code the name of the control group within the treatment variable.
#' @param active.code the name of the treatment/experimental group within the treatment variable.
#' @param outcome.class the outcome class of the 'outcome.var' variable. One of the 3 values - "survival", "binary", or "continuous".
#' @param alpha the confidence level (CI) for point estimate, i.e. 0.05 (default) for 95 percent CI.
#' @param covariate.var a vector specifying the covariate variables. This can be added to adjust for in the analysis for survival and continuous outcome.var variable classes. Default is NULL.
#' @param strat.factor.var a vector specifying the stratification variables. This can be added for the survival outcome.var variable class. Default is NULL.
#' @param return.fit if TRUE, returns a table of summary statistics. Default is FALSE.
#' @param fit.para a list of fitting parameters. Currently only \code{'prop.test.use.continuity.correction'} in use.
#' If \code{'prop.test.use.continuity.correction' = T} (default), the 'correct' parameter in \code{\link{prop.test}} will be set as TRUE.
#'
#' @return A named vector of following entries:
#' if binary - Effect.Size (Proportion Difference), Lower, Upper, P, Rsp.Placebo, Rsp.Active,  N.Placebo, N.Active, nRsp.Placebo, nRsp.Active;
#' if survival - [Events, N, Median Suvival Time] for each group, Effect.Size (Hazard Ratio), Lower, Upper, Wald P;
#' if continuous - Effect.Size (Mean Difference), Lower, Upper, P.
#'
#' @note This function requires "survival" package to call the coxph() function. Two treatment arms are required.
#' Treatment group variable can be forced into a factor. Censorship variable is 1 if an event happened, 0 if censored.
#'
#' @importFrom survival coxph Surv strata
#'
#' @examples
#' data(input)
#' SummaryTwoGroups(outcome.var = input$OS, treatment.var = input$Arm, placebo.code = "CTRL", active.code = "TRT", outcome.class = "continuous",surv.conf.type="plain")
#'
#' @export

SummaryTwoGroups <- function(outcome.var,
                        subgroup.var=NULL,
                        treatment.var,
                        placebo.code,
                        active.code,
                        outcome.class,
                        alpha=0.05,
			surv.conf.type="plain", ties="efron",
                        covariate.var=NULL,
                        strat.factor.var=NULL,
                        return.fit=FALSE,
                        fit.para = list('prop.test.use.continuity.correction'=T)) {


    # If subgroup.var is not defined, use all input data
    if (is.null(subgroup.var)) {
        subgroup.var <- 1:length(treatment.var)
    }

    # Make the treatment assingment a factor, with ordered levels
    treatment.var <- factor(treatment.var, levels=c(placebo.code, active.code))
    # Binary outcome - e.g., response
    if (outcome.class == "binary" ){

        y1 <- outcome.var[subgroup.var][treatment.var[subgroup.var] == active.code]
        n1 <- length(na.omit(y1))
        y2 <- outcome.var[subgroup.var][treatment.var[subgroup.var] == placebo.code]
        n2 <- length(na.omit(y2))

        r1 <- sum(y1, na.rm = TRUE)
        r2 <- sum(y2, na.rm = TRUE)

        mytest <- prop.test(c(r1, r2), c(n1, n2), conf.level = 1 - alpha,correct =  fit.para[['prop.test.use.continuity.correction']])

        ret <- c("Effect.Size" = r1/n1 - r2/n2
                 , "Lower" = mytest$conf.int[1]
                 , "Upper" = mytest$conf.int[2]
                 , "P" = mytest$p.value
                 , "Rsp.Placebo" = mytest$estimate[2]
                 , "Rsp.Active" = mytest$estimate[1]
                 , "N.Placebo" = n2
                 , "N.Active" = n1
                 , "nRsp.Placebo" = r2
                 , "nRsp.Active" = r1)
        names(ret)[5:6] <- c("Rsp.Placebo", "Rsp.Active")
    } # end binary

    # Continuous outcome.var - e.g., blood cholesterol level
    else if (outcome.class=="continuous"){

        # No covariate to adjust for
        if (missing(covariate.var) | is.null(covariate.var)){
            myfit <- lm(outcome.var[subgroup.var] ~ treatment.var[subgroup.var])
            # the last row (2nd) is the 'slope' estimate and its associated quantities
            # which corresponds to the treatment effect size
            coef.1 <- summary(myfit)$coef
            mytest <- coef.1[nrow(coef.1), ]
            mytest2 <- coef.1[1, ]
            myCI <- confint(myfit, level = 1-alpha)[nrow(coef.1), ]
        }

        # Covariates to adjust for
        else {

            # One variable?
            if (is.vector(covariate.var) | is.factor(covariate.var)) {
                myfit <- lm(outcome.var[subgroup.var] ~ covariate.var[subgroup.var] + treatment.var[subgroup.var])
                # the last row (3rd) is the 'slope' estimate and its associated quantities
                # which corresponds to the treatment effect size
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), ]
                mytest2 <- coef.1[1, ]
                myCI <- confint(myfit, level = 1-alpha)[nrow(coef.1), ]
            }

            # More than one variable?
            else {
                nCV <- ncol(covariate.var)
                mycov <- paste("covariate.var[subgroup.var,", 1:nCV, "]", collapse = " + ")
                myformula <- paste("lm(outcome.var[subgroup.var] ~"
                                   , mycov, "+ treatment.var[subgroup.var])")
                myfit <- eval(parse(text = myformula))
                # (nCV+2)th row (last) is the 'slope' estimate and its associated quantities
                # which corresponds to the treatment effect size
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), ]
                mytest2 <- coef.1[1, ]
                myCI <- confint(myfit, level = 1-alpha)[nrow(coef.1), ]
            }
        }

        ret <- c(mytest[1], myCI, mytest[4],mytest2[1], mytest2[1]+mytest[1])
        names(ret) <- c("Effect.Size","Lower","Upper","P","Mean.Placebo","Mean.Active")
    } # end continuous

    # Survival outcome.var - e.g., progression-free survival
    else if (outcome.class == "survival") {
        #require(survival)
        Response <- outcome.var[, 1]
        Event <- outcome.var[, 2]
        ####### Unstratified analysis, not adjusting for any covariate(s)
        if ((missing(covariate.var)||is.null(covariate.var)) & (missing(strat.factor.var)||is.null(strat.factor.var))){
            myfit <- coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~ treatment.var[subgroup.var],ties=ties)
            coef.1 <- summary(myfit)$coef[1:5]
            mytest <- coef.1
        }
        ####### Stratified analysis, not adjusting for any covariate(s)
        else if ((missing(covariate.var) || is.null(covariate.var)) & (!missing(strat.factor.var)&&!is.null(strat.factor.var))){
            # One stratification variable
            if (is.vector(strat.factor.var) | is.factor(strat.factor.var)){
                myfit <- coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~
                                   treatment.var[subgroup.var] + strata(strat.factor.var[subgroup.var]),ties=ties)
                coef.1 <- summary(myfit)$coef[1:5]
                mytest <- coef.1
            }
            # More than one stratification variable
            else {
                nSF <- ncol(strat.factor.var)
                mysf <- paste("strat.factor.var[subgroup.var,", 1:nSF, "]", collapse = ", ")
                myformula <- paste("coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~ treatment.var[subgroup.var] + strata("
                                   , mysf, "))")
                myfit <- eval(parse(text = myformula))
                coef.1 <- summary(myfit)$coef[1:5]
                mytest <- coef.1
            }
        }

        ####### Unstratified analysis, adjusting for covariate(s)
        else if ((!missing(covariate.var)&!is.null(covariate.var)) & (missing(strat.factor.var)||is.null(strat.factor.var))){
            # One covariate
            if (is.vector(covariate.var) | is.factor(covariate.var)){
                myfit <- coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~
                                   covariate.var[subgroup.var] + treatment.var[subgroup.var],ties=ties)
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), 1:5]
            }
            # More than one covariate
            else {
                nCV <- ncol(covariate.var)
                mycov <- paste("covariate.var[subgroup.var,", 1:nCV, "]", collapse = " + ")
                myformula <- paste("coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~"
                                   , mycov, "+ treatment.var[subgroup.var])")
                myfit <- eval(parse(text = myformula))
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), 1:5]
            }

        }
        ####### Stratified analysis adjusting for covariate(s)
        else {
            # One covariate and one stratification variable
            if ((is.vector(covariate.var) | is.factor(covariate.var)) &
                (is.vector(strat.factor.var) | is.factor(strat.factor.var))){
                myfit <- coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~
                                   covariate.var[subgroup.var] + treatment.var[subgroup.var] +
                                   strata(strat.factor.var[subgroup.var]), ties=ties)
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), 1:5]
            }
            # Multiple covariates and one stratification variable
            else if (is.vector(strat.factor.var) | is.factor(strat.factor.var)){
                nCV <- ncol(covariate.var)
                mycov <- paste("covariate.var[subgroup.var,", 1:nCV, "]", collapse = " + ")
                myformula <- paste("coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~"
                                   , mycov, "+ treatment.var[subgroup.var] +
                                   strata(strat.factor.var[subgroup.var]))")
                myfit <- eval(parse(text = myformula))
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), 1:5]
            }
            # One covariate and multiple stratification variables
            else if (is.vector(covariate.var) | is.factor(covariate.var)){
                nSF <- ncol(strat.factor.var)
                mysf <- paste("strat.factor.var[subgroup.var,", 1:nSF, "]", collapse = ", ")
                myformula <- paste("coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~
                                   covariate.var[subgroup.var] + treatment.var[subgroup.var] + strata("
                                   , mysf, "))")
                myfit <- eval(parse(text = myformula))
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), 1:5]
            }
            # Multiple covariates and multiple stratification variables
            else {
                nCV <- ncol(covariate.var)
                mycov <- paste("covariate.var[subgroup.var,", 1:nCV, "]", collapse = " + ")
                nSF <- ncol(strat.factor.var)
                mysf <- paste("strat.factor.var[subgroup.var,", 1:nSF, "]", collapse = ", ")

                myformula <- paste("coxph(Surv(Response[subgroup.var], Event[subgroup.var]) ~"
                                   , mycov, "+ treatment.var[subgroup.var] + strata("
                                   , mysf, "))")
                myfit <- eval(parse(text = myformula))
                coef.1 <- summary(myfit)$coef
                mytest <- coef.1[nrow(coef.1), 1:5]
            }
        }

        Effect.Size <- mytest[2]
        Lower <- exp(mytest[1] - qnorm(1 - alpha/2) * mytest[3])
        Upper <- exp(mytest[1] + qnorm(1 - alpha/2) * mytest[3])

        ret <- c(as.numeric(t(summary(survfit(Surv(Response[subgroup.var], Event[subgroup.var]) ~ treatment.var[subgroup.var],
					      conf.type=surv.conf.type))$table[,c("events","n.start","median"),drop=FALSE])),
                 Effect.Size, Lower, Upper, mytest[5])

        names(ret) <- c(paste(rep(levels(treatment.var),each=3), rep(c("events","n","MST"),2),sep=".")
                        ,"Effect.Size","Lower","Upper","P")
    } # end survival

    else {
        stop("Please input outcome.class as one of these three options: binary, continuous or survival...\n")
    }

    if(!return.fit) {
        return(ret)
    } else {
        return(list(ret, myfit))
    }
}
