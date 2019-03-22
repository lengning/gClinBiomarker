#' Create a STEPP (Subpopulation Treatment Effect Pattern Plot)
#'
#' This function creates a STEPP from the given point estimates and confidence intervals at desired percentiles.
#'
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param data input data frame. Rows are patients and columns are variables (e.g. demographics variables, time to event variables,
#' biomarker variables, treatment indicator, etc.). One patient per row.
#' @param outcome.var outcome variable. In case of a 'survival', variable, it will be a vector of two variables: 1) time to event 2) censorship
#' @param outcome.class outcome class of the 'outcome' variable. Can be either "continuous", "binary", or "survival".
#' @param trt name of the treatment variable.
#' @param var name of the biomarker variable.
#' @param covariate vector specifying the covariate variables.
#' This can be added to adjust for in the analysis for survival and continuous outcome variable classes.
#' @param strata vector specifying the stratification variables. This can be added for the survival outcome variable class.
#' @param placebo.code name of the control group within the treatment variable
#' @param active.code name of the treatment/experimental group within the treatment variable
#' @param quantile.type an integer between 1 and 9 selecting one of the nine quantile algorithms. See \code{\link{quantile}}. Default is 2.
#' @param alpha confidence level (CI) for point estimate, i.e. 0.05 for 95 percent CI. Default is 0.05.
#' @param min.pt minimum center.pt. Default is NULL.
#' @param max.pt maximum center.pt. Default is NULL.
#' @param window.width width of each window. Default is 0.25.
#' @param by size of 'slide', i.e. speed of the window that moves along the x-axis. Default is 0.05.
#' @param yrange.lower value of the lower y-axis. Default is NULL.
#' @param yrange.upper value of the upper y-axis. Default is NULL.
#' @param xticks x-tick marks. Default is NULL.
#' @param show.refline if TRUE, the reference line will be displayed. Default is TRUE.
#' @param refline.color color of the reference line. Default is "grey".
#' @param show.refline.ac if TRUE, the reference line for effect in All Comers will be displayed. Default is TRUE.
#' @param refline.color.ac color of the reference line for effect in All Comers. Default is "lightblue".
#' @param estimate.color color of the estimate line. Default is "blue".
#' @param estimate.lty type of the estimate line. Default is 1.
#' @param estimate.lwd width of the estimate line. Default is 2.
#' @param surv.conf.type confidence interval type. Default is "plain". see conf.type in survfit.
#' @param ci.color color of the CI lines. Default is "black".
#' @param ci.lty type of the CI lines. Default is 2.
#' @param ci.lwd width of the CI lines. Default is 1.
#' @param ci.shade if TRUE, the area between the lower and upper CI lines will be shaded in gray. Default is TRUE.
#' @param plot.title title of the plot. Default is "STEPP: Subgroup Treatment Effect Pattern Plot".
#' @param sub.title subtitle of the plot, this will be displayed on the bottom of the plot. Default is NULL.
#' @param xlabel x-axis label. Default is "Biomarker Percentile".
#' @param ylabel y-axis label. Default is NULL.
#' @param show.legend if TRUE, a legend will be displayed. Default is TRUE.
#' @param legend.loc location of the legend to be displayed. Default is "topright".
#' @param legend.text text to be displayed in the legend.
#' Default is \code{c("Point Estimate",paste((1-alpha)*100," percent Confidence Interval",sep=""))}
#' @param legend.col colors of the legend. Default is \code{c(estimate.color, ci.color)}.
#' @param legend.lty type of the legend lines. Default is \code{c(estimate.lty, ci.lty)}.
#' @param legend.lwd width of the legend lines. Default is \code{c(estimate.lwd, ci.lwd)}.
#' @param legend.bty type of box for the legend. Default is "n".
#' @param bm.digits digits to be displayed/used for the lower and upper confidence level estimates. Default is 2.
#' @param actual.scale if TRUE, it generates the figure using the actual scale instead of the scale in percentage. Default is FALSE.
#' @param equal.in.LL,equal.in.UL if both are TRUE, window is defined using >= and <= (in legend: "\code{[ ]}").
#' if (TRUE, FALSE), window is defined using >= and < (in legend: "\code{[ )}").
#' if (FALSE, TRUE), window is defined using > and <= (in legend: "\code{( ]}").
#' if both are FALSE, window is defined using > and < (in legend: "\code{()}"). Default is (TRUE, FALSE).
#' @param csv.name csv file name (includes numbers used in the graphs). If NULL (default), the function will return the result.
#' @param pdf.name name of output pdf file. If it's NULL (default), the plots will be displayed but not saved as pdf.
#' @param pdf.param a list of parameters that define pdf graphics device. See \code{\link{pdf}}. Default is \code{list(width=11, height=8.5)}.
#' @param par.param a list of parameters that define graphcial parameters. See \code{\link{par}}. Default is \code{list(mar=c(4,4,3,2))}.
#' @param ties Default is "efron". To match internal sas results, use "exact". See parameter "ties" in coxph.
#' @note Patient data without corresponding biomarker data are automatically removed.
#' For survival data, censorship variable is 1 if an event happened, 0 if censored.
#'
#' @importFrom survival coxph Surv strata
#'
#' @examples
#'
#' data(input)
#' PlotSTEPP(data = input,
#'          outcome.var = c("PFS", "PFS.event"),
#'          outcome.class = "survival",
#'          trt = "Arm",
#'          var = "KRAS.exprs",
#'          covariate = "Sex",
#'          strata = "Age",
#'          placebo.code = "CTRL",
#'          active.code = "TRT")
#'
#' @seealso
#' \link{SummaryTwoGroups}
#'
#' @export

PlotSTEPP <- function(data,
                     outcome.var,
                     outcome.class,
                     trt=NULL,
                     var,
                     covariate = NULL,
                     strata = NULL,
                     placebo.code = NULL,
                     active.code = NULL,
                     quantile.type = 1,
                     alpha = .05,
                     window.width = .25,
                     min.pt = NULL,
                     max.pt = NULL,
                     by = .05,
                     yrange.lower=NULL,
                     yrange.upper=NULL,
                     xticks = NULL,
                     show.refline = TRUE,
                     refline.color = "grey",
                     show.refline.ac = TRUE,
                     refline.color.ac = "lightblue",
                     estimate.color="blue",
                     estimate.lty = 1,
                     estimate.lwd = 2,
		     surv.conf.type="plain",
		     ties="efron",
                     ci.color = "black",
                     ci.lty = 2,
                     ci.lwd = 1,
                     ci.shade = TRUE,
                     plot.title = "STEPP: Subgroup Treatment Effect Pattern Plot",
                     sub.title = NULL,
                     xlabel = "Biomarker Percentile",
                     ylabel = NULL,
                     show.legend = TRUE,
                     legend.loc = "topright",
                     legend.text = c("Point Estimate", paste((1-alpha)*100," percent Confidence Interval",sep="")),
                     legend.col = c(estimate.color, ci.color),
                     legend.lty = c(estimate.lty, ci.lty),
                     legend.lwd = c(estimate.lwd, ci.lwd),
                     legend.bty = "n",
                     bm.digits = 2,
                     actual.scale = FALSE,
                     equal.in.LL = TRUE,
                     equal.in.UL = FALSE,
                     pdf.name = NULL,
                     pdf.param = list(width= 11, height=8.5),
                     par.param = list(mar=c(4,4,3,2)),
                     csv.name = NULL) {

  stopifnot(class(data)=="data.frame")
  outcome.class <- match.arg(outcome.class, c("survival", "binary", "continuous"))
  stopifnot(all(c(var, outcome.var, trt)%in%colnames(data)))
  if(any(is.na(data[[var]])))message("some NA in var column, will ignore NA entries")
  # Remove patient data without corresponding biomarker value
  data.ori <- data
  data <- data[which(!is.na(data[[var]])),]

    # Read off the data
    Outcome <- data[, outcome.var]
    Biomarker <- data[, var]


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
      if(!nArms%in%c(2))stop("only 2-arm is allowed")
      data[,trt] <- factor(data[,trt],levels=Arms)
      if(is.null(placebo.code))placebo.code <- Arms[1]
      if(is.null(active.code))active.code <- Arms[-1]
    }

    if (!is.null(covariate)) {
        Covariate <- data[, covariate]
    }
    if (!is.null(strata)) {
        Strat.factor <- data[, strata]
    }

    # If min.pt or max.pt are NULL
    if (is.null(min.pt)) {
        min.pt <- window.width/2
    }
    if (is.null(max.pt)) {
        max.pt <- 1 - window.width/2
    }

    nbins <- round((max.pt - min.pt)/by) + 1

    # X-axis ticks (if not specified by user)
    if (is.null(xticks)) {
        if (nbins > 4) {
            allticks <- seq(min.pt, max.pt, by = by)
            if (nbins %% 2 == 0) {
                xticks <- c(min.pt, quantile(allticks, c(.33, .67), type = quantile.type), max.pt)
            } else {
                xticks <- c(min.pt, quantile(allticks, c(.25, .5, .75), type = quantile.type), max.pt)
            }
        } else {
            xticks <- seq(min.pt, max.pt, by = by)
        }
    }

    # Calculate summary statistics for each bin
    cn <- c("Center.pt", "Effect.Size", "Lower", "Upper", "BMV.LL", "BMV.UL", "BMV.Center", "Left.pt", "Right.pt", "N")
    if (outcome.class == "survival") {
        cn <- c(cn, "Events")
    }
    sdata <- matrix(NA, nbins, length(cn))
    colnames(sdata) <- cn
    sdata[, "Center.pt"] <- seq(min.pt, max.pt, by)

    # Effect size in All Comers
    if (is.null(covariate) & is.null(strata)) {
        effect.ac <- SummaryTwoGroups(Outcome, 1:length(Biomarker), Treatment,
                               placebo.code, active.code, outcome.class,
                               alpha, surv.conf.type=surv.conf.type, ties=ties,
                               covariate.var = NULL,
                               strat.factor.var = NULL)["Effect.Size"]
    } else if (!is.null(covariate) & is.null(strata)) {
        effect.ac <- SummaryTwoGroups(Outcome, 1:length(Biomarker), Treatment,
                                 placebo.code, active.code, outcome.class,
                                 alpha,surv.conf.type=surv.conf.type, ties=ties,
                                 covariate.var = Covariate,
                                 strat.factor.var = NULL)["Effect.Size"]
    } else if (is.null(covariate) & !is.null(strata)) {
        effect.ac <- SummaryTwoGroups(Outcome, 1:length(Biomarker), Treatment,
                                 placebo.code, active.code, outcome.class,
                                 alpha,surv.conf.type=surv.conf.type, ties=ties,
                                 covariate.var = NULL,
                                 strat.factor.var = Strat.factor)["Effect.Size"]
    } else if (!is.null(covariate) & !is.null(strata)) {
        effect.ac <- SummaryTwoGroups(Outcome, 1:length(Biomarker), Treatment,
                                 placebo.code, active.code, outcome.class,
                                 alpha,surv.conf.type=surv.conf.type, ties=ties,
                                 covariate.var = Covariate,
                                 strat.factor.var = Strat.factor)["Effect.Size"]
    }
    options(scipen=999)

    for (i in 1:nbins) {
        # Define this window
        START_SEQ <- seq( 0 , 1 - min.pt*2 , by = by)
        start <- START_SEQ[[i]]
        end <- start + window.width

        # Save in sdata
        sdata[i, "Left.pt"] <- start
        sdata[i, "Right.pt"] <- end

        # Biomarker values corresponding to the percentiles
        # 'type'=2 is used to match the result from median()
        LL <- quantile(Biomarker, start, type = quantile.type)
        UL <- quantile(Biomarker, end, type = quantile.type)

        # Round LL and UL
        LL <- round(LL, bm.digits)
        UL <- round(UL, bm.digits)

        if (LL == 0 & UL == 0) {
            stop("No enough unique values, please consider to reduce number of windows (increase window.width) or increase number of digits when rounding (bm.digits)")
        }

        sdata[i, c("BMV.LL", "BMV.UL")] <- c(LL, UL)

        # The subgroup index for this window
        if (equal.in.LL == TRUE & equal.in.UL == FALSE) {
            sindex <- which(Biomarker >= LL  & Biomarker < UL)
        } else if (equal.in.LL == TRUE & equal.in.UL == TRUE) {
            sindex <- which(Biomarker >= LL  & Biomarker <= UL)
        } else if (equal.in.LL == FALSE & equal.in.UL == TRUE) {
            sindex <- which(Biomarker > LL  & Biomarker <= UL)
        } else if (equal.in.LL == FALSE & equal.in.UL == FALSE) {
            sindex <- which(Biomarker > LL  & Biomarker < UL)
        } else {
            stop("equal.in.LL and equal.in.UL can be equal to TRUE or FALSE only!")
        }

        # Calculate median for each bean for Biomarker
        sdata[i, "BMV.Center"] <- median(Biomarker[sindex])

        # Sample size in each bin
        sdata[i, "N"] <- length(sindex)

        if (outcome.class == "survival") {
            sdata[i, "Events"] <- sum(Outcome[sindex, 2], na.rm = TRUE)
        }

        if (is.null(covariate) & is.null(strata)) {
            sdata[i, c("Effect.Size", "Lower", "Upper")] <-
                SummaryTwoGroups(Outcome, sindex, Treatment,
                            placebo.code, active.code,
                            outcome.class, alpha,
                            covariate.var = NULL,
                            strat.factor.var = NULL)[c("Effect.Size", "Lower", "Upper")]
        } else if (!is.null(covariate) & is.null(strata)) {
            sdata[i, c("Effect.Size", "Lower", "Upper")] <-
                SummaryTwoGroups(Outcome, sindex, Treatment,
                            placebo.code, active.code,
                            outcome.class, alpha,surv.conf.type=surv.conf.type, ties=ties,
                            covariate.var = Covariate,
                            strat.factor.var = NULL)[c("Effect.Size", "Lower", "Upper")]
        } else if (is.null(covariate) & !is.null(strata)) {
            sdata[i, c("Effect.Size", "Lower", "Upper")] <-
                SummaryTwoGroups(Outcome, sindex, Treatment,
                            placebo.code, active.code,
                            outcome.class, alpha,surv.conf.type=surv.conf.type, ties=ties,
                            covariate.var = NULL,
                            strat.factor.var = Strat.factor)[c("Effect.Size", "Lower", "Upper")]
        } else if (!is.null(covariate) & !is.null(strata)) {
            sdata[i, c("Effect.Size", "Lower", "Upper")] <-
                SummaryTwoGroups(Outcome, sindex, Treatment,
                            placebo.code, active.code,
                            outcome.class, alpha,surv.conf.type=surv.conf.type, ties=ties,
                            covariate.var = Covariate,
                            strat.factor.var = Strat.factor)[c("Effect.Size", "Lower", "Upper")]
        }
    } # end of for (i in 1:nbins)

    ##### Print out sdata
    sdata2 <- sdata
    colnames(sdata2)[1:9] <- c("Window Center", "Effect.Size", "CI Lower", "CI Upper",
                               "BM Lower", "BM Upper", "BM Center", "Window Left", "Window Right")

    if (!is.null(ylabel)) {
        colnames(sdata2) <- gsub("Effect.Size", ylabel, colnames(sdata2))
    } else {
        if (outcome.class == "survival") {
            colnames(sdata2) <- gsub("Effect.Size", "Hazard Ratio", colnames(sdata2))
        } else if (outcome.class == "binary") {
            colnames(sdata2) <- gsub("Effect.Size", "Proportion Difference", colnames(sdata2))
        } else if (outcome.class == "continuous") {
            colnames(sdata2) <- gsub("Effect.Size", "Mean Difference", colnames(sdata2))
        }
    }

    if (!is.null(csv.name)) {
        write.csv(sdata2, file = csv.name, row.names = FALSE)
    }

    ##### Draw the plot
    PlotParam(pdf.name, pdf.param, par.param)

    if(is.null(sub.title))sub.title <- paste("soc:", placebo.code, "; trt:", active.code)

    center.pt <- sdata[, "Center.pt"]
    effect.size <- sdata[, "Effect.Size"]
    lower = sdata[, "Lower"]
    upper = sdata[, "Upper"]
    bml = sdata[, "BMV.LL"]
    bmu = sdata[, "BMV.UL"]
    if (outcome.class == "survival") {
        num.events = sdata[, "Events"]
    }
    num.pats = sdata[, "N"]

    if (actual.scale == FALSE) {
        # If center.pt is NULL
        if (is.null(center.pt)) {
            center.pt <- seq(min.pt, max.pt, by = by)
        }

        # x values should always be between 0 and 1
        if (sum(center.pt <= 0 | center.pt >= 1) > 0) {
            warning("Center percentile must be between 0 and 1 !")
        }
    }

    # Starting point should be smaller than Ending point
    if (min.pt < 0 | min.pt > max.pt | max.pt > 1 ) {
        warning("Inappropriate center percentile; Please check.")
        print(paste("window.width", window.width))
        print(paste("min.pt", min.pt))
        print(paste("max.pt", max.pt))
    }

    if (is.null(yrange.lower)) {
        yrange.lower= min(lower, na.rm=TRUE)
    }

    if (is.null(yrange.upper)) {
        yrange.upper= max(upper, na.rm=TRUE)
    }

    if (outcome.class == "binary") {
        if(is.null(ylabel)) {
            ylabel <- "Proportion Difference (%)"
        }
        if (yrange.lower < -1 | yrange.upper > 1) {
            warning("Input yrange.lower and yrange.upper should be between -1 and 1 !")
        }
    } else if (outcome.class == "survival") {
        if (is.null(ylabel)) {
            ylabel <- "Hazard Ratio"
        }
        if (yrange.lower < 0 | yrange.upper < 0) {
            warning("Input yrange.lower and yrange.upper should be >0 !")
        }
    } else if (outcome.class == "continuous") {
        if (is.null(ylabel)) {
            ylabel <- "Mean Difference"
        }
    }

    # Set an empty panel
    if (outcome.class == "binary") {
        plot(c(ifelse(actual.scale == FALSE, 0, sdata[1, "BMV.Center"]), ifelse(actual.scale == FALSE, 1, sdata[nbins, "BMV.Center"])),
             c(yrange.lower, yrange.upper),
             type = "n",
             axes = FALSE,
             xlim = c(ifelse(actual.scale == FALSE, min.pt, sdata[1, "BMV.Center"]), ifelse(actual.scale == FALSE, max.pt, sdata[nbins, "BMV.Center"])),
             xlab = xlabel,
             ylab = ylabel,
             main = plot.title,
             sub = sub.title)
        # y-axis in % scale
        axis(2,
             -10:10/10,
             paste(c(rep("", 11), rep("+", 10)), -10:10*10, sep = ""),
             las = 2)
        h_ <- 0
        text_ <- paste(round(effect.ac*100), "%", sep = "")
    } else if (outcome.class == "survival") {
        plot(c(ifelse(actual.scale == FALSE, 0, sdata[1, "BMV.Center"]), ifelse(actual.scale == FALSE, 1, sdata[nbins, "BMV.Center"])),
             c(yrange.lower, yrange.upper),
             type = "n",
             log = "y",
             xaxt = "n",
             xlim = c(ifelse(actual.scale == FALSE, min.pt, sdata[1, "BMV.Center"]), ifelse(actual.scale == FALSE, max.pt, sdata[nbins, "BMV.Center"])),
             xlab = xlabel,
             ylab = ylabel,
             las = 2,
             main = plot.title,
             sub = sub.title)
        h_ <- 1
        text_ <- round(effect.ac, 2)
    } else if (outcome.class == "continuous") {
        plot(c(ifelse(actual.scale == FALSE, 0, sdata[1, "BMV.Center"]), ifelse(actual.scale == FALSE, 1, sdata[nbins, "BMV.Center"])),
             c(yrange.lower, yrange.upper),
             type = "n",
             xaxt = "n",
             xlim = c(ifelse(actual.scale == FALSE, min.pt, sdata[1, "BMV.Center"]), ifelse(actual.scale == FALSE, max.pt, sdata[nbins, "BMV.Center"])),
             xlab = xlabel,
             ylab = ylabel,
             main = plot.title,
             sub = sub.title)
        h_ <- 0
        text_ <- round(effect.ac, 1)
    }

    #title(line = 0, xlab = xlabel)

    # Reference lines for no effect/all comers result
    if (show.refline) {
        abline(h = h_, col = refline.color)
    }
    if (show.refline.ac) {
        abline(h = effect.ac, col = refline.color.ac)
        mtext(text = text_,
              side = 4,
              line = 0,
              at = effect.ac,
              las = 1)
    }

    # X-axis ticks for actual scale
    if (actual.scale == TRUE) {
        xticks2 <- sdata[sdata[, "Center.pt"] %in% xticks, "BMV.Center"]
        old.xticks = xticks
        xticks = xticks2
    }

    # Annotate windows
    abline(v = xticks, col = refline.color, lty = 2)
    if (actual.scale == TRUE) {
        xticks = old.xticks
    }

    # ticks on x-axis: percentiles
    lefts <- paste((xticks - window.width/2)*100, "%", sep = "")
    rights <- paste((xticks + window.width/2)*100, "%", sep = "")

    # ticks on x-axis: biomarker values
    bmlefts <- bml[match(as.numeric(xticks), as.numeric(center.pt))]
    bmrights <- bmu[match(as.numeric(xticks), as.numeric(center.pt))]

    if (actual.scale == TRUE) {
        xticks = xticks2
        #mtext(paste("Median value:", xticks[1], sep = ""), side = 1, line = 0.2, adj=1, at = xticks[1], cex = 0.6)
        #for (i in 2:length(xticks)) {
        #    mtext(xticks[i], side = 1, line = 0.2, at = xticks[i], cex = 0.6)
        #}
        mtext("Median value:            ", side = 1, line = 0.2, adj=1, at = xticks[1], cex = 0.6)
        mtext(xticks, side = 1, line = 0.2, at = xticks, cex = 0.6)
    }

    if (equal.in.LL == TRUE & equal.in.UL == FALSE) {
        for (i in 1:(length(xticks)-1)) {
            axis(1, xticks[i], paste("[", lefts[i], ", ", rights[i], ")", sep = ""), cex.axis = 0.6)
            mtext(paste("[", bmlefts[i], ", ", bmrights[i], ")", sep = ""), side = 1, line = 1.7, at = xticks[i], cex = 0.6)
        }
        axis(1, xticks[length(xticks)], paste("[", lefts[length(xticks)], ", ", rights[length(xticks)], "]", sep = ""), cex.axis = 0.6)
        mtext(paste("[", bmlefts[length(xticks)], ", ", bmrights[length(xticks)], "]", sep = ""), side = 1, line = 1.7, at = xticks[length(xticks)], cex = 0.6)
    } else if (equal.in.LL == TRUE & equal.in.UL == TRUE) {
        axis(1, xticks, paste("[", lefts, ", ", rights, "]", sep = ""), cex.axis = 0.6)
        mtext(paste("[", bmlefts, ", ", bmrights, "]", sep = ""), side = 1, line = 1.7, at = xticks, cex = 0.6)
    } else if (equal.in.LL == FALSE & equal.in.UL == TRUE) {
        axis(1, xticks[1], paste("[", lefts[1], ", ", rights[1], "]", sep = ""), cex.axis = 0.6)
        mtext(paste("[", bmlefts[1], ", ", bmrights[1], "]", sep = ""), side = 1, line = 1.7, at = xticks[1], cex = 0.6)
        for (i in 2:length(xticks)) {
            axis(1, xticks[i], paste("(", lefts[i], ", ", rights[i], "]", sep = ""), cex.axis = 0.6)
            mtext(paste("(", bmlefts[i], ", ", bmrights[i], "]", sep = ""), side = 1, line = 1.7, at = xticks[i], cex = 0.6)
        }
    } else {
        axis(1, xticks[1], paste("[", lefts[1], ", ", rights[1], ")", sep = ""), cex.axis = 0.6)
        mtext(paste("[", bmlefts[1], ", ", bmrights[1], ")", sep = ""), side = 1, line = 1.7, at = xticks[1], cex = 0.6)
        for (i in 2:(length(xticks)-1)) {
            axis(1, xticks[i], paste("(", lefts[i], ", ", rights[i], ")", sep = ""), cex.axis = 0.6)
            mtext(paste("(", bmlefts[i], ", ", bmrights[i], ")", sep = ""), side = 1, line = 1.7, at = xticks[i], cex = 0.6)
        }
        axis(1, xticks[length(xticks)], paste("(", lefts[length(xticks)], ", ", rights[length(xticks)], "]", sep = ""), cex.axis = 0.6)
        mtext(paste("(", bmlefts[length(xticks)], ", ", bmrights[length(xticks)], "]", sep = ""), side = 1, line = 1.7, at = xticks[length(xticks)], cex = 0.6)
    }

    if (actual.scale == TRUE) {
        xticks = old.xticks
    }

    # add events and patients
    num.pats.rights <- num.pats[match(as.numeric(xticks), as.numeric(center.pt))]
    if (outcome.class == "survival") {
        num.events.lefts <- num.events[match(as.numeric(xticks), as.numeric(center.pt))]
        if (actual.scale == TRUE) {
            xticks = xticks2
        }
        mtext("Events/N:            ", side = 1, line = 2.35, at = xticks[1], adj = 1, cex = 0.6)
        mtext(paste(num.events.lefts, "/", num.pats.rights, sep = ""), side = 1, line = 2.35, at = xticks, cex = 0.6)
    } else {
        if (actual.scale == TRUE) {
            xticks = xticks2
        }
        mtext(paste("N: ", num.pats.rights, sep = ""), side = 1, line = 2.35, at = xticks, cex = 0.6)
    }

    if (actual.scale == TRUE) {
        xticks = old.xticks
    }

    # Wrap around with a box
    box()

    # Define center points as median values
    if (actual.scale == TRUE) {
        center.pt <- sdata[, "BMV.Center"]
    }

    # Actual STEPP lines are drawn at last on top of all the reference lines
    lines(center.pt, effect.size, lwd = 2, col = estimate.color)
    lines(center.pt, lower, lty = 2,  col = ci.color)
    lines(center.pt, upper, lty = 2,  col = ci.color)

    # Add shading
    if (ci.shade) {
        x <- c(center.pt[1], center.pt)
        xp <- c(x, rev(x))
        yp <- c(upper[1], upper, rev(lower), lower[length(lower)])
        polygon(xp, yp, density=NULL, col= rgb(0.6,0.6,0.6,0.1), lwd=0.5, border=NA)
    }

    # Add legend
    if (show.legend) {
        # Display All Comers reference line in legend
        if (show.refline.ac) {
            legend.text <- c(legend.text, "Biomarker Population")
            legend.col <- c(legend.col, refline.color.ac)
            legend.lty <- c(legend.lty, 1)
            legend.lwd <- c(legend.lwd, 1)
        }
        legend(legend.loc, legend = legend.text, col = legend.col, lty = legend.lty, lwd = legend.lwd, bty = legend.bty)
    }

    PlotParam()

    return(sdata2)

}
