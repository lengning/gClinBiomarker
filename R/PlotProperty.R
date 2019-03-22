#' Plot Distribution
#'
#' Plot biomarker or clinical variables properties.
#'
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}, and previous team members (see DESCRIPTION)
#'
#' @param data input data frame. Rows are patients and columns are variables (e.g. demographics variables, time to event variables,
#' biomarker variables, treatment indicator, etc.). One patient per row.
#' @param biomarker.var name of the biomarker variable. Should be in colnames of \code{data}.
#' @param biomarker.class can be either "numeric" or "categorical". If NULL (default), the class will be defined automatically.
#' @param var name of a clinical variable. It can be a vector of variables. Should be in colnames of \code{data}. Default is NULL.
#' @param var.class can be either "numeric" or "categorical". It can be a vector of classes. If NULL (default) and var is not NULL, then \code{var.class} will be defined automatically.
#' @param log2 if TRUE, computes binary (i.e. base 2) logarithm. It can be a vector if there are several numeric variables. The \code{log2} transofrmation can be applied to numeric variables only. Default is FALSE.
#' @param col the color of the line segments or dots. Default is "blue" with 30 percent transparency, i.e. \code{rgb(0, 0, 1, alpha=0.3)}.
#' @param add.num the constant to add to all values. Helps to avoid applying log transformation on 0 or negative values. Will be ignored if covariate is categorical. Default is 0.
#' @param text.font legend text font size. Default is 3.
#' @param main the main title. Default is \code{"Distribution of"}.
#' @param xlab x axis label. Default is "".
#' @param add.lab an additional text to y axis label. Default is "".
#' @param border a vector of colors for the outlines of the boxplots. Default is NULL.
#' @param add.cor add the correlation coefficient to the boxplot. Default is FALSE.
#' @param cor.method which correlation coefficient to compute. One of "pearson", "kendall", or "spearman" (default) can be abbreviated.
#' @param lowess.line performs the computations for the \code{LOWESS} smoother which uses locally-weighted polynomial regression. See \code{\link{lowess}}.
#' @param lowess.line.col the smoother color. Default is "deepskyblue".
#' @param f the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness. Default is 0.3.
#' @param show.biomarker.uni,show.clinical.uni,show.association indicate whether to show biomarker uni-variate plot, clinical variable uni-variate plot, biomarker-clinical variable association plot, respectively. Default is TRUE for all but show.clinical.uni.
#' @param las create a barplot with labels parallel (horizontal) to bars if \code{las=2}. Default is 1.
#' @param pdf.name name of output pdf file. If it's NULL (default), the plots will be displayed but not saved as pdf.
#' @param pdf.param a list of parameters that define pdf graphics device. See \code{\link{pdf}}. Default is \code{list(width=6, height=4.5)}.
#' @param par.param a list of parameters that define graphcial parameters. See \code{\link{par}}. Default is \code{list(mar=c(4,4,3,2))}.
#' @param ... other arguments passed on to the individual functions, like hist(), boxplot(), etc.
#'
#' @return If only a biomarker variable is given, it will crete a density plot for a numeric variable or bar plot for a categorical variable.
#' If only a vactor of clinical variables is provided, it will create a density plot for each numeric variable and a bar plot for each categorical variable.
#' If both a biomarker variable and a vector of clinical variables are given, it will create: a scatter plot for a numeric pair of variables;
#' a bar plot for each categorical variable; a bar plot for a pair of categorical variables; a density plot for a numeric clinical variable;
#' a boxplot for a pair of numeric and categorical variables.
#'
#' @examples
#' data(input)
#' PlotProperty(data=input, biomarker.var="KRAS.exprs", biomarker.class="numeric", log2=TRUE)
#' PlotProperty(data=input, biomarker.var="KRAS.exprs", biomarker.class="numeric", log2=TRUE, breaks=5)
#' PlotProperty(data=input, biomarker.var="KRAS.exprs", biomarker.class="numeric", var="OS", var.class="numeric", log2=c(TRUE, FALSE))
#' PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical", var=c("Arm","OS"), var.class=c("categorical", "numeric"), par.param=list(mfrow=c(3, 2)))
#' PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical", var=c("Country", "Age"), var.class=c("categorical", "numeric"), col=rgb(0, 0, 1, 0.2), par.param=list(mfrow=c(3,2)))
#'
#' @export

PlotProperty <- function(data,
                         biomarker.var,
                         biomarker.class=NULL,
                         var=NULL,
                         var.class=NULL,
                         log2=FALSE,
                         col=rgb(0, 0, 1, alpha=0.3),
                         add.num=0,
                         text.font=3,
                         main="Distribution of",
                         xlab="",
                         add.lab="",
                         border=NULL,
                         add.cor=FALSE,
                         cor.method="spearman",
                         lowess.line=FALSE,
                         lowess.line.col="deepskyblue",
                         show.biomarker.uni=TRUE, show.clinical.uni=FALSE, show.association=TRUE,
                         f=0.3,
                         las=1,
                         pdf.name=NULL,
                         pdf.param=list(width=6, height=4.5),
                         par.param=list(mar=c(4,4,3,2)), ...) {

    if (all(c(show.biomarker.uni, show.clinical.uni, show.association)==FALSE)) {
        stop("show.biomarker.uni, show.clinical.uni, show.association cannot all be FALSE!")
    }

    # Check the data
    stopifnot(class(data) == "data.frame")

    if (!is.data.frame(data)) {
        stop("An input is not a data frame!")
    }

    if (is.null(biomarker.var) & is.null(var))  {
        stop("At least one biomarker variable or clinical variable should be provided!")
    }

    if (length(biomarker.var) > 1 ) {
        stop("Only one biomarker variable should be given!")
    }

    if (!is.null(biomarker.class) & length(biomarker.class) > 1 ) {
        stop("Only one class of biomarker variable should be given!")
    }

    sel <- which(!biomarker.var %in% colnames(data))
    if (length(sel) > 0) {
        stop(paste("A covariate", biomarker.var[sel], "is not in the data frame! "))
    }

    sel <- which(!var %in% colnames(data))
    if (length(sel) > 0) {
        stop(paste("A covariate", var[sel], "is not in the data frame! "))
    }

    if (length(biomarker.class) != 0) {
        if (!all(biomarker.class %in% c("numeric", "categorical"))) {
            stop("The class of the biomarker variables can be only 'numeric' or 'categorical'!")
        }
    }

    if (length(var.class) != 0) {
        if (!all(var.class %in% c("numeric", "categorical"))) {
            stop("The class of the clinical variables can be only 'numeric' or 'categorical'. Please check spelling!")
        }
    }

    possible.class <- c("categorical", "numeric")
    if (is.null(biomarker.class)) {
        if (class(data[, biomarker.var]) %in% c("numeric", "integer")) {
            biomarker.class <- "numeric"
        }
        if (class(data[, biomarker.var]) %in% c("logical")) {
            class(data[, biomarker.var]) <- "character"
        }
        if (class(data[, biomarker.var]) %in% c("character", "factor")) {
            biomarker.class <- "categorical"
        }
    }

    if (!is.null(var) & is.null(var.class)) {
        for (i in 1:length(var)) {
            if (class(data[, var[i]]) %in% c("numeric", "integer")) {
                var.class[i] <- "numeric"
            }
            if (class(data[, var[i]]) %in% c("logical")) {
                class(data[, var[i]]) <- "character"
            }
            if (class(data[, var[i]]) %in% c("character", "factor")) {
                var.class[i] <- "categorical"
            }
        }
    }
    if (length(log2) == 1) {
        log2 <- rep(log2, length(c(biomarker.var, var)))
    }

    if (length(log2) < length(c(biomarker.var, var))) {
        stop(paste('the length of parameter log2 should be either 1 or length of c(biomarker.var, var)'))
    }

    PlotParam(pdf.name, pdf.param, par.param)

    # Case 1: Biomarker variable property or Clinical variables property
    if ((!is.null(biomarker.var) & is.null(var)) | (is.null(biomarker.var) & !is.null(var))) {
        if (!is.null(biomarker.var) & is.null(var)) {
            vars <- biomarker.var
            classes <- biomarker.class
        } else if (is.null(biomarker.var) & !is.null(var)) {
            vars <- var
            classes <- var.class
        }

        # Create a plot for each variable
        for (i in 1:length(vars)) {
            # if continues, then density plot
            if (classes[i] == "numeric") {
                V <- data[, vars[i]] + add.num

                if (log2[i] == TRUE) {
                    if (any(data[, vars[i]] <= 0, na.rm=T)) {
                        stop(paste(vars[i], " contains values less than or equal 0!
                       No log transformation possible! You may set add.num to add a constant to all values."))
                    }
                    V <- log2(V)
                }

                sd.v <- round(sd(V, na.rm=TRUE), 2)
                s.v <- summary(V)
                leg.text <- paste(c("Mean", "SD", "Median", "Range"),
                                  c(round(s.v, 2)["Mean"],
                                    sd.v,
                                    round(s.v, 2)["Median"],
                                    paste("(", round(s.v["Min."], 2), " - ", round(s.v["Max."], 2), ")", sep="")),
                                  sep=": ")

                y2 <- max(density(V, na.rm=T)$y, na.rm=TRUE)

                old.xlab <- xlab

                if (xlab == "") {
                    xlab <- paste(vars[i], ifelse(log2[i] == TRUE, "(log2 scale)", ""))
                }

                hist(V, prob=T, main=paste(main, vars[i], sep=" "), xlab=xlab, col="grey", ylim=c(0, y2), ...)

                xlab <- old.xlab

                lines(density(V, na.rm=T), lwd=2, col=col)
                legend("topright", leg.text, bty="n", text.font=text.font)
                mtext(side=3, line=0, paste("N=", sum(!is.na(V))))
                box()
                # if categorical, then barplot
            } else {
                tab <- table(data[, vars[i]])
                if (length(names(tab)) > 0) {
                    for (p in 1:length(names(tab))) {
                        if (nchar(names(tab)[p]) >= 8) {
                            names(tab)[p] <- paste0(trimws(substr(names(tab)[p], 1, 8)), "..")
                        }
                    }
                }
                freqs <- paste("(", round(100*tab/sum(tab), 2), "%)", sep="")
                barplot(tab, names.arg=paste(names(tab), freqs), main=paste(main, vars[i], sep=" "), las=las, ...)
            }
        }


        # Case 2: Biomarker variable vs Clinical variables
    } else {
        # if biomarker variable is numeric
        if (biomarker.class == "numeric") {

            BV <- data[, biomarker.var] + add.num
            j <- 1

            if (log2[j] == TRUE) {
                if (any(data[, biomarker.var] <= 0, na.rm=T)) {
                    stop(paste(biomarker.var, " contains values less than or equal 0!
                       No log transformation possible! You may set add.num to add a constant to all values."))
                }
                BV <- log2(BV)
            }
            if (show.biomarker.uni == TRUE) {
                V <- BV
                sd.v <- round(sd(V, na.rm=TRUE), 2)
                s.v <- summary(V)
                leg.text <- paste(c("Mean", "SD", "Median", "Range"),
                                  c(round(s.v, 2)["Mean"],
                                    sd.v,
                                    round(s.v, 2)["Median"],
                                    paste("(", round(s.v["Min."], 2), " - ", round(s.v["Max."], 2), ")", sep="")),
                                  sep=": ")

                y2 <- max(density(V, na.rm=T)$y, na.rm=TRUE)

                old.xlab <- xlab

                if (xlab == "") {
                    xlab <- paste(biomarker.var, ifelse(log2[j] == TRUE, "(log2 scale)", ""))
                }

                hist(V, prob=T, main=paste(main, biomarker.var, sep=" "), xlab=xlab, col="grey", ylim=c(0, y2), ...)

                xlab <- old.xlab

                lines(density(V, na.rm=T), lwd=2, col=col)
                legend("topright", leg.text, bty="n", text.font=text.font)
                mtext(side=3, line=0, paste("N=", sum(!is.na(V))))
                box()
            }

            for (i in 1:length(var)) {
                # if clinical variable is numeric, then scatterplot
                if (var.class[i] == "numeric") {
                    j <- j+1
                    V <- data[, var[i]] + add.num

                    if (log2[j] == TRUE) {
                        if (any(data[, var[i]] <= 0, na.rm=T)) {
                            stop(paste(var[i], " contains values less than or equal 0!
                                       No log transformation possible! You may set add.num to add a constant to all values."))
                        }
                        V <- log2(V)
                    }
                    if(show.clinical.uni == TRUE){
                        sd.v <- round(sd(V, na.rm=TRUE), 2)
                        s.v <- summary(V)
                        leg.text <- paste(c("Mean", "SD", "Median", "Range"),
                                          c(round(s.v, 2)["Mean"],
                                            sd.v,
                                            round(s.v, 2)["Median"],
                                            paste("(", round(s.v["Min."], 2), " - ", round(s.v["Max."], 2), ")", sep="")),
                                          sep=": ")

                        y2 <- max(density(V, na.rm=T)$y, na.rm=TRUE)

                        old.xlab <- xlab

                        if (xlab == "") {
                            xlab <- paste(var[i], ifelse(log2[j] == TRUE, "(log2 scale)", ""))
                        }

                        hist(V, prob=T, main=paste(main, var[i], sep=" "), xlab=xlab, col="grey", ylim=c(0, y2), ...)

                        xlab <- old.xlab

                        lines(density(V, na.rm=T), lwd=2, col=col)
                        legend("topright", leg.text, bty="n", text.font=text.font)
                        mtext(side=3, line=0, paste("N=", sum(!is.na(V))))
                        box()
                    }

                    x <- V
                    y <- BV
                    if(show.association==TRUE){
                        plot(x, y, ylab=paste(biomarker.var, add.lab, ifelse(log2[1] == TRUE, "(log2 scale)", "")),
                             xlab=paste(var[i], ifelse(log2[j] == TRUE, "(log2 scale)", "")),
                             main=paste(biomarker.var, "by", var[i]), col=col, ...)
                        grid(nx=NULL, ny=NULL)

                        if(lowess.line) {
                            lines(lowess(x, y, f=f), lwd=2, col=lowess.line.col)
                        }
                    }

                    # if clinical variable is categorical, then boxplot
                } else if (var.class[i] == "categorical") {
                    x <- factor(data[, var[i]])
                    xx <- jitter(as.numeric(x))
                    yy <- BV
                    nlev <- nlevels(x)
                    ylim <- range(yy, na.rm=T)

                    if(show.clinical.uni==TRUE) {
                        tab <- table(data[, var[i]])
                        if (length(names(tab)) > 0) {
                            for (p in 1:length(names(tab))) {
                                if (nchar(names(tab)[p]) >= 8) {
                                    names(tab)[p] <- paste0(trimws(substr(names(tab)[p], 1, 8)), "..")
                                }
                            }
                        }
                        freqs <- paste("(", round(100*tab/sum(tab), 2), "%)", sep="")
                        barplot(tab, names.arg=paste(names(tab), freqs), main=paste(main, var[i], sep=" "), las=las, ...)
                    }

                    if(show.association==TRUE){
                        if (is.null(col)){
                            col <- colorRampPalette(c("deepskyblue", "tomato"))(nlev)
                            col <- col[as.numeric(x)]
                        }

                        if (na.exclude(col)[1] == FALSE) {
                            col <- NULL
                        }

                        bx <- boxplot(as.formula(paste("yy ~ factor(", var[i], ")")), data=data,
                                      main=paste(biomarker.var, "by", var[i]),
                                      border=border, ylim=ylim, outline=F, axes=F,
                                      ylab=paste(biomarker.var, add.lab, ifelse(log2[1] == TRUE, "(log2 scale)", "")), ...)
                        points(xx, yy, col=col)

                        for (v in 1:length(bx$names)) {
                            if (nchar(bx$names[v]) >= 8) {
                                bx$names[v] <- paste0(trimws(substr(bx$names[v], 1, 8)), "..")
                            }
                        }

                        if (add.cor) {
                            mycor <- cor(xx, yy, method=cor.method, use="pairwise.complete")
                            legend("bottomright", paste("spear cor =", round(mycor, 2), sep=""), text.font=3)
                        }

                        axis(2)
                        box()
                    }
                    if (par("srt") != 0) {
                        sp <- ylim[2]-ylim[1]
                        axis(1, labels=F, at=1:length(bx$n), lwd=0)
                        text(1:length(bx$n), par("usr")[3] - 0.03*sp, labels=bx$names)
                    }

                    if (par("srt") == 0) {
                        axis(1, labels=bx$names, at=1:length(bx$names), lwd=0, font=2)
                    }

                    axis(3, line=-1, lwd=0, at=1:length(bx$n), paste("N=", bx$n, sep=""), cex=0.8)
                    grid(nx=NULL, ny=NULL)
                }
            }

            # if biomarker variable is categorical
        } else {
            if (show.biomarker.uni == TRUE) {
                tab <- table(data[, biomarker.var])
                if (length(names(tab)) > 0) {
                    for (p in 1:length(names(tab))) {
                        if (nchar(names(tab)[p]) >= 8) {
                            names(tab)[p] <- paste0(trimws(substr(names(tab)[p], 1, 8)), "..")
                        }
                    }
                }
                freqs <- paste("(", round(100*tab/sum(tab), 2), "%)", sep="")
                barplot(tab, names.arg=paste(names(tab), freqs), main=paste(main, biomarker.var, sep=" "), las=las, ...)
            }

            for (i in 1:length(var)) {
                j <- 1
                # if clinical variable is numeric, then hist(clin:num) and boxplot (biomar:cat + clin:num)
                if (var.class[i] == "numeric") {

                    V <- data[, var[i]] + add.num

                    if (log2[j] == TRUE) {
                        if (any(data[, var[i]] <= 0, na.rm=T)) {
                            stop(paste(var[i], " contains values less than or equal 0!
                       No log transformation possible! You may set add.num to add a constant to all values."))
                        }
                        V <- log2(V)
                    }

                    # hist(clin:num)
                    if (show.clinical.uni == TRUE) {
                        sd.v <- round(sd(V, na.rm=TRUE), 2)
                        s.v <- summary(V)
                        leg.text <- paste(c("Mean", "SD", "Median", "Range"),
                                          c(round(s.v, 2)["Mean"],
                                            sd.v,
                                            round(s.v, 2)["Median"],
                                            paste("(", round(s.v["Min."], 2), " - ", round(s.v["Max."], 2), ")", sep="")),
                                          sep=": ")

                        y2 <- max(density(V, na.rm=T)$y, na.rm=TRUE)

                        old.xlab <- xlab

                        if (xlab == "") {
                            xlab <- paste(var[i], ifelse(log2[j] == TRUE, "(log2 scale)", ""))
                        }

                        hist(V, prob=T, main=paste(main, var[i], sep=" "), xlab=xlab, col="grey", ylim=c(0, y2), ...)

                        xlab <- old.xlab

                        lines(density(V, na.rm=T), lwd=2, col=col)
                        legend("topright", leg.text, bty="n", text.font=text.font)
                        mtext(side=3, line=0, paste("N=", sum(!is.na(V))))
                        box()
                    }
                    # boxplot (biomar:cat + clin:num)
                    if (show.association == TRUE) {
                        x <- factor(data[, biomarker.var])
                        xx <- jitter(as.numeric(x))
                        yy <- V
                        nlev <- nlevels(x)
                        ylim <- range(yy, na.rm=T)

                        if (is.null(col)){
                            col <- colorRampPalette(c("deepskyblue", "tomato"))(nlev)
                            col <- col[as.numeric(x)]
                        }

                        if (na.exclude(col)[1] == FALSE) {
                            col <- NULL
                        }

                        bx <- boxplot(as.formula(paste("yy ~ factor(", biomarker.var, ")")), data=data,
                                      main=paste(var[i], "by", biomarker.var),
                                      border=border, ylim=ylim, outline=F, axes=F,
                                      ylab=paste(var[i], add.lab, ifelse(log2[j] == TRUE, "(log2 scale)", "")), ...)
                        points(xx, yy, col=col)

                        for (v in 1:length(bx$names)) {
                            if (nchar(bx$names[v]) >= 8) {
                                bx$names[v] <- paste0(trimws(substr(bx$names[v], 1, 8)), "..")
                            }
                        }

                        if (add.cor) {
                            mycor <- cor(xx, yy, method=cor.method, use="pairwise.complete")
                            legend("bottomright", paste("spear cor =", round(mycor, 2), sep=""), text.font=3)
                        }

                        axis(2)
                        box()

                        if (par("srt") != 0) {
                            sp <- ylim[2]-ylim[1]
                            axis(1, labels=F, at=1:length(bx$n), lwd=0)
                            text(1:length(bx$n), par("usr")[3] - 0.03*sp, labels=bx$names)
                        }

                        if (par("srt") == 0) {
                            axis(1, labels=bx$names, at=1:length(bx$names), lwd=0, font=2)
                        }

                        axis(3, line=-1, lwd=0, at=1:length(bx$n), paste("N=", bx$n, sep=""), cex=0.8)
                        grid(nx=NULL, ny=NULL)

                        j <- j+1
                    }
                    # if clinical variable(s) is (are) categorical, then barplot for each + barplot of interaction
                } else if (var.class[i] == "categorical") {
                    # if (show.biomarker.uni == TRUE) {
                    #     tab <- table(data[, biomarker.var])
                    #     freqs <- paste("(", round(100*tab/sum(tab), 2), "%)", sep="")
                    #     barplot(tab, names.arg=paste(names(tab), freqs), main=paste(main, biomarker.var, sep=" "), las=las)
                    # }

                    if (show.clinical.uni == TRUE) {
                        tab <- table(data[, var[i]])
                        if (length(names(tab)) > 0) {
                            for (p in 1:length(names(tab))) {
                                if (nchar(names(tab)[p]) >= 8) {
                                    names(tab)[p] <- paste0(trimws(substr(names(tab)[p], 1, 8)), "..")
                                }
                            }
                        }
                        freqs <- paste("(", round(100*tab/sum(tab), 2), "%)", sep="")
                        barplot(tab, names.arg=paste(names(tab), freqs), main=paste(main, var[i], sep=" "), las=las, ...)
                    }

                    if (show.association == TRUE) {
                        tab <- table(data[, biomarker.var], data[, var[i]])
                        if (length(names(tab[1, ])) > 0) {
                            for (p in 1:length(names(tab[1, ]))) {
                                if (nchar(names(tab[1, ])[p]) >= 8) {
                                    colnames(tab)[p] <- paste0(trimws(substr(names(tab[1, ])[p], 1, 8)), "..")
                                }
                            }
                        }
                        freqs <- paste("(", round(100*tab[1, ]/sum(tab[1, ]), 2), "%)", sep="")
                        barplot(tab, names.arg=paste(names(tab[1, ]), freqs),
                                main=paste(main, biomarker.var, "by", var[i], sep=" "),
                                beside=TRUE, las=las, ...)
                        legend("topleft", legend = names(tab[, 1]), fill=c("black", "grey"), cex=0.8)
                    }
                }
            }
        }
    }

    PlotParam()
}
