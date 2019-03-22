#' Box/Box-Percentile Plot with additional vertical stripchart, mean(s), and trend-lines.
#'
#' This function produces boxplots as function \code{boxplot} and box-percentile plots as function bpplot (package Hmisc),
#' offering a formula interface, data specification via a vector of variable names, as numeric matrix where all columns will be interpreted as
#' variables or by specifying one or multiple numeric vectors (see examples below). It additionally/optionally
#' adds vertical jitterplots to each boxplot, adds the number of observations, adds mean-values, and allows to add a trend-line connecting either
#' mean-values (trend="mean") or median-values (trend="median"). Furthermore, this function filters for groups with less than 'threshold'
#' observations. For these (sub-)groups no boxplots will be drawn but a jitterplot will be produced. Labels for each group (X-axis labels) can be customized via
#' 'xlab.srt', 'xlab.cex', 'xlab.col', 'xlab.font'. Note: When using the formula-interface, one cannot specify formulas as in function 'boxplot'
#' i.e. using something like 'a + b' to get factor crossing, one has to use the correct formula instead as used in e.g. 'lm', i.e. 'a:b'.
#'
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}, Vinzent Rolny \email{vinzent.rolny@@roche.com}, christina Rabe \email{rabe.christina@@gene.com}
#'
#' @param ...           a numeric matrix, multiple numeric vectors or a (possibly named) list with multiple numeric vectors to be plotted and additional
#'                      graphical parameters (identified by names) to be passed on.
#' @param obj           (data.frame, matrix) corresponding to the dataset to be used.
#' @param form          (formula) such as 'y~grp', where 'y' is a numeric vector of data values to be split into
#'                      groups according to the grouping variable 'grp'.
#' @param var           (character) vector specifying the columns in 'obj' to be used for plotting. The order of elements is retained in the boxplot.
#' @param box.type      (character) "b" = for regular boxes, "bp" = for box-percentile boxes
#' 						(see ? bpplot of package \code{Hmisc}).
#' @param horizontal	(logical) TRUE = the boxes are drawn horizontally, FALSE = vertical boxplot. Note, this has no effect when box.type="bp", since
#' 						the underlying function bpplot of the \code{Hmisc} package cannot do it. Also note, that \code{horizontal=TRUE} automatically
#'                      changes the meaning of X- and Y-axis. All arguments referring to Y-axis will then implicitely mean the X-axis and vice versa.
#' @param col           (character) string(s) specifying colors to be used to color the bodies of the box plots. By default they are in the background color.
#' @param Xaxis         (list) passed to function 'axis' which allows to fully specify the X-axis labelling (see ?axis). Default is to determine group labels
#'                      automatically and plotting group labels below the plot. Set to NULL for omitting group labels.
#'                      For custom group-labels use the list-element "labels". The order of these labels must correspond to the order of group-levels in case the
#'                      formula interface was used, i.e. check the order of sort(unique(obj$grp)). To obtain your desired ordering, specify classification variables as
#'                      factor-objects using 'factor(dat$trt, levels=c("Level_1", "Level_2", ..., "Level_n"))' where "trt" is the classification variable and Level_1 to Level_n
#'                      represents your desired order of factor levels (which is used in the plot).
#' @param Xaxis2        (list) passed to function 'axis' which allows to fully specify the X-axis labelling above the plot (see ?axis). This is used to specify the number
#'                      of observations per sub-group according to the group-labels. Default is to plot these numbers above the plot. To omit set to NULL.
#'  					Note, if 'horizontal=TRUE' this axis-labelling will appear in the right margin and labels will be oriented perpendicular to the axis.
#'                      Set 'Xaxis2=list(las=0)' to overwrite the default for usual parallel orientation.
#' @param XaxisTab		(list) or NULL, if a (possibly empty) list AND 'form' is specified, 'Xaxis' will be set to NULL and a table will be used a X-axis label
#' 						representing the combination of factor-levels for all factors appearing in 'form' (see example document and see \code{\link{addXlabTable}}).
#'                      If 'horizontal=TRUE' the table will appear in the left margin as any other setting originally intended to apply to an X-axis element.
#' 						Note, there has to be enough space in the bottom margin for adding the table (use "mai" or "mar" in \code{\link{par}}).
#' 						Your may specify 'XaxisTab' as list with two sub-lists "Label" and "Text", which will then be evaluated for rownames of the table ("Label")
#'                      and the text in the cells of the table ("Text") separately (see examples of function \code{\link{addXlabTable}} for details).
#' @param Yaxis         (list) passed to function 'axis' which allows to fully specify the the appearence of the Y-axis (see ?axis).
#' @param Ylabel        (list) passed to function 'mtext' which can be used to fully specify the Y-axis label.
#' @param Xlabel        (list) passed to function 'mtext' which can be used to fully specify the X-axis label or even group labels (see example document).
#' @param Title         (list) passed to function 'title'.
#' @param Grid          (list) passed to function 'grid'. Set to NULL to omit (default). For adding a grid simply set 'Grid=TRUE' which applies default settings
#'                      for the added grid or fully specify the grid to be added by specifying each argument (see ?addGrid for details).
#' @param transf        (function) name of a function to be used to tranform data before plotting, e.g. log, log10 or user-defined functions,
#'                      i.e. func=function(x){x[x==0]<-min(x, na.rm=TRUE)/2} and setting 'transf=func'.
#' @param sc.col        (character) string or vector of strings specifying color(s) of plotting symbols in the stripchart.
#' 						By default "black" with 80\% transparency (alpha=.2). In case of a vector of strings, it is best practise to provide as
#' 						many elements as there are elements in 'obj', which can be used to highlight an additional grouping factor within a single box.
#' 						Note, it is the responsibility of the user to provide 'sc.col', 'sc.pch' or 'sc.cex' in such a way, that the graphical output
#'                      is meaningful. In case of providing less than 'nrow(obj)' elements, 'sc.col' will be replicated to the required length.
#' @param sc.pch        (integer) value or vector of integers specifying plotting smybols for the stripchart. The same rules as for 'sc.col' and 'sc.cex' apply here.
#' @param sc.cex        (numeric) value or vector of numeric values specifying the magnification of plotting symbols in the stripchart. The same rules as for 'sc.col' and 'sc.cex' apply here.
#' @param sc.jitter     (numeric) specifying the amount of jittering in the stripchart.
#' @param trend         (character) "mean" = mean values are connected by a line emphasizing the dynamics (especially for time course data)
#'                                  "median" = median values are used, set to NULL to omit.
#' @param trend.lty     (integer) line type of the trend line.
#' @param trend.lwd     (numeric) line width of the trend line.
#' @param trend.col     (character) color of the trend line.
#' @param threshold     (integer) minimum number of points required for plotting a boxplot, otherwise only the stripchart will be plotted.
#' @param border        (character) string specifying the border color(s) of boxes. This can be a vector with different colors for multiple boxes.
#' @param mean.pch      (integer) plotting symbol for mean-values, which are added to the plot. Use '-1' to prevent plotting of mean-values.
#' @param mean.cex      (numeric) specifying the magnification of mean-value plotting symbols.
#' @param mean.col      (character) specifying the color of mean-value plotting symbols.
#' @param mean.lwd      (integer) specifying the line width of mean-value plotting symbols (cross, plus, asterisk etc.).
#' @param vline         (numeric) value(s) specifying vertical lines added to the plot.
#' @param vl.lwd        (integer) line width of vertical lines.
#' @param vl.lty        (integer) line type of vertical lines.
#' @param vl.col        (character) color of vertical lines.
#' @param Box           (logical) TRUE = a box is plotted surrounding the plot, FALSE = no box.
#'
#' @examples
#' data(input)
#'
#' ## specify variables to plot
#' BoxPlot(input, KRAS.exprs~Arm, transf = log)
#'
#' ## same plot, now horizontally plotted
#' BoxPlot(input, KRAS.exprs~Arm, transf = log, horizontal=TRUE,
#'        Xaxis=list(las=2, hadj=2), Xaxis2=list(las=2, hadj=-.25))
#'

#' @export

BoxPlot <- function(..., obj, form=NULL, var=NULL, box.type="b",
		horizontal=FALSE, col="white",
		sc.col="#00000040", sc.pch=15L, sc.cex=1, sc.jitter=.1,
		Xaxis=list(side=1, mgp=c(3, .5, 0), font=1, tick=TRUE),
		Xaxis2=list(side=3, font=2, mgp=c(3, .5, 0), tick=TRUE),
		XaxisTab=NULL,
		Yaxis=list(side=2, mgp=c(3,1,0)),
		Ylabel=list(text="", side=2, line=2.5, font=1, cex=1),
		Xlabel=list(text="", side=1, line=2.5, font=1, cex=1),
		Title=list(line=2.5), Grid=NULL, transf=NULL,
		trend=NULL, trend.lty=1L, trend.lwd=1, trend.col="blue", threshold=5L,
		border="black", mean.pch=3L, mean.cex=1.5, mean.col="salmon", mean.lwd=2L,
		vline=NULL, vl.lwd=1L, vl.lty=1L, vl.col="black", Box=TRUE )
{
	box.type <- match.arg(tolower(box.type), c("b", "bp"))

	stopifnot(is.logical(horizontal))
	if(box.type == "bp" && horizontal)
		stop("Box-percentile plots cannot be plotted horizontally!")

	args <- list(...)

	if(length(args) != 0)                                                       # additional arguments specified
	{
		if(is.null(names(args)) || any(names(args) == ""))                      # split args into named and unnamed arguments
		{
			if(is.null(names(args)))                                            # all arguments unnamed
			{
				unargs <- args
				args <- NULL
			}
			else
			{
				unargs <- args[which(names(args) == "")]                        # note which args are unnamed
				args <- args[which(names(args) != "")]
				if("xlab" %in% names(args))
					args <- args[-which(names(args) == "xlab")]
				if("ylab" %in% names(args))
					args <- args[-which(names(args) == "ylab")]
			}
		}
		else
			unargs <- NULL

		if(!"mar" %in% names(args))												# do not overwrite user-specification of mar
		{
			if(horizontal)
				Mar <- c(4.1, 4.1, 4.1, 4.1)
			else
				Mar <- par("mar")
		}
		else
		{
			Mar <- args[["mar"]]
			args <- args[-which(names(args) == "mar")]
		}

		if(horizontal)
			old.par <- par(mar=Mar, xaxs="r", yaxs="i")
		else
			old.par <- par(mar=Mar, xaxs="i", yaxs="r")

		cls <- unlist(lapply(unargs, class))

		if(any(cls == "data.frame"))                                            # data.frame specified (arg 'obj')
		{
			obj <- unargs[[which(cls == "data.frame")]]
			if(any(cls == "formula") && is.null(form))
				form <- unargs[[which(cls == "formula")]]
		}
		else if(any(cls == "matrix"))                                           # numeric/integer matrix specified
		{
			mat <- unargs[[which(cls == "matrix")]]
			if(class(mat[1,1]) %in% c("integer", "numeric"))                    # numeric matrix --> use columns as groups
			{
				if(is.null(colnames(mat)))
					colnames(mat) <- 1:ncol(mat)
				obj <- as.data.frame(mat)
				var <- colnames(obj)
			}
		}
		else if(any(cls %in% c("integer", "numeric")))                          # multiple numeric/integer vectors specified
		{
			ind <- which(cls %in% c("integer", "numeric"))
			len <- unlist(lapply(unargs[ind], length))
			tmp <- matrix(nrow=max(len), ncol=length(unargs[ind]))
			for(i in 1:ncol(tmp))
			{
				vec <- unargs[[ind[i]]]
				tmp[1:length(vec),i] <- vec
			}
			colnames(tmp) <- 1:length(ind)
			obj <- as.data.frame(tmp)
			var <- colnames(obj)
		}
		else if(any(cls == "list"))
		{
			tmpL <- unargs[[which(cls == "list")]]
			namesL <- names(tmpL)
			if( all( sapply(tmpL, class) %in% c("integer", "numeric")) )
			{
				len <- sapply(tmpL, length)
				tmp <- matrix(nrow=max(len), ncol=length(tmpL))
				for(i in 1:ncol(tmp))
				{
					vec <- tmpL[[i]]
					tmp[1:length(vec), i] <- vec
				}
				if(is.null(namesL))
					colnames(tmp) <- 1:ncol(tmp)
				else
					colnames(tmp) <- namesL
				obj <- as.data.frame(tmp)
				var <- colnames(tmp)
			}
		}

		if(any(cls == "formula"))
			form <- unargs[[which(cls == "formula")]]
	}
	else
		old.par <- par()

	if(!is.null(form) && class(form)=="formula")
	{
		if(grepl("\\+", as.character(form)[3]))
			stop("One must not use the '+' operator in formulas! Use ':' instead for crossing multiple factors!")
		tmp.var <- apply(attributes(terms(form))$factors, 1, sum)           # get var names + status
		dep.var <- names(tmp.var[which(tmp.var==0)])
		tmp.var <- names(tmp.var[-which(tmp.var==0)])                       # removes dependent variable (0)
		if(!is.null(transf))
		{
			if(class(transf) == "function")
				obj[,dep.var] <- do.call(transf, list(obj[,dep.var]))       # apply data transformation/manipulation
			else
				warning("Data transformation could no be applied! 'transf' was not correctly specified!")
		}
		obj[,tmp.var] <- lapply(obj[,tmp.var, drop=FALSE], factor)          # convert grouing variables to factors
		bp <- boxplot(form, data=obj, plot=FALSE)
		fit <- lm(formula(paste(deparse(form), "-1", sep="")), obj)
		Means <- coef(fit)
		if(length(tmp.var) == 1)                                            # just one grouping variable
		{
			names(Means) <- gsub(tmp.var, "", names(Means))                 # names now refer to factor levels
			tmp <- rep(NA, length(levels(obj[,tmp.var])))
			names(tmp) <- levels(obj[,tmp.var])
			tmp[names(Means)] <- Means
			Means <- tmp
		}

		if( !is.null(XaxisTab) && is.list(XaxisTab) )
		{
			Xaxis <- NULL
			mat.Xtab <- getMatrix(form, bp$names)
		}
	}
	else
	{
		if(is.null(var))
		{
			warning("Boxplot cannot be drawn because there is neither formula 'form' nor variable names 'var' provided!")
			return(1)
		}
		stopifnot(all(var %in% colnames(obj)))
		if(!is.null(transf))                                                # apply data-transformation
		{
			if(class(transf) == "function")
				obj[, var] <- lapply(obj[,var, drop=FALSE], transf)         # apply data transformation/manipulation
			else
				warning("Data transformation could no be applied! 'transf' was not correctly specified!")
		}
		bp <- boxplot(obj[,var], plot=FALSE)
		Means <- apply(obj[,var, drop=FALSE], 2, mean, na.rm=TRUE)
	}

	if(is.null(args) || !"ylim" %in% names(args))
	{

		args$ylim <- range( c(c(bp$stats), bp$out), na.rm=TRUE )            # NAs have to be removed for 'ylim'
	}


	Nbox <- length(bp$n)
	if(Nbox != length(border))                                              # replicate border colors
	{
		border <- rep(border, ceiling(Nbox/length(border)))[1:Nbox]
		col <- rep(col, ceiling(Nbox/length(col)))[1:Nbox]
	}
	if(any(bp$n < threshold))                                               # prevent boxes from being plotted if n < 'threshold'
	{
		border[bp$n < threshold] <- "white"
		col[bp$n < threshold] <- "white"
	}
	if(box.type == "b")
	{
		ARGS <- list(bp, show.names=FALSE, axes=FALSE, main=NA,             # combine required arguments with additional graphical args from '...'
				border=border, boxfill=col, outline=FALSE,
				horizontal=horizontal)
		ARGS <- c(ARGS, args)
		do.call(bxp, ARGS)                                                  # actually plot the boxplot
	}
	else
	{
		ARGS <- list(obj=obj, form=form, var=var, col=col, ylab=NA,
				labels=NA, add=FALSE, main=NA, line.col=border,
				add.xlab=FALSE)
		ARGS <- c(ARGS, args)
		do.call(BPplot, ARGS)
	}
	if(!is.null(Grid))
	{
		grid.default <- list(x=NULL, y=NULL, col="lightgray", lty=2L, lwd=1L)
		grid.default[names(Grid)] <- Grid
		Grid <- grid.default
		if(horizontal)
		{
			tmp <- Grid$x
			Grid$x <- Grid$y
			Grid$y <- tmp
		}
		do.call("addGrid", Grid)
	}

	if(!is.null(Xlabel))
	{
		xlab.default <- list(text="", side=1, line=2.5, font=1, cex=1)              # change specified parameters
		xlab.default[names(Xlabel)] <- Xlabel
		Xlabel <- xlab.default
		if(horizontal)
			Xlabel$side=2
		do.call(mtext, Xlabel)
	}
	if(!is.null(Ylabel))
	{
		ylab.default <- list(text="", side=2, line=2.5, font=1, cex=1)   # change specified parameters
		ylab.default[names(Ylabel)] <- Ylabel
		Ylabel <- ylab.default
		if(horizontal)
			Ylabel$side <- 1
		do.call(mtext, Ylabel)
	}

	if(!is.null(Xaxis))
	{
		xaxis.default <- list(side=1, at=1:length(bp$n), labels=bp$names, mgp=c(3,.5,0), font=1, tick=TRUE)
		xaxis.default[names(Xaxis)] <- Xaxis
		Xaxis <- xaxis.default
		if(horizontal)
			Xaxis$side=2
		do.call("axis", Xaxis)
	}

	if(!is.null(Xaxis2))
	{
		xaxis2.default <- list(side=3, at=1:length(bp$n), labels=paste("N", bp$n, sep="="),
				mgp=c(3,.5,0), font=2, tick=TRUE, las=ifelse(horizontal, 1, 0))
		xaxis2.default[names(Xaxis2)] <- Xaxis2
		Xaxis2 <- xaxis2.default
		if(horizontal)
			Xaxis2$side=4
		do.call("axis", Xaxis2)
	}

	if(!is.null(XaxisTab))
	{
		Label <- Text <- NULL

		if("Text" %in% names(XaxisTab) && is.list(XaxisTab$Text))
			Text <- XaxisTab$Text

		if("Label" %in% names(XaxisTab) && is.list(XaxisTab$Label))
			Label <- XaxisTab$Label

		addTableToMargin(mat=mat.Xtab, margin=ifelse(horizontal, "left", "bottom"),
				Label=Label, Text=Text, reorder=ifelse(horizontal, TRUE, FALSE))
	}

	if(!is.null(Yaxis))
	{
		yaxis.default <- list(side=2, mgp=c(3,1,0))
		yaxis.default[names(Yaxis)] <- Yaxis
		Yaxis <- yaxis.default
		if(horizontal)
			Yaxis$side <- 1
		do.call("axis", Yaxis)
	}

	if(!is.null(Title))
	{
		title.default <- list(main=ifelse(box.type=="b", "Box Plot", "Box-Percentile Plot"), line=2.5)
		title.default[names(Title)] <- Title
		Title <- title.default
		do.call("title", Title)
	}

	if(Box)
		box()

	if (!is.null(form) && class(form) == "formula") 				### stripchart-argument evaluations
		sc.exprs <- "stripchart(form, data=obj.temp"				# formula provided
	else
		sc.exprs <- "stripchart(obj.temp[, var]"					# non-formula

	lcol <- length(sc.col)
	lpch <- length(sc.pch)
	lcex <- length(sc.cex)

	if(	(lcol != lpch && all( c(lcol, lpch) > 1) ) || 				# any pair has more than one level but is of differnt length --> do not fit together
			(lcol != lcex && all( c(lcol, lcex) > 1) ) ||
			(lpch != lcex && all( c(lpch, lcex) > 1) ) )
	{
		warning("Parameter settings of 'sc.col', 'sc.pch' and 'sc.cex' probably result in missleading graphical output!")
	}

	sc.col <- rep(sc.col, ceiling(nrow(obj)/lcol))					# replicate arguments to sufficient length
	sc.pch <- rep(sc.pch, ceiling(nrow(obj)/lpch))
	sc.cex <- rep(sc.cex, ceiling(nrow(obj)/lcex))

	sc.combi <- data.frame( col=sc.col, pch=sc.pch, cex=sc.cex,		# determine unique combinations of all 3 sc-arguments
			stringsAsFactors=FALSE)

	sc.combi <- unique(sc.combi)

	sc.exprs <- paste(	sc.exprs,
			"col=sc.combi$col[i]",
			"pch=sc.combi$pch[i]",
			"cex=sc.combi$cex[i]", sep=", ")

	sc.exprs <- paste(	sc.exprs,
			"method = \"jitter\", vertical = !horizontal",
			"jitter = sc.jitter, add = TRUE)", sep=", ")

	for(i in 1:nrow(sc.combi))
	{
		obj.temp <- obj[sc.col == sc.combi$col[i] &					# sc-grouping according to combination of 3 sc-arguments
						sc.pch == sc.combi$pch[i] &
						sc.cex == sc.combi$cex[i], , drop=FALSE]
		eval(parse(text=sc.exprs))									# color or use symbols according to user-specification
	}

	if(!is.null(trend))
	{
		if(tolower(trend) == "mean")
		{
			if(horizontal)
				lines(na.omit(Means), which(!is.na(Means)), lty=trend.lty, col=trend.col, lwd=trend.lwd)
			else
				lines(which(!is.na(Means)), na.omit(Means), lty=trend.lty, col=trend.col, lwd=trend.lwd)
		}

		if(tolower(trend) == "median")
		{
			if(horizontal)
				lines(na.omit(bp$stats[3,]), which(!is.na(bp$stats[3,])), lty=trend.lty, col=trend.col, lwd=trend.lwd)
			else
				lines(which(!is.na(bp$stats[3,])), na.omit(bp$stats[3,]), lty=trend.lty, col=trend.col, lwd=trend.lwd)
		}
	}
	if(!is.null(vline) && class(vline) %in% c("integer", "numeric"))
	{
		if(horizontal)
			abline(h=vline, lwd=vl.lwd, lty=vl.lty, col=vl.col)
		else
			abline(v=vline, lwd=vl.lwd, lty=vl.lty, col=vl.col)
	}

	if (is.na(Means[1]) & names(Means[1]) == "") {
	    Means <- Means[2:length(Means)]
	}

	if(horizontal)
		points(Means, 1:length(bp$n), pch=mean.pch, col=mean.col, cex=mean.cex, lwd=mean.lwd)
	else
		points(1:length(bp$n), Means, pch=mean.pch, col=mean.col, cex=mean.cex, lwd=mean.lwd)
	suppressWarnings(par(old.par))
	dif <- nrow(obj) - sum(bp$n)
	if (dif > 0) {
	    message("Number of samples removed due to missing Y: ", dif)
	}
	invisible(bp)
}
