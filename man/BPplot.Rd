% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BPplot.R
\name{BPplot}
\alias{BPplot}
\title{Enhanced Box-Percentile Plots.}
\usage{
BPplot(obj, form = NULL, var = NULL, main = "Box-Percentile Plot",
  add = FALSE, col = NULL, labels = NULL, ylab = NULL,
  ylab.line = 2.5, ylab.font = 1L, swag = 0.9, line.col = "black",
  line.lwd = 1L, line.lty = 1L, add.xlab = TRUE, ...)
}
\arguments{
\item{obj}{(data.frame, matrix) corresponding to the dataset to be used}

\item{form}{(formula) such as \code{y~grp}, where \code{y} is a numeric vector of data values to be split into
groups according to the grouping variable \code{grp}.}

\item{var}{(character) vector specifying the columns in 'obj' to be used for plotting. The order of elements is retained in the boxplot.}

\item{main}{(character) string giving the main title of the plot, left out if \code{add=TRUE}}

\item{add}{(logical) TRUE=add box-percentile plot(s) to an existing plot}

\item{col}{(character) vector specifying the colors of percentile-boxes, defaults to no color}

\item{labels}{(character) string(s) for group labels, which are drawn under each boxplot. Note: The order must correspond to the order of group-levels
if the formula interface was used, i.e. check the order of sort(unique(obj$grp)).}

\item{ylab}{(character) character string specifying the optional label of the Y-axis (vertically centered)}

\item{ylab.line}{(numeric) specifying the line 'ylab' is put (useful for custom figure margins).}

\item{ylab.font}{(integer) specifying the font to be used for the Y-axis label}

\item{swag}{(numeric) values within ]0,1] specifying the swaging of percentile boxes, for better separating boxes from each other}

\item{line.col}{(character) color of boundry lines of the percentile boxes}

\item{line.lwd}{(integer) line width of the boundry lines of the percentile boxes}

\item{line.lty}{(integer) line type of the boundry lines of the percentile boxes}

\item{add.xlab}{(logical) TRUE = automically determined group-labels are plotted below each box/percentile-box}

\item{...}{additional graphical parameters passed on}
}
\description{
Please use function \code{BoxPlot} with \code{box.type="bp"} for box-percentile plots, since this function here
was designed as helper function for function \code{BoxPlot}.
It uses the R-code of function \code{bpplot} which generates box-percentile plots as implemented in 
package \code{Hmisc}. This R-code is changed and augmented to implement an formula-interface and to be able to 
add box-percentile plots to an existing plot. It is assumed that \code{form} is a simple formula object, 
e.g. \code{y~grp:time} with just simple interactions or even simpler with just a single grouping factor, e.g. \code{y~grp}.
}
\seealso{
\link{BoxPlot}, function \code{bpplot} in package \code{Hmisc}
}
\author{
Andre Schuetzenmeister (using source-code of function \code{bpplot} from package \code{Hmisc})
}
