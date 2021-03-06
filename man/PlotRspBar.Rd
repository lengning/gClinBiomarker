% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotRspBar.R
\name{PlotRspBar}
\alias{PlotRspBar}
\title{Generate bar plot to summarize response outcome by subgroups}
\usage{
PlotRspBar(data, outcome.var, binary = FALSE, rsp.response = c("CR", "PR"),
  rsp.nonresponse = c("SD", "PD", "NON CR/PD", "NE"), rsp.levels = c("CR",
  "PR", "SD", "NON CR/PD", "PD", "NE"), rsp.name = "Rsp", col = NULL,
  plot.count = FALSE, digits = 1, trt = NULL, trt.name = NULL,
  show.combine.trt = TRUE, compare.bep.itt = FALSE, bep = NULL,
  bep.name = NULL, itt.name = "Full", bep.indicator = 1,
  compare.var = FALSE, var = NULL, var.name = NULL,
  show.combine.var = TRUE, horiz = FALSE,
  main = "Association of response rate", sub = NULL, cex = 1,
  pdf.name = NULL, pdf.param = NULL, par.param = NULL)
}
\arguments{
\item{data}{Input data frame. Rows are patients and columns are variables (e.g. demographics variables, time to event variables, 
biomarker variables, treatment indicator, etc.). One patient per row.}

\item{outcome.var}{column name of the response variable (e.g. best overall response).
entries with missing rsp value will be ignored in analysis.}

\item{binary}{whether summarize the response rate by responder/nonresponder (if binary=TRUE),
or by its original category (e.g. PD/SD/PR/CR). If binary is TRUE, responder categories
and nonresponder categories should be specified in rsp.response and rsp.nonresponse (all values in the outcome column
should be included in rsp.response and rsp.nonresponse),
at the same time rsp.levels will be ignored. If binary is FALSE, rsp.levels should be
specified to provide order of the categories. At the same time rsp.response and
rsp.nonresponse will be ignored.}

\item{rsp.response}{categories that should be considered as responder.}

\item{rsp.nonresponse}{categories that should be considered as non responder.}

\item{rsp.levels}{vector that indicates how to sort the response categories. This parameter will be ignored if binary is TRUE.}

\item{rsp.name}{Display name for the responders. Default is "rsp".}

\item{col}{color for different categories}

\item{plot.count}{default is FALSE. By default percentages will be shown.
If it is TRUE, will show counts instead}

\item{digits}{see \code{\link{pdf}}. Number of digits to display for the response rate.}

\item{trt}{column name that indicates treatment variable. If this is not NULL,
the program will summarize the response outcome within each treatment arm.}

\item{trt.name}{treatment name to display}

\item{show.combine.trt}{if this is TRUE, will show summarization of combined data across treatment arms as well}

\item{compare.bep.itt}{whether want to generate two groups of bars to compare the summary statistics
in full population (e.g. ITT) vs in BEP. If this is TRUE, paramemeter bep should be specified. If this is FALSE,
parameters bep, bep.name, itt.name, bep.indicator will be ignored}

\item{bep}{name of the column which indicates biomarker evaluable population. If it is null, patients who have non NA records
in biomarker variable will be used as BEP.}

\item{bep.name}{preferred display name of the biomarker evaluable population.
If it is NULL, bep will be used.}

\item{itt.name}{preferred display name of the full population (e.g. ITT).
If it is NULL, "All" will be used.}

\item{bep.indicator}{In the subpopulation column, which value is used
to define the biomarker evaluable population.}

\item{compare.var}{whether want to generate multiple groups of bars to compare the summary statistics
in subgroups defined by var (e.g. categorical demographics or biomarker variable)
. If this is TRUE, paramemeter var should be specified. If this is FALSE,
parameters var, var.name, show.combine.var will be ignored}

\item{var}{a vector of covariate names - the clinical covariate to test}

\item{var.name}{preferred display names of the clinical covariates 
If it is NULL, var will be used.}

\item{show.combine.var}{if this is TRUE, will show summarization of combined data across var levels as well}

\item{horiz}{whether show barplot horizontally}

\item{main}{main title of the barplot. Default is "Association of biomarker expression with response rate".}

\item{sub}{footnote under the bar plot. Default is NULL.}

\item{cex}{cex}

\item{pdf.name}{name of output pdf file. If it's NULL, the plots will be displayed but not saved as pdf. Default is "rsp.barplot.pdf".}

\item{pdf.param}{list of parameters that define pdf graphics device. See \code{\link{pdf}}. Default is \code{list(width=6, height=4.5)}.}

\item{par.param}{list of parameters that define graphcial parameters. See \code{\link{par}}. Default is \code{list(mar=c(4,4,3,2))}.}
}
\value{
a bar plot that summarizes response outcome
}
\description{
This function generates bar plots to compare response outcome summarization (e.g. response rate) across different groups
(e.g. treatment arm, biomarker variable, demographics variable, full population vs bep, etc.)
}
\examples{
data(input)
PlotRspBar (input, outcome.var="Response",
           binary=TRUE,
           rsp.response=c("CR","PR"),
           rsp.nonresponse=c("SD", "PD","NON CR/PD","NE",NA),
           trt="Arm",
           compare.var=TRUE,
           var="KRAS.mutant")

}
\author{
Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}, and previous team members (see DESCRIPTION)
}
