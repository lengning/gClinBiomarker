% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotTabForestBiomarker.R
\name{PlotTabForestBiomarker}
\alias{PlotTabForestBiomarker}
\title{Generate forest plot and summarization table for 2-arm or within-arm comparison}
\usage{
PlotTabForestBiomarker(data, outcome.class = c("survival", "binary"),
  outcome.var, trt = NULL, var = NULL, var.class = NULL,
  var.name = NULL, percentile.cutoff = 0.5, numerical.cutoff = NULL,
  greater = TRUE, less = FALSE, greater.by.less = FALSE,
  across.and.within = FALSE, equal.in.high = TRUE, within.bin = FALSE,
  show.itt = TRUE, show.bep = TRUE, bep = NULL, bep.name = "BEP",
  itt.name = "All", bep.indicator = 1, covariate = NULL, strata = NULL,
  tabforest = FALSE, quantile.type = 2, digits = 2, cutoff.digits = 2,
  placebo.code = NULL, active.code = NULL, rsp.cat = TRUE,
  rsp.response = c("CR", "PR"), rsp.nonresponse = c("SD", "PD", "NON CR/PD",
  "NE", NA), var.code = NULL, surv.conf.type = "plain", ties = "efron",
  alpha = 0.05, main = NULL, main.prefix = NULL, sub = NULL,
  clip = NULL, xticks = NULL, xticks.digits = 1, xlab = NULL,
  cex.headings = 1.1, cex.note = 0.8, cols = NULL, only.stat = FALSE,
  pdf.name = NULL, pdf.param = list(width = 12, height = 4.5),
  par.param = list(cex = 1, cex.main = 1, cex.sub = 1, cex.axis = 1))
}
\arguments{
\item{data}{input data frame. Rows are patients and columns are variables (e.g. demographics variables, time to event variables,
biomarker variables, treatment indicator, etc.). One patient per row.}

\item{outcome.class}{type of the outcome variable. Default is \code{c("survival", "binary", "continuous")}. Continuous is not available now}

\item{outcome.var}{name of the outcome varible. If the outcome.class is binary or coutinuous, only one value should be provided.
If the outcome.class is survival, two values should be provided - name of the 'time to event' variable and 'censorship' variable
 For the censoring variable, 1 indicates event and 0 indicates censoring. In all cases but when outcome.class=binary and rsp.cat=TRUE,
 patients with missing outcome variable (NA) will be excluded from BEP.}

\item{trt}{name of the treatment variable. If this is NULL, within-arm analysis will be performed}

\item{var}{name of the biomarker variable. only one variable should be specified.}

\item{var.class}{class of the variable. valid categories are "numeric", "categorical". If the class is continuous,
user needs to specify percentile.cutoff or numerical.cutoff to dichotomize the continuous measure into subgroups}

\item{var.name}{display name for the biomarker variable}

\item{percentile.cutoff}{percentile to dichotomize continuous biomarker measure. This could be a vector with multiple elements.
Values should be between 0 and 1}

\item{numerical.cutoff}{raw value to dichotomize continuous biomarker measure. numerical.cutoff and percentile.cutoff
cannot be both specified}

\item{greater}{whether calculate summary statistics within the subgroup whose biomarker value is greater than or equal to
cutoff value. If this is TRUE, in 2-arm study, across-arm HR within biomarker high group will be calculated.
In single arm study HR of biomarker high vs low will be calculated.}

\item{less}{whether calculate summary statistics within the subgroup whose biomarker value is less than the cutoff value.
greater and less can both be TRUE}

\item{greater.by.less}{whether show "greater" bin and "less" bin in consecutive rows. Default is FALSE. If it is TRUE,
parameters greater and less will be both set as TRUE}

\item{across.and.within}{whether show across- and within- arm results in the same figure. Default is FALSE. This parameter
will be ignored if number of arm is 1. If it is TRUE, within-arm analysis results will be shown below the across-arm results.}

\item{equal.in.high}{whether include equal in high group. Default is TRUE. If it is TRUE, ">=" and "<" will be
applied. Otherwise "<=" and ">" will be applied.}

\item{within.bin}{whether calculate summary statistics within bin (e.g. > cutoff1 and <= cutoff2). If within.bin is TRUE,
greater and less will be set as FALSE.}

\item{show.itt}{whether calculate summary statistics using all patients in full population (e.g. ITT). This will be ignored in 1arm case}

\item{show.bep}{whether calculate summary statistics using all patients in BEP (biomarker evaluable population). This will be ignored in 1arm case}

\item{bep}{name of the column which indicates biomarker evaluable population. If it is null, patients who have non NA records
in biomarker variable will be used as BEP.}

\item{bep.name}{preferred display name of the biomarker evaluable population.
If it is NULL, bep will be used.}

\item{itt.name}{preferred display name of the full population (e.g. ITT).
If it is NULL, "All" will be used.}

\item{bep.indicator}{In the subpopulation column, which value is used
to define the biomarker evaluable population.}

\item{covariate}{a vector specifying the covariate variables to be adjusted in the model. Default is set to NULL, meaning no adjustment.}

\item{strata}{name of the stratification variables. Default is set to NULL, meaning no stratification.}

\item{tabforest}{Default is FALSE. If it is FALSE, forest plot will be generated using forestplot::forestplot() function.
If it is TRUE, a table will be generated with forest plots incorpriated}

\item{quantile.type}{an integer between 1 and 9 selecting one of the nine quantile algorithms. See \code{\link{quantile}}. Default is 2.}

\item{cutoff.digits, digits}{cutoff.digits:number of digits for rounding when calculating cutoff. will only be used when percentile.cutoff is specified. digits: number of digits for the summary statistics display}

\item{placebo.code}{name of the control arm of the treatment variable. If you want to specify placebo code using this parameter, both placebo.code and active.code need to be provided.}

\item{active.code}{of the treatment/experimental arm of the treatment variable. If you want to specify active code using this parameter, both placebo.code and active.code need to be provided.}

\item{rsp.cat}{whether the response outcome variable is coded as binary (1 as responder and 0 as non-responder),
If rsp.cat is TRUE, responder categories
and nonresponder categories should be specified in rsp.response and rsp.nonresponse (all values in the outcome column
should be included in rsp.response and rsp.nonresponse)
. If rsp.cat is FALSE, the response outcome variable should be coded as binary (0/1).
At the same time rsp.response and rsp.nonresponse will be ignored.}

\item{rsp.response}{categories that should be considered as responder.}

\item{rsp.nonresponse}{categories that should be considered as non responder.}

\item{var.code}{ordered levels of the biomarker variable. This will be ignored for continuous biomarker.
If the biomarker is categorical and this is NULL, biomarker subgroups will be ordered by the order from factor() function}

\item{surv.conf.type}{confidence interval type. Default is "plain". see conf.type in survfit}

\item{ties}{Default is "efron". To match internal sas results, use "exact". See parameter "ties" in coxph.}

\item{alpha}{type I error rate. Default is 0.05.}

\item{main, main.prefix}{main title (prefix of title) of the forest plot. Default is "Association of biomarker effect within treatment arms".}

\item{sub}{sub title under the forest plot. Default is NULL.}

\item{clip}{range of the x-axis of the forest plot. Default is NULL.}

\item{xticks, xticks.digits}{x axis tick marks for the forest plot}

\item{xlab}{xlab for forest plot}

\item{cex.headings}{amount of magnification of headings of the forest plot relative to cex. Default is 1.1.}

\item{cex.note}{amount of magnification of the note. Default is 1.}

\item{cols}{Color of the 'effect size' displayed in the forest plot.}

\item{only.stat}{if it is TRUE, only summary statistics will be generated. No figure will be generated}

\item{pdf.name}{name of output pdf file. If it's NULL, the plots will be displayed but not saved as pdf. Default is "forestplot::forestplot.pdf".}

\item{pdf.param}{a list of parameters that define pdf graphics device. See \code{\link{pdf}}. Default is \code{list(width=6, height=4.5)}.}

\item{par.param}{a list of parameters that define graphcial parameters. See \code{\link{par}}. Default is \code{list(mar=c(4,4,3,2))}.}
}
\description{
This function creates a forest plot along with table with summary statistics to infer biomarker effects, within a single arm
or across two treatment arms.  The outcome could be survival, binary or continuous. This function can be used to summarize a single
biomarker variable
}
\examples{
data(input)
PlotTabForestBiomarker(data=input,
                      outcome.class=c("survival"),
                      outcome.var=c("PFS","PFS.event"),
                      trt="Arm",
                      var="KRAS.mutant",
                      var.class="categorical")
}
\author{
Ning Leng \email{leng.ning@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}, and previous team members (see DESCRIPTION)
}
