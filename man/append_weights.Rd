% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sas-convenience-funcs.R
\name{append_weights}
\alias{append_weights}
\title{Append argument list with calculated weights as frequency}
\usage{
append_weights(emmeans_args, data = NULL, verbose = FALSE)
}
\arguments{
\item{emmeans_args}{arguments with which emmeans will be called}

\item{data}{a dataset over which the weights should be calculated as level
frequencies.}

\item{verbose}{whether additional information should be printed to console}
}
\value{
The same arguments passed through \code{emmeans_args} with the
  weights parameter specified by the frequencies of occurence in the
  \code{data} variable, or the \code{emmeans_args$data} list item if no data
  was included (default behavior of \code{emmeans}).
}
\description{
emmeans expects accepts a weights argument to specify the
  denominator for mean calculations. This function will calculate these
  weights from the frequency of the emmeans levels from a new dataset.

  Arguments with which emmeans is to be called are passed to this function
  which will call \code{emmeanss()} with \code{weights = 'show.levels'} in
  order to determine the levels over which weights are calculated. It will
  then use the passed \code{data} to evaluate frequencies for those weights.
}
\author{
Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
}
