% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sas-convenience-funcs.R
\name{clean_emmeans_specs}
\alias{clean_emmeans_specs}
\title{Helper function to clean one-sided formulas for emmeans}
\usage{
clean_emmeans_specs(specs)
}
\arguments{
\item{specs}{the specs input to process}
}
\value{
a two-sided formula with the LHS being 'pairwise' if the input is a
  formula, otherwise return original input
}
\description{
The \code{emmeans::emmeans()} function broke a golden rule: they
  return a different class of data based on unrelated parameters. In this
  case, if the \code{specs} parameter is a one-sided formula it returns a
  \code{emmGrid} whereas if it's two-sided for a list, it returns a
  \code{emm_list}. This makes it quite hard to work with reliably. To address
  this, this helper function cleans the inputs before handing them off to
  \code{specs} to make a one sided formula (like \code{~ x + y}) into an
  acceptable two-sided formula (like \code{pairwise ~ x + y}), allowing
  emmeans to be used and always producing a single class of output.
}
