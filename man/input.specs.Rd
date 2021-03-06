% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.specs.R
\name{input.specs}
\alias{input.specs}
\title{Data specs example}
\format{A data frame with 18 rows and 2 variables:
\describe{
     \item{\code{Variable}}{Variable names.}
     \item{\code{Type}}{Describes how the column types must be specified:\cr\cr
                 \code{character:} used to label samples (e.g. Sample ID);\cr
                 \code{numeric:} continuous numeric quantities (e.g. Weight, Age);\cr
                 \code{categorical:} specifying categories of patients (e.g. Sex, CD8.ihc, Response);\cr
                 \code{time:} specifying time to event (e.g. PFS, OS);\cr
                 \code{event:} specifying if time is real or censored (1=actual, 0=censored).\cr
                 }
}}
\usage{
data(input.specs)
input.specs
}
\description{
Data specs define the type of a variable for each column.
}
\author{
Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}
}
