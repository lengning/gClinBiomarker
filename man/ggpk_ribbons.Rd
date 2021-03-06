% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot-ggpk-ribbons.r
\name{ggpk_ribbons}
\alias{ggpk_ribbons}
\title{ggplot wrapper for plotting ribbon, taking a stats string or function set as
an argument}
\usage{
ggpk_ribbons(mapping = NULL, data = NULL, show.counts = FALSE,
  fun.data = "mean_se", fun.args = list(), id = "ggpk_ribbons", ...)
}
\arguments{
\item{mapping}{a ggplot aesthetic mapping}

\item{data}{data to pass to ggplot2::ggplot2::stat_summary functions}

\item{show.counts}{True, False, "label" or "table" to display either as
labels at each datapoint or as a table at the top of the plot}

\item{fun.data}{argument to be passed to \code{\link[ggplot2]{stat_summary}}
calls}

\item{fun.args}{list of argumnets to be passed to
\code{\link[ggplot2]{stat_summary}} calls}

\item{id}{\code{ggpackets} id to indentify arguments which should be passed
to all \code{ggpk_ribbons} ggplot2 layers. Defaults to
\code{'ggpk_ribbons'}.}

\item{...}{Arguments prefixed with the following \code{id}s will be parsed
and passed to the following calls:
\itemize{
\item \strong{ggpk_ribbons}, \strong{ribbon} ids: Passed to call to
  \code{\link[ggplot2]{geom_ribbon}}. \code{color} and \code{alpha} fixed.
\item \strong{ggpk_ribbons}, \strong{ribbon} ids: Passed to call to
  \code{\link[ggplot2]{geom_ribbon}}. \code{fun.data}, \code{color} and
  \code{alpha} fixed.
\item \strong{ggpk_ribbons}, \strong{point} ids: Passed to call to
  \code{\link[ggplot2]{geom_point}}.
\item \strong{ggpk_ribbons}, \strong{line} ids: Passed to call to
  \code{\link[ggplot2]{geom_line}}.
\item \strong{ggpk_ribbons}, \strong{label} ids: Passed to call to
  \code{\link[ggrepel]{geom_label_repel}}. \code{direction},
  \code{nudge_y}, \code{label.size}, \code{fill} and \code{alpha} fixed.
\item \strong{ggpk_ribbons}, \strong{label} ids: Passed to call to
  \code{\link[ggplot2]{geom_label}}. \code{label.size}, \code{fill} and
  \code{alpha} fixed.
\item \strong{ggpk_ribbons}, \strong{label} ids: Passed to call to
  \code{\link{geom_text_table}}. \code{show.legend} fixed.
}}
}
\value{
a ggpacket object that can be added to any ggplot to draw a
  collection of layers
}
\description{
ggplot wrapper for plotting ribbon, taking a stats string or function set as
an argument
}
\note{
When specifying arguments for wrapped ggplot geometries, id the
  parameter with geometry specifiers: "ribbon" for ribbon geometry "line" for
  central line geometry, or "label" for count labels. e.g.
  \code{ggpk_stat_ribbon(ribbon.alpha = 0.8, line.alpha = 0.6)}
}
\examples{
library(dplyr)
library(ggplot2)

nasa \%>\%
  as.data.frame() \%>\%
  mutate(date = as.Date(sprintf("\%04d\%02d01", year, month), "\%Y\%m\%d")) \%>\%
  mutate(hemisphere = ifelse(lat > 0,"Northern", "Southern")) \%>\%
  mutate(temperature = temperature - 273) \%>\%

  # plotting
  ggplot() +
    aes(x = date, y = temperature) +
    ggpk_ribbons(fun.data = "deciles") +
    facet_grid(hemisphere ~ .) +
    labs(title = "Temperature by Hemisphere")

}
\author{
Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
}
