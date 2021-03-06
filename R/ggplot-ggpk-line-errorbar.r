#' ggplot wrapper for plotting ribbon, taking a stats string or function set as
#' an argument
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @note When specifying arguments for wrapped ggplot geometries, prefix the
#'   parameter with geometry specifiers: "ribbons" for ribbon geometry "line"
#'   for central line geometry, or "label" for count labels. e.g.
#'   \code{ggpk_stat_line_errorbar(ribbons.alpha = 0.8, line.alpha = 0.6)}
#'
#' @param mapping a ggplot aesthetic mapping
#' @param data data to pass to ggplot2::ggplot2::stat_summary functions
#' @param show.counts True, False, "label" or "table" to display either as
#'   labels at each datapoint or as a table at the top of the plot
#' @param fun.data argument to be passed to \code{\link[ggplot2]{stat_summary}}
#'   calls
#' @param fun.args list of argumnets to be passed to
#'   \code{\link[ggplot2]{stat_summary}} calls
#' @param position PositionDodge object (constructed using
#'   \code{\link[ggplot2]{position_dodge}()}) to be used for positioning
#'   errorbars
#' @param id \code{ggpackets} id to indentify arguments which should be passed
#'   to all \code{ggpk_ribbons} ggplot2 layers. Defaults to
#'   \code{'ggpk_ribbons'}.
#' @param ...  Arguments prefixed with the following \code{id}s will be parsed
#'   and passed to the following calls:
#'   \itemize{
#'   \item \strong{ggpk_line_errorbar}, \strong{errorbar} ids: Passed to call
#'     to \code{\link[ggplot2]{geom_errorbar}}. \code{width} fixed.
#'   \item \strong{ggpk_line_errorbar}, \strong{errorbar} ids: Passed to call
#'     to \code{\link[ggplot2]{geom_errorbar}}. \code{alpha}, \code{linetype},
#'     \code{position}, \code{width}, \code{fun.data} and \code{size} fixed.
#'   \item \strong{ggpk_line_errorbar}, \strong{point} ids: Passed to call to
#'     \code{\link[ggplot2]{geom_point}}.
#'   \item \strong{ggpk_line_errorbar}, \strong{line} ids: Passed to call to
#'     \code{\link[ggplot2]{geom_line}}.
#'   \item \strong{ggpk_line_errorbar}, \strong{label} ids: Passed to call to
#'     \code{\link[ggrepel]{geom_label_repel}}. \code{direction},
#'     \code{nudge_y}, \code{label.size}, \code{fill} and \code{alpha} fixed.
#'   \item \strong{ggpk_line_errorbar}, \strong{label} ids: Passed to call to
#'     \code{\link[ggplot2]{geom_label}}. \code{label.size}, \code{fill} and
#'     \code{alpha} fixed.
#'   \item \strong{ggpk_line_errorbar}, \strong{label} ids: Passed to call to
#'     \code{\link{geom_text_table}}. \code{show.legend} fixed.
#'   }
#'
#' @return a ggpacket object that can be added to any ggplot to draw a
#'   collection of layers
#'
#' @examples
#' library(dplyr)
#' library(ggplot2)
#'
#' nasa %>%
#'   as.data.frame() %>%
#'   mutate(date = as.Date(sprintf("%04d%02d01", year, month), "%Y%m%d")) %>%
#'   mutate(hemisphere = ifelse(lat > 0,"Northern", "Southern")) %>%
#'   mutate(temperature = temperature - 273) %>%
#'
#'   # plotting
#'   ggplot() +
#'     aes(x=date, y=temperature) +
#'     ggpk_line_errorbar(fun.data = 'quartiles') +
#'     facet_grid(hemisphere ~ .) +
#'     labs(title="Temperature by Hemisphere")
#'
#' @export
#'
ggpk_line_errorbar <- function(mapping = NULL, data = NULL, show.counts = FALSE,
  fun.data = 'mean_se', fun.args = list(),
  position = ggplot2::position_dodge(width = ggplot2::rel(0.5)),
  id = 'ggpk_line_errorbar', ...) {

  fun.data <- stat_summary_funs(fun.data, fun.args)

  # prep dots arguments with default values, overwritten with ellipses args
  dots <- substitute(...())

  # add some conditional arguments to handle stat behavior
  if (!('identity' %in% dots[c(paste0(id, '.stat'), 'point.stat')]))
    dots$point.fun.data <- last(fun.data)

  if (!('identity' %in% dots[c(paste0(id, '.stat'), 'line.stat')]))
    dots$line.fun.data <- last(fun.data)

  if (!('identity' %in% dots[c(paste0(id, '.stat'), 'label.stat')]))
    dots$label.fun.data <- last(fun.data)

  ## pack ## ribbons
  # use standard errorbar geom if errorbar.stat is 'identity'
  (if ('identity' %in% dots[c(paste0(id, '.stat'), 'errorbar.stat')])
    ggpack(ggplot2::geom_errorbar, id = c(id, 'errorbar'),
           alpha = 0.5,
           size = ggplot2::rel(2),
           linetype = 1,
           position = position,
           dots = dots,
           width = 0)

    # reduce through list of errorbar geoms and collect sum
    else
      Reduce(`+`, mapply(function(f, i) {
        ggpack(ggplot2::geom_errorbar, id = c(id, 'errorbar'),
               stat = 'summary',
               dots = dots,
               alpha = 0.5,
               linetype = 1,
               position = position,
               width = 0,
               fun.data = f,
               size = dots$errorbar.size %||% ggplot2::rel(2) *
                 (1 - (i - 1) / length(fun.data)))

      }, f = fun.data, i = 1:length(fun.data))) ) +

    ## pack ## point
    # plot point along stat y
    ggpack(ggplot2::geom_point, id = c(id, 'point'),
           stat = 'summary',
           position = position,
           dots = dots) +

    ## pack ## line
    # plot line along stat y
    ggpack(ggplot2::geom_line, id = c(id, 'line'),
           stat = 'summary',
           position = position,
           dots = dots) +

    ## pack ## label
    # add labels of group counts
    if (isTRUE(show.counts) || show.counts == 'label') {
      tryCatch({
        geom_label_repel <- getExportedValue('ggrepel', 'geom_label_repel')
        ggpack(geom_label_repel,
               id = c(id, 'label'),
               stat = 'summary',
               dots = dots,
               direction = "y",
               nudge_y = 0.1,
               label.size = 0,
               fill = 'white',
               alpha = 0.85)
      }, error = function(e) {
        ggpack(ggplot2::geom_label, id = c(id, 'label'),
               stat = 'summary',
               dots = dots,
               label.size = 0,
               fill = 'white',
               alpha = 0.85)
      })

    } else if (show.counts == 'table') {
      ggpack(geom_text_table, id = c(id, 'label'),
             stat = 'summary',
             dots = dots,
             show.legend = FALSE)

    } else NULL

}
