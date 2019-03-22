#' Wrapper for fitting model and adjusting to predicted values and passing newly
#' fitted model to ggplot geom_stat_ribbon
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @param data data to use for model fitting and plotting
#' @param mapping ggplot aesthetic mapping to use for plotting
#' @param model model function to use for fitting; defaults to lm
#' @param model.per grouping variables by which to isolate individual models
#' @param model.formula model formula to use for model fitting
#' @param facet.fun function to use for ggplot faceting in ggplot2::facet_grid
#' @param plot.style one of 'errorbars' or 'ribbons'
#' @param ... Additional arguments allow for an assortion of additional
#'   functionality.
#'   \itemize{
#'   \item{"Aesthetics"}{
#'   First, arguments which can be used as any of ggplot's
#'   default aesthetics will be pulled from \code{...} args.
#'   }.
#'   \item{"ggplot layers"}{
#'   Arguments are passed to underlying ggplot components by prefixing value
#'   with ggplot call name. ("ribbons", "line", "text", "facet", "xlab", "ylab",
#'   "labs" or "theme" - e.g. \code{ribbons.color = 'red')}).
#'   }
#'   \item{"Model fitting"}{
#'   Arguments with prefixes of "model" will be passed to the model function if
#'   one is provided.
#'   }
#'   \item{"ggpackets"}{
#'   Additional arguments passed to
#'   \code{\link[gClinBiomarker]{ggpk_ribbons}} or
#'   \code{\link[gClinBiomarker]{ggpk_line_errorbar}}.
#'   }
#'   }
#'
#' @return a ggplot object
#'
#' @examples
#' # load data
#' nasa.data <- as.data.frame(dplyr::nasa)
#'
#' # default representation
#' PlotLong(nasa.data, x = month, y = temperature)
#'
#' # a more appropriate representation for large n
#' PlotLong(nasa.data, x = month, y = temperature, fun.data = 'quartiles')
#'
#' # include linear adjustment accounting for amount of ozone
#' PlotLong(nasa.data, x = month, y = temperature,
#'          formula = temperature ~ ozone, fun.data = 'deciles',
#'          show.counts = T)
#'
#' # adjusting by independent models for the northern and southern hemispheres
#' library(dplyr) # needed for %>%
#' PlotLong(nasa.data %>% mutate(hemi=ifelse(lat>0, "North", "South")),
#'          x = month, y = temperature, formula = temperature ~ ozone,
#'          model.per = ~ hemi, facet.fun = ~ hemi, fun.data = 'deciles',
#'          xlab = "Month", ylab = "Temperature Adjusted for Ozone",
#'          labs.title = "Temperature by Hemisphere",
#'          labs.caption = "*idependent models fit per hemisphere")
#'
#' # including a table of value counts and subsetting value data to specific
#' # months
#' library(dplyr) # needed for %>%
#' PlotLong(nasa.data %>% mutate(hemisphere=ifelse(lat > 0, "North", "South")),
#'          x = month, y = temperature, group = hemisphere,
#'          color = hemisphere, fill = hemisphere,
#'          formula = temperature ~ ozone,
#'          model.per = ~ hemisphere, fun.data = 'deciles',
#'          plot.style = 'errorbars',
#'          show.counts = 'table',
#'          label.data = . %>% filter(month %in% c(1, 6, 12)),
#'          label.hjust = 'inward',
#'          xlab = "Month", ylab = "Temperature Adjusted for Ozone",
#'          labs.title = "Temperature by Hemisphere",
#'          labs.caption = "*idependent models fit per hemisphere")
#'
#' @export
#'
#' @import ggplot2 dplyr
PlotLong <- function(data, mapping = NULL, model = lm, model.per = NULL,
    model.formula = NULL, facet.fun = NULL, plot.style = 'ribbons', ...) {

  if (is.null(mapping)) {
    args <- split_aes_from_dots(...)
    mapping <- args$aes
    .dots   <- args$not_aes
  } else .dots <- list(...)

  # collapse linetype to group to allow for ggpack overrides (errorbar color)
  mapping <- flatten_aesthetics_to_group(mapping, 'linetype')
  if (all(c('ymin', 'ymax') %in% names(mapping)))
    .dots <- modifyList(.dots, list(plotlong.stat = 'identity'))

  # accommodate model fitting
  if (!is.null(model.formula))  {
    # ensure model variables reflect plotted variables
    if (gsub("^~", "", deparse(mapping$y)) != model.formula[[2]])
      stop('Independent model.formula variable must be the same as y aesthetic')

    # add model fit output to data
    data <- do.call(augment_predict,
        c(list(data, model, model.per, model.formula=model.formula), .dots))

    # change aesthetic y to be fitted values
    mapping$y <- as.name(".fitted")
    .dots$ylab <- .dots$ylab %||% paste("Adjusted", deparse(model.formula[[2]]))
  }

  # choose our plotting functions
  if (plot.style == 'ribbons') ggcall.plot <- ggpk_ribbons
  else if (plot.style == 'errorbars') ggcall.plot <- ggpk_line_errorbar

  # choose facetting functions
  if (is.null(facet.fun)) ggcall.facets <- ggplot2::facet_null()
  else
    ggcall.facets <- ggpack(ggplot2::facet_grid, id = 'facet',
        dots = .dots, facets = facet.fun)

  # plot using geom_stat_ribbons, passing extra arguments to geom
  data %>% ggplot2::ggplot() + mapping +
    do.call(ggcall.plot, c(.dots, list(id = 'plotlong'))) +
    ggcall.facets +
    do.call(ggpk_decorators, .dots)
}
