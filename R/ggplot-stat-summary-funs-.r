#' helper function to interpret strings and functions passed to
#' geom_stat_ribbons
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @param format either a function returning a named list of stat
#'   transformations or a string matching to one of 'mean_se', 'median_iqr',
#'   'tukey_hinges', 'tukey_whiskers' or 'quantile_contour', or a list of
#'   strings from this group. The y calculation of the first function in the
#'   list returned is used for placing labels.
#' @param args additional arguments to pass to target function
#'
#' @return a list of functions to apply to pass to fun.data in stat_summary
#'
#' @export
#'
stat_summary_funs <- function(format, args) {
  if (is.character(format)) {
    unlist(lapply(format, function(f) switch(f,

    mean_sd = function(d) c(
      y = mean(d),
      ymin = mean(d) - sd(d),
      ymax = mean(d) + sd(d),
      label = length(d)),

    mean_se = function(d) c(
      y = mean(d),
      ymin = mean(d) - sd(d)/sqrt(length(d)),
      ymax = mean(d) + sd(d)/sqrt(length(d)),
      label = length(d)),

    median_iqr = function(d) c(
      y = median(d),
      ymin = as.numeric(quantile(d, 0.25)),
      ymax = as.numeric(quantile(d, 0.75)),
      label = length(d)),

    tukey_hinges = function(d) c(
      y = median(d),
      ymin = as.numeric(quantile(d, 0.25)),
      ymax = as.numeric(quantile(d, 0.75)),
      label = length(d)),

    tukey_notches = function(d) c(
      y = median(d),
      ymin = as.numeric(quantile(d, 0.25)) - 1.58 * IQR(d) / sqrt(length(d)),
      ymax = as.numeric(quantile(d, 0.75)) + 1.58 * IQR(d) / sqrt(length(d)),
      label = length(d)),

    tukey_whiskers = function(d) c(
      y = median(d),
      ymin = max(min(d), as.numeric(quantile(d, 0.25)) - 1.5 * IQR(d)),
      ymax = min(max(d), as.numeric(quantile(d, 0.75)) + 1.5 * IQR(d)),
      label = length(d)),

    tukey = list(stat_summary_funs('tukey_hinges'),
                stat_summary_funs('tukey_whiskers')),

    deciles = stat_summary_iles(seq(0.1, 0.5, 0.1)),

    quartiles = stat_summary_iles(c(0.25, 0.5)),

    quantiles = do.call(stat_summary_iles, args)

    )))
  } else if (is.function(format)) {
    if (is.null(args))
      list(format)
    else
      list(do.call(format, args))
  } else if (is.list(format) && is.function(format[[1]])) format
  else stop("fun.data must be either a string, function or list of functions")
}

stat_summary_iles <- function(s) {
  lapply(s, function(s) {
    function(d) c(
      y = median(d),
      ymin = as.numeric(quantile(d, 0.5 - s)),
      ymax = as.numeric(quantile(d, 0.5 + s)),
      label = length(d)
    ) })
}

stat_summary_from_seq <- function(s) {
  mapply(function(q_l, q_u) {
    function(d) c(
      y = median(d),
      ymin = as.numeric(quantile(d, q_l)),
      ymax = as.numeric(quantile(d, q_u)),
      label = length(d)
    ) },
    q_l = head(s, length(s) %/% 2),
    q_u = rev(tail(s, length(s) %/% 2)) )
}
