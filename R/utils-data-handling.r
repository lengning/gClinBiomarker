#' Adds prediction outputs to data based on model and formula
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @param .data data to use for model fitting
#' @param model model function to use (e.g `lm`)
#' @param model.per a formula representing all terms to group data by to fit
#'   individual models (e.g. \code{~ am + vs})
#' @param ... additional arguments to be passed to the model function, prefixed
#'   with 'model' (eg. \code{model.singular.ok = FALSE})
#'
#' @return a data frame with new columns based on model fit
#'
#' @importFrom broom augment
#' @importFrom dplyr rename rename_ group_by group_by_
#'
#' @export
#' @importFrom dplyr group_by_ rename
#' @importFrom broom augment
#' @importFrom stats setNames
augment_predict <- function(.data, model, model.per = NULL, ...) {
  .dots <- filter_args('model', list(...))
  ## Temporarily commenting out warning suppression in broom::augment
  ## Likely can be deleted permenantly after a quick test
  # withCallingHandlers({
    .data %>%
      dplyr::group_by_(.dots = model.per) %>%
      do(broom::augment(do.call(model, c(list(data=.), .dots)), .))
  # }, warning = function(w) {
  #   if (grep('Deprecated.*purrr::possibly()', conditionMessage(w)))
  #     invokeRestart('muffleWarning')
  # })
}


#' Copy attributes of dataframe variables from a source dataframe
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @description particularly useful for copying variable label attributes, which
#'   are often not preserved during data transformation operations.
#'
#' @param obj target dataframe
#' @param attr_src_obj source dataframe
#'
#' @return the target dataframe with attributes copied from the source dataframe
#'   for each variable.
#'
copy_variable_attributes <- function(obj, attr_src_obj) {
  for (n in intersect(names(obj), names(attr_src_obj)))
    attributes(obj[[n]]) <- attributes(attr_src_obj[[n]])
  obj
}


#' Produce summary table of unary value of specified levels
#'
#' @description Given a dataframe and set of grouping variables, a dataset will
#'   be returned with all groups summarized containing only variables with one
#'   unique value per group, with a value of that single unique value.
#'
#' @param data the dataframe to summarized
#' @param ... character vectors or strings containing variable names to group
#'   over
#' @param .dots a character vector ccontaining variables names to group over
#'
#' @return a dataframe with a single value per group combination with only
#'   variables with a unique value per group.
#'
#' @importFrom dplyr group_by select summarize_all select_if slice funs "%>%"
summarize_unary_vars <- function(data, ..., .dots = c()) {
  groups <- c(unlist(list(...)), .dots)
  data %>%
    dplyr::group_by(.dots = groups) %>%
    dplyr::select(c(
      groups,
      dplyr::summarize_all(., dplyr::funs(n_distinct)) %>%
        dplyr::select_if(dplyr::funs(all(. == 1))) %>%
        names
    )) %>% dplyr::slice(1)
}
