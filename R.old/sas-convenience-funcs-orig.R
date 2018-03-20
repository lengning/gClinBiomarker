#' Compute LSMEANS to mirror SAS functionality
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}, Christina Rabe \email{rabe.christina@gene.com}
#'
#' @param model.obj A model object on which to calculate lsmeans
#' @param specs specs to be passed to lsmeans()
#' @param mode mode to be passed to lsmeans()
#' @param as.data.frame whether to return a lsm.list object or a data.frame
#' @param quietly whether to hide message informing of changes to default functionality
#' @param ... additional arguments to be passed to lsmeans
#'
#' @importFrom lsmeans lsmeans
#' @importFrom broom tidy
#' @import dplyr
#'
#' @examples
#' library(dplyr) # for %>% operator
#'
#' # create model on which to calculate lsmeans
#' gls(mpg ~ hp + carb + wt,
#'     data = mtcars %>%
#'              mutate(carb = as.factor(carb)), # make sure contrast is factor
#'     na.action = na.exclude) %>%
#'
#' # call lsmeans, specifying factors over which means should be calculated
#' sas.lsmeans(~ carb)
#'
sas.lsmeans.orig <- function(model.obj, specs, data = NULL, mode = 'kenward-roger',
                        quietly = FALSE, verbose = FALSE, ...) {

  if (!quietly) message(paste(
    "To mirror SAS functionality:",
    " - mode of 'kenward-roger' will be used unless otherwise specified.",
    " - model.obj will be fit with only complete response variables.",
    " - grouping variables for lsmeans converted to factor.",
    "     (creating a new factor column where necessary)",
    " - continuous means calculated using rows where response is not NA.",
    " - categorical level weights from all rows of dataset\n",
    sep = "\n"))

  # get model.obj call
  model.call  <- as.list(model.obj$call)

  ## Defensive handling of inputs
  # try to get model formula from model object, defaulting to first arg
  # if a 'model' or 'formula' argument are not found
  if ('model' %in% names(model.call)) model.form <- as.formula(model.call$model)
  else if ('formula' %in% names(model.call)) model.form <- as.formula(model.call$formula)
  else if (length(model.call) > 1) model.form <- as.formula(model.call[[2]])
  else stop('Model formula parameter not found in model object')
  model.form.trms <- terms(model.form)

  # collect variables defined from the specs input
  if (inherits(specs, "formula")) {
    specs.dep.vars <- all.vars(delete.response(terms(specs)))
    # edit specs to preface with 'pairwise' if not provided to lsmeans
    # specifying a response makes the function type stable, always returning a lsm.list
    if (!attr(terms(specs), 'response'))
      specs <- as.formula(paste0('pairwise ~ ', as.character(specs)[[2]]))
  } else if (class(specs) == 'character') specs.dep.vars <- specs
  else if (class(specs) == 'list') specs.dep.vars <- unlist(specs, use.names = F)
  else stop('Unable to parse variables from provided specs')

  # ensure model has a response variable
  if (length(model.form) < 3)
    stop('model must be a formula of the form "resp ~ pred"')


  ## argument handling
  # modify original model call data to use subsetted data, extracting full
  # dataset for use in lsmeans
  model.call.data <- model.call$data
  if (is.null(data))
    data <- eval(substitute(model.call.data), envir = sys.frame(-1))
  specs.dep.nf <- Map(class, Filter(Negate(is.factor), data[specs.dep.vars]))

  # prep model formula dissection terms for modifying lsmeans call
  resp <- all.vars(model.form[[2]])
  covs <- all.vars(model.form[[3]])

  # convert any numeric grouping variables to factor
  data <- data[which(complete.cases(data[resp])),] %>%
    dplyr::mutate_at(.vars = names(specs.dep.nf), .funs = as.factor)

  # refit model with new, complete cases data for reponse variables & factor
  # grouping variables
  model.call$data <- data

  if (!quietly & verbose) message(sprintf(
    '%s model is being refit with arguments: \n%s',
    as.character(model.call[[1]]),
    show(model.call[-1])))

  sas.model.obj <- do.call(as.character(model.call[[1]]), model.call[-1])

  # evaluate variable classes to determine lsmeans behavior
  model.form.vars <- as.character(attr(model.form.trms, "variables")[-1])
  model.form.type <- Map(class, data[model.form.vars])
  non.spec.type <- model.form.type[!(names(model.form.type) %in% c(specs.dep.vars, resp))]

  # split non-grouping model covariates into numeric and non-numeric lists
  covs.f <- Filter(function(t) t != 'numeric', non.spec.type)
  covs.n <- Filter(function(t) t == 'numeric', non.spec.type)

  # give a bit of info about how data was subset
  orig.data.len <- nrow(data)
  new.data.len <- nrow(eval(model.call$data, envir = sys.frame(-1)))
  if (!quietly & orig.data.len != new.data.len) message(paste0(
    orig.data.len - new.data.len,
    " records with incomplete data removed in model fit. (",
    round((orig.data.len - new.data.len) / orig.data.len * 100, digits=2),
    "%)\n"))

  ## lsmeans calculation
  # prep lsmeans levels for weighting lsmeans call (sink cat outputs to null)
  .dots <- modifyList(
    list(object = sas.model.obj,
         specs = specs,
         mode = mode,
         # mean continuous variables calculated with proportional weights and filtered dataset
         # TODO: check if complete cases should only be for reponse variable or also for numeric covariates
         # TODO: if also for numeric covariates, remove na.rm=T in summarise_all
         at = data %>%
           dplyr::do( .[which(complete.cases(.[resp])),] ) %>%
           dplyr::select(names(covs.n)) %>%
           dplyr::summarise_all(dplyr::funs(mean(.))) %>%
           as.list,
         # weights for categorical terms using marginal frequencies
         weights = "show.levels"),
    list(...))

  sink(tempfile()); on.exit(suppressWarnings(sink()))
  lsmeans.levels <- do.call(lsmeans::lsmeans, .dots)
  sink()

  if (!quietly & verbose) {
    message('Levels for lsmeans weights parameter order:')
    message(show(lsmeans.levels))
  }

  # prep sas-style arguments with which to call lsmeans
  if (is.null(names(lsmeans.levels)) < 1 || length(covs.f) < 1)
    .dots$weights = 'proportional'
  else
    .dots$weights = dplyr::left_join(
      lsmeans.levels,
      broom::tidy(table(data[names(covs.f)],
                        dnn =names(covs.f))),
      by = names(lsmeans.levels)) %>%
    dplyr::pull(Freq)

  if (!quietly & verbose) {
    message('Running lsmeans::lsmeans() with parameters: ')
    message(show(list(object = summary(.dots$object))))
    message(show(.dots[!which(names(.dots) %in% 'object')]))
  }

  lsmlist <- do.call(lsmeans::lsmeans, .dots)

  # modify grid to revert back to numeric values that had to be coerced to
  # factors (retaining factorized form used for lsmeans)
  lsmlist$lsmeans@grid <- lsmlist$lsmeans@grid %>%
    dplyr::mutate_(.dots = setNames(
      names(specs.dep.nf),
      paste(names(specs.dep.nf), 'lsmeans.factor', sep = '.'))) %>%
    dplyr::mutate_(.dots =
      Map(function(n, c) sprintf("as.%s(levels(%s)[%s])", c, n, n),
      names(specs.dep.nf), specs.dep.nf))

  # add some extra information to retain unary values from aggregated original dataframe
  lsmlist$unary.vars.gdf <- data %>%
    dplyr::mutate_(.dots =
      Map(function(n, c) sprintf("as.%s(levels(%s)[%s])", c, n, n),
      names(specs.dep.nf), specs.dep.nf)) %>%
    dplyr::group_by_(.dots = specs.dep.vars) %>%
    (function(gdf) {
      dplyr::left_join(
        gdf %>% dplyr::select_(.dots = specs.dep.vars) %>%
          dplyr::slice(1),
        gdf %>% dplyr::select_(.dots = c(specs.dep.vars,
           dplyr::summarise_all(., funs(n_distinct)) %>%
           dplyr::select_if(funs(all(. == 1))) %>%
           names)) %>%
          dplyr::slice(1),
        by = specs.dep.vars)
    })

  lsmlist
}

#' @method as.data.frame lsm.list
#' @export
as.data.frame.lsm.list <- function(lsm) {
  if ('unary.vars.gdf' %in% names(lsm)) {
    dplyr::left_join(
      as.data.frame(lsm$lsmeans),
      lsm$unary.vars.gdf %>% dplyr::ungroup(),
      by = attr(lsm$unary.vars.gdf, 'vars'))
  } else
    as.data.frame(lsm$lsmeans)
}

#' @method as.data.frame lsmobj
#' @export
as.data.frame.lsmobj <- function(lsmo) {
  as.data.frame(summary(lsmo)) %>%
    left_join(lsmo@grid %>% select(c(names(lsmo@levels), '.wgt.')),
              by = names(lsmo@levels))
}

#' @method tidy lsm.list
#' @export
tidy.lsm.list <- function(lsm) { broom::tidy(lsm@lsmeans) }
