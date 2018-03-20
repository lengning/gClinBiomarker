#' Compute emmeans to mirror SAS emmeans functionality
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}, Christina Rabe
#'   \email{rabe.christina@gene.com}
#'
#' @param model.obj A model object on which to calculate emmeans
#' @param specs specs to be passed to emmeans()
#' @param data specified data to be used for model fitting. (defaults to data
#'   bound to model, though this can not always been inferred from the model
#'   object and scope of function call)
#' @param mode mode to be passed to emmeans()
#' @param quietly whether to hide message informing of changes to default
#'   functionality
#' @param verbose whether all debug information should be printed to console
#' @param confidence.level the confidence interval to use for upper and lower
#'   confidence level outputs
#' @param ... additional arguments to be passed to emmeans
#'
#' @importFrom emmeans emmeans
#' @importFrom broom tidy
#' @import dplyr
#'
#' @examples
#' # create model on which to calculate emmeans
#' # for sake of reusable example lm is used, but can accommodate all models
#' model <- lm(mpg ~ hp + carb + wt,
#'   data = mtcars, na.action = na.exclude)
#'
#' # call emmeans, specifying factors over which means should be calculated
#' sas.emm <- sas.emmeans(model, ~ carb, data = mtcars, confidence.level = 0.9)
#'
#' # update the confidence level used
#' # NOTE: if you're examining confidence intervals, it's much faster to udpate
#' # rather than re-run the model and means
#' sas.emm <- update(sas.emm, level = 0.95)
#'
#' # using the means as a dataframe
#' as.data.frame(sas.emm)
#'
#' @export
#' @importFrom dplyr "%>%" filter_at vars all_vars mutate_at
sas.emmeans <- function(model.obj, specs, data = NULL, mode = 'kenward-roger',
    quietly = FALSE, verbose = FALSE, confidence.level = NULL, ...) {

  if (!quietly) message(paste(
    "To mirror SAS functionality:",
    " 1. mode of 'kenward-roger' will be used unless otherwise specified.",
    " 2. data used for initial model fit...",
    "   a. is filtered for only complete cases in the response variable",
    "   b. will have any numeric variables in specs converted to categorical",
    " 3. grouping variables for emmeans converted to factor.",
    "      (creating a new factor column where necessary)",
    " 4. continuous means calculated using rows where response is not NA.",
    " 5. categorical level weights from all rows of dataset\n",
    sep = "\n"))

  envir <- parent.frame()
  if (is.call(q <- as.list(match.call(expand.dots = TRUE)[-1])$model.obj)) {
    model.call <- as.list(match.call(eval(q[[1]], envir), q, expand.dots = TRUE))
    model.call[-1] <- lapply(model.call[-1], eval, envir)
    model      <- append_split_terms(as.formula(get_model_formula(model.call)))
    model.env  <- envir
  } else {
    model.call <- as.list(model.obj$call)
    model      <- append_split_terms(as.formula(get_model_formula(model.obj)))
    model.env  <- attr(model$terms, ".Environment")
  }

  # prevent model that doesn't store model construction without explicit data
  if (is.null(data) && identical(model.env, environment()))
    stop(paste(
      "\nModel environment was not found within model object. ",
      "To resolve, provide value for data to sas.emmeans() as well.", sep = "\n"))

  # filter dataset down to complete cases, coerce spec vars to factor
  if (is.null(data)) data <- eval(model.call$data, envir=model.env)

  specs         <- clean_emmeans_specs(specs) # ~ x | y => pairwise ~ x | y
  specs.pred    <- get_emmeans_specs_predictors(specs) # c('x', 'y')
  specs.pred.nf <- Map(class, Filter(Negate(is.factor), data[specs.pred]))

  # convert non-factor variables to factor
  cleaned_data <- data %>%
    dplyr::filter_at(dplyr::vars(model$vars), dplyr::all_vars(!is.na(.))) %>%
    dplyr::mutate_at(as.character(names(specs.pred.nf)), as.factor)

  if (!quietly & nrow(data) != nrow(cleaned_data)) message(sprintf(
    "%d records with incomplete data removed in model fit. (%1.2f%%)\n",
     nrow(data) - nrow(cleaned_data),
    (nrow(data) - nrow(cleaned_data)) / nrow(data) * 100))

  if (verbose) message(sprintf(
    '%s model is being refit with arguments: \n%s\n',
    as.character(model.call[[1]]),
    print_to_string(model.call[-1], condense = TRUE)))

  model.call$data <- cleaned_data
  sas.model.obj <- do.call(as.character(model.call[[1]]), model.call[-1])

  non.spec.vars <- setdiff(model$covs, specs.pred)
  emmeans_args <- modifyList(
    list(object = sas.model.obj, specs = specs, mode = mode, data = cleaned_data,
         at = cleaned_data %>% numeric_means(non.spec.vars, verbose = verbose)
         ) %>% append_weights(data = data, verbose = verbose),
    list(...))

  if (verbose) {
    message('Running emmeans::emmeans() with parameters: ')
    message(print_to_string(list(object = summary(emmeans_args$object))))
    message(print_to_string(emmeans_args[setdiff(names(emmeans_args), 'object')]))
  }

  emmlist <- do.call(emmeans::emmeans, emmeans_args)
  emmlist$emmeans@grid <-
    coerce_from_factor(emmlist$emmeans@grid, specs.pred.nf,
      suffix = '.emmeans.factor', verbose = verbose) %>%
    copy_variable_attributes(data)
  emmlist$emmeans@misc$unary.vars.gdf <- summarize_unary_vars(data, specs.pred)
  if (!is.null(confidence.level)) emmlist <- update(emmlist, level = confidence.level)
  emmlist
}


#' List means over numeric variables
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @param data A dataframe upon which means will be calculated
#' @param potential_vars A character vector of potential variable names to
#'   consider
#' @param verbose include verbose output of numeric means calculations
#'
#' @return A named list with names as the numeric subset of the variable names
#'   specified in potential_vars and values of the means
#'
#' @importFrom dplyr summarize_at funs "%>%"
numeric_means <- function(data, potential_vars = names(data), verbose = FALSE) {
  numeric_vars <- names(Filter(is.numeric, data[potential_vars]))

  if (verbose && length(sd <- setdiff(potential_vars, numeric_vars)) > 0)
    message('While calculating "at" means, some variables omitted due to ',
            'non-numeric class: \n', paste(sd, collapse = ", "), '\n')

  data %>%
    dplyr::summarize_at(as.character(numeric_vars), dplyr::funs(mean)) %>%
    as.list
}


#' Append argument list with calculated weights as frequency
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @description emmeans expects accepts a weights argument to specify the
#'   denominator for mean calculations. This function will calculate these
#'   weights from the frequency of the emmeans levels from a new dataset.
#'
#'   Arguments with which emmeans is to be called are passed to this function
#'   which will call \code{emmeanss()} with \code{weights = 'show.levels'} in
#'   order to determine the levels over which weights are calculated. It will
#'   then use the passed \code{data} to evaluate frequencies for those weights.
#'
#' @param emmeans_args arguments with which emmeans will be called
#' @param data a dataset over which the weights should be calculated as level
#'   frequencies.
#' @param verbose whether additional information should be printed to console
#'
#' @return The same arguments passed through \code{emmeans_args} with the
#'   weights parameter specified by the frequencies of occurence in the
#'   \code{data} variable, or the \code{emmeans_args$data} list item if no data
#'   was included (default behavior of \code{emmeans}).
#'
#' @importFrom emmeans emmeans
#' @importFrom utils modifyList
#' @importFrom dplyr left_join mutate_at group_by summarize ungroup pull "%>%"
append_weights <- function(emmeans_args, data = NULL, verbose = FALSE) {
  levels_args <- utils::modifyList(emmeans_args, list(weights = 'show.levels'))
  levels <- sink_to_temp(do.call(emmeans::emmeans, levels_args))

  # even with 'show.levels', when no groups are defined emmeans returns emm_list
  if (any(class(levels) %in% c('emm_list', 'emmGrid')) || is.null(names(levels)))
    return(utils::modifyList(emmeans_args, list(weights = 'proportional')))

  if (verbose) message(sprintf('Calculating weights over levels: %s',
    paste(names(levels), collapse = ', ')))

  level_weights <- dplyr::left_join(levels,
      (data %||% emmeans_args$data) %>%
        dplyr::mutate_at(vars(names(levels)), as.factor) %>%
        dplyr::group_by(.dots=names(levels)) %>%
        dplyr::summarize(wgt=n()) %>% dplyr::ungroup(),
      by = names(levels))

  if (verbose) message(print_to_string(level_weights), '\n')

  emmeans_args$weights <- level_weights %>% dplyr::pull(wgt)
  emmeans_args
}


#' Coerce dataframe variables to new type
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @description Coerce many variables to new classes based on a named list of
#'   variable names with values of the character string of their target class.
#'
#' @param data the dataframe whose variables should be coerced to a new class
#' @param class_list a named list of classes to be coerced to from class factor,
#'   with the list names being the variable names (e.g. \code{list(vs =
#'   'character', carb = 'numeric')}). Classes are all coerced using the
#'   function \code{as.<class>(levels(<var>)[<var>])}.
#' @param preserve whether original factor variables should be preserved
#' @param suffix a suffix to add to preserved variable names
#' @param verbose whether additional information should be printed to console
#'
#' @return the original dataframe with variables coerced to the target class and
#'   optionally preserved copies of the original variables.
#'
#' @importFrom dplyr mutate_ "%>%"
#' @importFrom stats setNames
coerce_from_factor <- function(data, class_list, preserve = TRUE,
    suffix = '.old', verbose = FALSE) {

  if (length(class_list) == 0) return(data)

  if (verbose) message(
    'Coercing variables according to: ', print_to_string(class_list), '\n')

  data %>%
    dplyr::mutate_(.dots =
      if (!preserve) list()
      else stats::setNames(names(class_list), paste0(names(class_list), suffix))
    ) %>%
    dplyr::mutate_(.dots =
      Map(function(n, c) sprintf("as.%s(levels(%s)[%s])", c, n, n),
      names(class_list), class_list))
}


#' Get model formula from a model object
#'
#' @description Search within a model object for model formula information and
#'   derive the formula If the model object contains a named element matching
#'   any name in the \code{params} parameter, that item will be returned,
#'   presumably containing the environment information in which it was
#'   instantiated.
#'
#'   If the model object doesn't store this information, the model call will be
#'   inspected to search for the formula argument used and derive formula from a
#'   new formula produced from this input. As a last resort, it will take the
#'   first passed argument and try to coerce it to a formula. A formula derived
#'   in this way will be created in the local function environment and will not
#'   be able to be used to evaluate the symbols stored within them.
#'
#' @param model.obj the model object with which to derive the formula
#' @param params model attributes or parameters to consder for deriving formula
#'
#' @return A formula object or item which may be coerced to a formula
#'
get_model_formula <- function(model.obj, params = c('terms', 'model', 'formula')) {
  # search in the model object itself (WILL contain model env)
  mdl.ind <- first(na.omit(c(NA, match(params, names(model.obj[-1])))))
  if(!is.na(mdl.ind)) return(model.obj[[mdl.ind+1]])

  # pull from model call (will NOT contain model env)
  model.call  <- as.list(model.obj$call)
  if (length(model.call) <= 1)
    stop('Model object does not contain any arguments. Cannot discern model formula')
  mdl.ind <- first(na.omit(c(NA, match(params, names(model.call[-1])))))
  if (!is.na(mdl.ind)) return(model.call[[mdl.ind+1]])

  # try to cast the first argument to formula (will NOT contain model env)
  tryCatch({ as.formula(model.call[[2]]); model.call[[2]] },
    error = function(e) {
      stop("Attempted to coerce first argument to formula but failed.")
    })
}


#' Wrap a formula in a helper list splitting response and predictor terms
#'
#' @param t the formula or \code{terms()} object to construct from
#'
#' @return a list with four elements, "terms" containing formula terms, "resp"
#'   containing a vector of reponse variable names, "covs" containing predictor
#'   covariate variable names and "vars" containing both response and covariate
#'   variable names.
#'
#' @importFrom stats terms
append_split_terms <- function(t) {
  m <- list(terms = stats::terms(t))
  if (attr(m$terms, 'response')) {
    m$resp <- all.vars(m$terms[[2]])
    m$covs <- all.vars(m$terms[[3]])
  } else m$covs <- all.vars(m$terms[[2]])
  m$vars <- c(m$resp, m$covs)
  m
}


#' Helper function to clean one-sided formulas for emmeans
#'
#' @description The \code{emmeans::emmeans()} function broke a golden rule: they
#'   return a different class of data based on unrelated parameters. In this
#'   case, if the \code{specs} parameter is a one-sided formula it returns a
#'   \code{emmGrid} whereas if it's two-sided for a list, it returns a
#'   \code{emm_list}. This makes it quite hard to work with reliably. To address
#'   this, this helper function cleans the inputs before handing them off to
#'   \code{specs} to make a one sided formula (like \code{~ x + y}) into an
#'   acceptable two-sided formula (like \code{pairwise ~ x + y}), allowing
#'   emmeans to be used and always producing a single class of output.
#'
#' @param specs the specs input to process
#'
#' @return a two-sided formula with the LHS being 'pairwise' if the input is a
#'   formula, otherwise return original input
#'
clean_emmeans_specs <- function(specs) {
  if (inherits(specs, 'formula') && !attr(terms(specs), 'response'))
    as.formula(paste('pairwise ~', as.character(specs)[[2]]))
  else specs
}


#' Return predictors from specs input
#'
#' @description \code{emmeans::emmeans()} accepts a bunch of different variable
#'   types to its \code{specs} parameter. This function handles that variety of
#'   variable types to figure out what variables were specified.
#'
#' @param specs a specs parameter to be passed to \code{emmeans::emmeans()}
#'
#' @return a character vector of the variables from the specs object
#'
#' @importFrom stats delete.response terms
get_emmeans_specs_predictors <- function(specs) {
  if (inherits(specs, "formula"))
    all.vars(stats::delete.response(stats::terms(specs)))
  else if (class(specs) == 'character') specs
  else if (class(specs) == 'list')
    stop('Passing list as specs is currently unsupported in sas.emmeans. Please file an issue if this is functionality you require.')
  else stop('Unable to parse variables from provided specs')
}


#' @method as.data.frame emm_list
#' @export
#' @importFrom dplyr left_join "%>%"
as.data.frame.emm_list <- function(x, row.names, optional, ...) {
  if ('emeans' %in% names(x) && 'unary.vars.gdf' %in% names(x$emmeans@misc)) {
    dplyr::left_join(
      as.data.frame(x$emmeans),
      x$emmeans@misc$unary.vars.gdf %>% dplyr::ungroup(),
      by = attr(x$emmeans@misc$unary.vars.gdf, 'vars'),
      suffix = c('.emmeans', ''))
  } else
    as.data.frame(x$emmeans)
}


#' @method as.data.frame emmGrid
#' @export
#' @importFrom dplyr left_join select "%>%"
as.data.frame.emmGrid <- function(x, row.names, optional, ...) {
  as.data.frame(summary(x)) %>%
    dplyr::left_join(
      x@grid %>% dplyr::select(c(names(x@levels), '.wgt.')),
      by = names(x@levels))
}


#' @method update emm_list
#' @export
#' @import emmeans
update.emm_list <- function(object, ..., silent = FALSE) {
  object$emmeans <- safe_private_export('emmeans', 'update.emmGrid')(object$emmeans, ...)
  object
}


#' @method tidy emm_list
#' @export
#' @importFrom broom tidy
tidy.emm_list <- function(emm, ...) { broom::tidy(emm@emmeans) }
