#' Handling of argument passing to ggplots
#'
#' A helper function to wrap ggplot calls, allowing for passing of a shared list
#' of arguments prefixed by an id string. For example, a call to a custom ggplot
#' template may allow for arguments 'line.color' and 'bar.color' to specify the
#' line and bar colors seperately. These arguments are parsed appropriately and
#' passed to the appropriate sub-function.
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gene.com}
#'
#' @param _call the function to call
#' @param id the prefix string to subset arguments by
#' @param dots arguments to subset
#' @param ... additional arguments to use when calling the provided function
#' @param null.empty return NULL if no arguments are received
#' @return a call to the specified function with arguments subset for only those
#'   which match the specified prefix
#'
#' @note Use `substitute(...())` when collecting parent function ellipses args
#'   to pass into `ggpack()` to accommodate unevaluated symbols.
#' @note Ellipses arguments and the arguments passed through the `dots` args
#'   will overwrite previous arguments of the same name.
#'
#' @examples
#' library(ggplot2)
#'
#' # wrap layer calls in a ggpack to package them
#' ggpack(geom_bar(na.rm = TRUE))
#'
#' # or more explicitly define layer properties
#' ggpack(geom_bar,              # layer function to bundle
#'   id = 'bar',                 # name of layer - used for filter listed args
#'   fill = 'red',               # ...'s before 'args' used as defaults
#'   dots = list(fill = 'blue'), # named 'args' param will overwrite prior ...'s
#'   fill = 'green')             # ...'s after 'args' will overwrite any others
#'
#' @importFrom utils modifyList head tail
#' @export
ggpack <- function(`_call` = NULL, ..., id = NULL, dots = NULL, null.empty = FALSE) {
  if (is.null(`_call`) && length(list(...)) == 0) return(ggpacket())

  # unpack function call to split ... names into names before and after 'dots'
  ggpack.named.args <- names(formals())
  .arg.names <- Map(function(s) '['(s, !(s %in% ggpack.named.args)),
                    Map(names, split_on_name(as.list(sys.call()[-1])[-1], 'dots')))

  # split dots; 'default' (.dflt) & 'overwrite' (.over), left & right of 'dots'
  .dots   <- list(...)
  .dflt <- utils::head(.dots, length(.arg.names$left))  # default
  .dflt <- .dflt[!duplicated(names(.dflt), fromLast = TRUE)]
  .over <- utils::tail(.dots, length(.arg.names$right)) # overwrite
  .over <- .over[!duplicated(names(.over), fromLast = TRUE)]

  fdots <- filter_args(id, dots)
  thru <- c(fdots[names(fdots) == ''], Reduce(utils::modifyList,
                                              list(.dflt, fdots[names(fdots) != ''], .over)))

  if (null.empty && length(thru) == 0) return(ggpacket(NULL))

  callfname <- deparse(as.list(match.call())$'_call')
  # call ggproto construction with stripped args to determine geom
  # might cause issue if aesthetics are dependent on additional arguments
  ggproto_tmp <- tryCatch({
      do.call(`_call`, thru[names(thru) %in% c('geom', 'stat')])
    }, error = function(e) NULL)
  geom <- if ('ggproto' %in% class(ggproto_tmp)) ggproto_tmp$geom else NULL
  stat <- if ('ggproto' %in% class(ggproto_tmp)) ggproto_tmp$stat else NULL

  # account for mismatched aesthetics when mapping being passed through
  # (code largely made to mirror ggplot2 layer.r's compute_aesthetics())
  if (!is.null(geom)) {
    if ('mapping' %in% names(thru)) {
      thru$mapping <- filter_aesthetics(geom, thru$mapping)
    } else if (any(names(thru) %in% allowed_aesthetics(geom))) {
      thru_aes_idx <- names(thru) %in% allowed_aesthetics(geom) &
        sapply(thru, aes_arg_is_uneval)
      aes_string_args <- lapply(thru[thru_aes_idx], deparse)
      thru <- thru[!thru_aes_idx]
      thru$mapping <- do.call(ggplot2::aes_string, aes_string_args)
    }
  }

  ggpacket_name <- if (is.null(id)) NULL else paste(id, collapse = '.')
  if (all(class(`_call`) == 'function'))
    ggpacket(do.call(`_call`, thru), ggpacket_name)
  else
    ggpacket(`_call`, ggpacket_name)
}



#' Helper function to split a list on a named args
#'
#' Specifically used for splitting a call's arguments based on the location
#' of a specified argument name to use the placement of the arguments to
#' represent useful information.
#'
#' @param l list to split
#' @param split_name name of element in list to split on. If the element is not
#' found in the list, the returned list will contain the entirety of the passed
#' list, l, in the name 'left' and an empty list in the name 'right'.
#' @param omit.named logical indicating whether the list value specified by
#' split_name should be included in the output lists.
#'
#' @return A named list with two named values, "left" and "right", indicating
#' arguments to the left and right of the split_name argument.
#'
split_on_name <- function(l, split_name, omit.named = TRUE) {
  idx <- match(split_name, names(l))
  if (is.na(idx)) return(list(left = l, right = list()))
  list(left = l[1:(idx - omit.named)],
       right = l[(idx + omit.named):length(l)])
}


#' Helper function for ggpack to filter arguments based on a prefix
#'
#' @param prefix a string specifying the prefix to be used to filter args
#' @param args a list of arguments to be subset by the prefix match
#' @param sep a regex joining string between prefix and args.
#' Defaults to \code{"\\\\."}.
#'
#' @return a list of arguments that originally were prefaced by the
#' specified prefix, now with that prefix removed.
#'
#' @importFrom stats setNames
filter_args <- function(prefix, args, sep = '\\.') {
  if (is.null(prefix) || is.null(args)) return(args %||% list())
  unnamed_args <- unlist(args[which(names(args) %in% prefix)], use.names = FALSE)
  named_args <- args[grep(paste0('^',prefix,'.', collapse = '|'), names(args))]
  named_args <- stats::setNames(named_args,
                                # e.g. ^(prefix1|prefix2|prefix3)\\.(.*)$
                                gsub(paste0('^(',paste0(prefix, collapse='|'),')',sep,'(.*)$'),
                                     '\\2', names(named_args)))
  as.list(c(named_args, unnamed_args))
}

#' A class for wrapping ggplot layers.
#'
#' ggplot ggproto objects can not be added together as they
#' are and must be made reusable by first adding to a list.
#' This class aggregates ggproto objects into a list and handles
#' adding into a ggplot construction.
#'
#' @author Doug Kelkhoff \email{kelkhoff.douglas@gnee.com}
#'
#' @examples
#' library(ggplot2)
#'
#' # a ggpacket can be created as an object by not passing arguments
#' # to the ggpack function. Any ggproto layers or ggpack'd ggpacket
#' # layers can be added directly
#' bar_error_counts <- ggpacket() +
#'   stat_summary(fun.y = mean, geom = 'bar') +
#'   stat_summary(fun.data = mean_se, width = 0.2, geom = "errorbar") +
#'   stat_summary(fun.data = function(d) c(
#'                  y = mean(d) + sd(d)/sqrt(length(d)),
#'                  label = length(d)),
#'                vjust = -1, geom = 'text')
#'
#' ggplot(mtcars, aes(x = gear, y = mpg)) +
#'   bar_error_counts
#'
#' # easier functionalization of subcomponents of a ggplot
#' custom_errorbars <- function(error_function, ...) {
#'   ggpack(stat_summary, prefix = 'bar', args = substitute(...()),
#'          fun.y = mean, geom = 'bar') +
#'   ggpack(stat_summary, prefix = 'errorbar', args = substitute(...()),
#'          fun.data = error_function, width = 0.2, geom = "errorbar")
#' }
#'
#' ggplot(mtcars, aes(x = gear, y = mpg)) +
#'   custom_errorbars(mean_se, bar.fill = gear) +
#'   facet_grid(. ~ am, scales = 'free')
#'
#'
#' @slot ggcalls a list of ggproto layer objects
#'
#' @import methods
#' @export ggpacket
ggpacket <- setClass("ggpacket",
  slots = c(ggcalls = "list"),
  prototype = list(ggcalls = list())
)

#' Initialize empty ggpacket
#' @param .Object ggpacket object to be initialized
#' @importFrom methods signature
setMethod(f = "initialize", methods::signature(.Object = "ggpacket"),
  function(.Object) { .Object })

#' Initialize new object by adding ggproto to list with name label
#' @param ggproto_obj an object of class ggproto to add to the ggpacket's internal list
#' @param label an optional label to be able to index this component within the ggpacket
setMethod("initialize", "ggpacket", function(.Object, ggproto_obj = NULL, label = NULL) {
  if (is.null(ggproto_obj)) return(.Object)
  else if (all(class(ggproto_obj) == 'list'))
    .Object <- .Object + stats::setNames(
      ggproto_obj,
      names(ggproto_obj) %||% label)
  else
    .Object <- .Object + stats::setNames(
      list(ggproto_obj),
      label)

  .Object
})

#' Overload show method to print ggpacket
#' @param object the ggpacket object to show
setMethod("show", "ggpacket", function(object) {
  ## attempt to build plot by calling packet with `ggplot() + ggpacket_obj`
  plt_output <- try(silent = TRUE, {
    plt <- ggplot2::ggplot() + object
    utils::capture.output(type = 'message', {
      bld <- ggplot2::ggplot_build(plt)
    })
  })

  if (length(plt_output) == 0 && all(Map(length, bld$data) > 0))
    return(show(plt))

  cat("ggpacket\nA container for multiple ggplot ggproto objects\n\n")
  cat('Errors preventing standalone plotting: \n')
  if (length(plt_output))
    message(attr(plt_output, 'condition')$message)
  else if (!all(Map(length, bld$data) > 0))
    message('Not all layers have been passed sufficient data')

  cat('\n')
  if (length(object@ggcalls) == 0) cat("empty\n\n")
  else
    mapply(
      function(n, name, ggp) cat(sprintf("[[%s]] %s\n%s\n\n", n, name, ggp)),
      n = 1:length(object@ggcalls),
      ggp = Map(function(ggc)
        paste(capture.output(print(ggc)), collapse = "\n", sep = "\n"),
        object@ggcalls),
      name = names(object@ggcalls) %||% ""
    )
})

#' Overload names method to print ggpacket ggcall list names
#' @param x the ggpacket object to show
setMethod("names", "ggpacket", function(x) { names(x@ggcalls) })

#' Overload [ generic to index into ggcalls list
#' @param x ggpacket object
#' @param i index of the ggcall in the ggpacket object
setMethod("[", "ggpacket", function(x, i) { x@ggcalls[i] })

#' Overload [ generic to index into ggcalls list
#' @param x ggpacket object
#' @param i index of the ggcall in the ggpacket object
setMethod("[[", "ggpacket", function(x, i) { x@ggcalls[[i]] })

#' Primitive methods for adding ggpackets to various ggplot objects
#' @param e1 ggpacket object to collect additional terms
#' @param e2 anything which would be added to a ggplot::gg call for evaluation
setMethod("+", c("ggpacket", "ANY"), function(e1, e2) {
  if      ("ggpacket" %in% class(e2)) e1@ggcalls <- append(e1@ggcalls, e2@ggcalls)
  else if ("theme" %in% class(e2))    e1@ggcalls <- append(e1@ggcalls, list(e2))
  else                                e1@ggcalls <- append(e1@ggcalls, e2)
  e1
})

#' Add NULL to a ggpacket to return the ggpacket
#' @param e1 NULL value
#' @param e2 the ggpacket object
#' @return when NULL is added to a ggpacket object, the original object is returned
setMethod("+", c("NULL", "ggpacket"), function(e1, e2) e2)

# register unexported gg class from ggplot2 so signature is accepted
setOldClass("gg")

#' Combine ggpacket contents with gg call to incorporate seamlessly with ggplot calls
#' @param e1 ggplot gg class object to add to ggpacket internal list
#' @param e2 ggpacket object
setMethod("+", c("gg", "ggpacket"), function(e1, e2) Reduce("+", e2@ggcalls, e1))

#' Get list of allowed aesthetics for a given geometry class
#'
#' @param geom a ggplot2 Geom object
#'
#' @return a character vector of the names of aesthetics permitted by this
#' geometry object (including Americanized names and base R equivalent names
#' also accepted by ggplot2)
#'
allowed_aesthetics <- function(geom) {
  aes_names <- c('x', 'y', 'group', geom$required_aes, names(geom$default_aes))
  add_eqv_aes(aes_names)
}


#' Add equivalent Americanized and base equivalent names to ggplot aesthetic
#' list
#'
#' @param aes_names a character vector of aesthetic names
#'
#' @return a character vector of aesthetic names including any Americanized or
#' base R equivalent argument names accepted by ggplot2.
#'
#' @import ggplot2
add_eqv_aes <- function(aes_names) {
  base_eqv_idx <- unlist(safe_private_export('ggplot2', '.base_to_ggplot')) %in% aes_names |
    names(safe_private_export('ggplot2', '.base_to_ggplot'))  %in% aes_names
  base_eqv <- safe_private_export('ggplot2', '.base_to_ggplot')[base_eqv_idx]
  aes_names_new <- unique(c(aes_names, names(base_eqv), unlist(base_eqv, use.names = F)))
  # attempt to add new equivalent mappings if anything new added (e.g. (color => colour) => fg)
  if (!(all(aes_names_new %in% aes_names))) add_eqv_aes(aes_names_new)
  else aes_names_new
}


#' Helper function to filter aesthetic mappings based on geometry
#'
#' @param geom a ggplot2 Geom object
#' @param mapping a ggplot2 aesthetic mapping
#'
#' @return the mapping filtered by accepted aesthetics for the given Geom
#'
#' @import ggplot2
#' @export
filter_aesthetics <- function(geom, mapping) {
  allowed_aes <- allowed_aesthetics(geom)
  mapping_aes_names <- names(safe_private_export('ggplot2', 'rename_aes')(mapping))
  disallowed_aes <- setdiff(mapping_aes_names, allowed_aes)
  do.call(remove_aesthetics, c(list(mapping), disallowed_aes))
}


#' Helper function to flatten specific aesthetics to group
#'
#' @param mapping a ggplot2 aesthetic mapping
#' @param ... aesthetic handles to flatten into the group aesthetic as
#' interaciton terms
#'
#' @return the mapping aesthetics with specified aesthetics filtered out and
#' instead stored as interaction terms in the group aesthetics with any
#' preexisting group aesthetics
#'
#' @export
flatten_aesthetics_to_group <- function(mapping, ...) {
  .dots <- add_eqv_aes(unlist(list(...)))
  mapped_vars <- mapping[!(names(mapping) %in% c('x', 'y'))]
  if (length(.dots) > 0)
    mapped_vars <- mapped_vars[names(mapped_vars) %in% c(list('group'), .dots)]
  mapped_vals <- unique(unlist(mapped_vars, use.names=FALSE))
  if (length(mapped_vals) > 0)
    mapping$group <- as.call(c(list(as.symbol("interaction")), mapped_vals))
  mapping
}


#' Helper function to cull aesthetic mapping from parameter inputs in ellipses
#'
#' @param ... ellipses containing aesthetics mapping
#' @param geom if specified, a ggplot2 Geom object can be passed to filter
#' aesthetics for only a specified Geom.
#'
#' @return a list containing two fields: 'aes' which contains the aesthetic
#' mappings within the arguments passed and 'not_aes' which contains a list of
#' all other arguments
#'
#' @export
split_aes_from_dots <- function(..., geom = NULL) {
  if (is.null(geom)) aes_list <- safe_private_export('ggplot2', '.all_aesthetics')
  else aes_list <- allowed_aesthetics(geom)
  aes_args     <- substitute(...())
  not_aes_args <- aes_args[!names(aes_args) %in% aes_list]
  aes_args     <- structure(aes_args[names(aes_args) %in% aes_list], class = 'uneval')
  list(aes = do.call(ggplot2::aes, aes_args),
       not_aes = not_aes_args)
}


#' Helper function to filter out aesthetics from mappings
#'
#' @param mapping aesthetic mapping to use for filtering
#' @param ... mapping labels to filter out
#'
#' @return an aesthetic mapping with filtered mappings removed
#' and group mapping set as interaction terms of all non axial
#' terms.
#'
#' @export
remove_aesthetics <- function(mapping, ...) {
  mapping <- flatten_aesthetics_to_group(mapping, ...)
  mapping[!(names(mapping) %in% list(...))]
}


#' Determine if an argument is an unevaluated expression
#' @param i the input object to be tested
#'
#' @return TRUE if the value is of type 'name', 'call' or 'expression'
aes_arg_is_uneval <- function(i) { any(is.name(i), is.call(i), is.expression(i)) }

#' Wrapper for common decorators to package
#'
#' @param ... arguments to be passed to any of the bundled
#' decorator plot functions, blank, xlab, ylab, labs or theme
#' prefaced by any of those strings ('labs.title').
#'
#' @return a ggpacket object containing a collection of
#' decorating plot functions.
#'
#' @export
ggpk_decorators <- function(...) {
  dots <- substitute(...())
  ggpack(ggplot2::xlab,   id = 'xlab',    dots = dots, null.empty = TRUE) +
  ggpack(ggplot2::ylab,   id = 'ylab',    dots = dots, null.empty = TRUE) +
  ggpack(ggplot2::labs,   id = 'labs',    dots = dots, null.empty = TRUE) +
  ggpack(ggplot2::theme,  id = 'theme',   dots = dots, null.empty = TRUE) +
  ggpack(ggplot2::ggtitle,id = 'ggtitle', dots = dots, null.empty = TRUE)
}



