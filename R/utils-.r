`%||%` <- function (a, b) {
  if (!is.null(a)) a
  else b
}

`%>%` <- dplyr::`%>%`

#' A function for safely trying to get private namespace functions
#'
#' Predominately, this is used to get private exports from ggplot2 for ggpackets
#' functions. This behavior is more responsibly handled in ggpackets
#' (https://www.github.com/dgkf/ggpackets). Once published, gClinBiomarker will
#' instead include ggpackets as a dependency and this behavior will not be
#' needed.
#'
#' @param pkg package name as character
#' @param f package function name as character
#'
safe_private_export <- function(pkg, f) {
  if (require(pkg, quietly = TRUE, character.only = TRUE) &&
      f %in% names(getNamespace(pkg))) {
    getNamespace(pkg)[[f]]
  } else
    stop(sprintf('Invalid namespace value `%s:::%s` requested', pkg, f))
}
