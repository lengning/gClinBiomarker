#' Helper function to hide side effect output to console
#'
#' @description Some packages make it very difficult to hide console output of
#'   their functions. In situations where you want to use those functions
#'   without cluttering your console output, this function can be used to wrap
#'   code for evaluation, suppressing all outputs.
#'
#' @param expr the expression to be evaluated
#' @return the result of the expression.
#'
sink_to_temp <- function(expr) {
  sink(tempfile()); on.exit(suppressWarnings(sink()))
  invisible(eval.parent(expr))
}


#' Capture output and return as string
#'
#' @description Take output from a segment of code and return it as a single
#'   string.
#'
#' @param expr the expression to be evaluated
#' @param condense whether empty lines should be stripped from the output
#'
#' @return a single string of the captured output
#'
#' @importFrom utils capture.output
print_to_string <- function(expr, condense = FALSE) {
  out <- utils::capture.output(print(eval.parent(expr)))
  if (condense) out <- Filter(nzchar, out)
  paste(out, collapse = '\n')
}
