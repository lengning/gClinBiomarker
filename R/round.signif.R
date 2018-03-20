#' @keywords internal

round.signif <- function(x, digits=NULL, p=NULL)
{
	  if(is.null(digits))digits <- p
  ifelse(abs(x)>=1, round(x, digits), signif(x, digits))
}
