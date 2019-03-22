#' Convert color-names or RGB-code to possibly semi-transparent RGB-code.
#' 
#' Function takes the name of a color and converts it into the rgb space. Parameter "alpha" allows
#' to specify the transparency within [0,1], 0 meaning completey transparent and 1 meaning completey
#' opaque. If an RGB-code is provided and alpha != 1, the RGB-code of the transparency adapted color 
#' will be returned.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}
#' 
#' @param col (character) name of the color to be converted/transformed into RGB-space (code). Only
#'               those colors can be used which are part of the set returned by function colors(). Defaults
#'               to "black".
#' @param alpha (numeric) value specifying the transparency to be used, 0 = completely transparent, 
#'               1 = opaque.
#' 
#' @return RGB-code
#' 
#' @examples 
#' # convert character string representing a color to RGB-code using alpha-channel of .25 (75\% transparent)
#'      as.rgb("red", alpha=.25)
#' 
#' # same thing now using the RGB-code of red (alpha=1, i.e. as.rgb("red"))
#'      as.rgb("#FF0000FF", alpha=.25)
#'      
#' @export

as.rgb <- function(col="black", alpha=1)
{
    if(length(col) > 1 && (length(alpha) == 1 || length(alpha) < length(col)))         # unclear which alpha to use or only one alpha specified
    {
        if(length(alpha) < length(col) && length(alpha) > 1)
            warning("Multiple (but too few) 'alpha' specified! Only use 'alpha[1]' for each color!")
        return(sapply(col, as.rgb, alpha=alpha[1]))
    }
    if(length(col) > 1 && length(col) <= length(alpha))                                 # process each color separately
    {
        res <- character()
        for(i in 1:length(col))
            res <- c(res, as.rgb(col[i], alpha[i]))
        return(res)
    }
    if( col %in% colors() )
        return( rgb(t(col2rgb(col))/255, alpha=alpha) )
    else
    {
        col <- sub("#", "", col)
        R <- as.numeric(paste("0x", substr(col, 1,2), sep=""))
        G <- as.numeric(paste("0x", substr(col, 3,4), sep=""))
        B <- as.numeric(paste("0x", substr(col, 5,6), sep=""))
        return( rgb(R/255, G/255, B/255, alpha=alpha, maxColorValue=1) )
    }        
}