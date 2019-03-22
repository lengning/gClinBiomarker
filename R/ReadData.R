#' Read data
#'
#' Read two csv files and perform concordance check. 
#' First file contains data.
#' Second file provides specs for data.
#'
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}
#' 
#' @importFrom utils read.csv
#'
#' @param input A csv file with data.
#' @param input.specs A csv file with specs for \code{input}. Default is NULL.
#'
#' @return A list with two data frames: data and data.specs.
#' 
#' @note The ReadData() function reads in the data file and spec file from csv files. 
#' The function also perform concordance check between the data and its spec. 
#' Data column names and spec entry names should match. No duplicate names are allowed.
#' If input.specs is NULL, it will be created based on the input data and then returned without checks.
#' 
#' @examples
#' \dontrun{
#' ReadData("input.csv", "input_specs.csv")
#' ReadData("input.csv")
#' }
#'
#' @export

ReadData <- function(input, input.specs=NULL) {
    
    # Check the input file exists.
    if (file.exists(input)) {
        input <- read.csv(input, stringsAsFactors=FALSE)
    } else {
        stop("The input file does not exist!")
    }
    
    # Check for duplicate columns in the input file.
    dup.input <- abs(length(colnames(input)) - length(unique(colnames(input))))
    if (dup.input != 0) {
        stop("There are duplicate name columns in the input file!")
    }

    # If input.specs is not NULL. Do checks.
    if (!is.null(input.specs)) {
        # Check the spec file exists.
        if (file.exists(input.specs)) {
        input.specs <- read.csv(input.specs, stringsAsFactors=FALSE)
        } else {
            stop("The specs file does not exist!")
        }
        # Check data types in specs are allowed.
        types.allowed <- c("character", "numeric", "categorical", "time", "event")
        if (!(all(input.specs$Type %in% types.allowed))) {
            stop(paste("Allowed data types are:",
                       "character, numeric, categorical, time, and event.",
                       "Please check the types in the specs file!"))
        }
        # Check for duplicate columns in the spec file.
        dup.input.specs <- abs(length(input.specs$Variable) - length(unique(input.specs$Variable)))
        if (dup.input.specs != 0) {
            stop("There are duplicate rows in the specs file!")
        }
        # Check colnames(input) match input.specs$Variable.
        col.dif <- length(setdiff(colnames(input), input.specs$Variable))
        if (col.dif != 0) {
            stop("The column names between input and input.specs files do not match!")
        }
        # Check class of variables.
        for (i in 1:nrow(input.specs)) {
            # Convert categorical to factor.
            if (input.specs$Type[i] == "categorical") {
                input[, i] <- as.factor(input[, input.specs$Variable[i]])
            }
            # Check the class of numeric variables is numeric.
            if (input.specs$Type[i] %in% c("numeric", "time", "event")) {
                if (is.numeric(input[, input.specs$Variable[i]]) == FALSE) {
                    stop("The variable ", input.specs$Variable[i],
                         " should be numeric but identified as ",
                         class(input.specs$Variable[i]), "!")
                }
            }
        }
    } else { #If input.specs is NULL, create it.
        Variable <- names(input)
        Type <- sapply(input, mode)    
        input.specs <- data.frame(Variable, Type, row.names=NULL)
    }
  
    return(list(data=input, data.specs=input.specs))
}