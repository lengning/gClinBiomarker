#' Read a VAD
#'
#' Read a VAD and create specs for it. 
#'
#' @author Alexey Pronin \email{pronin.alexey@gene.com}, Ning Leng \email{leng.ning@gene.com}
#' 
#' @importFrom haven read_sas
#' @importFrom utils write.csv
#'
#' @param file A SAS dataset.
#' @param write.csv.to Defines where both data and data.specs csv files will be written to. 
#' The csv files will be written in the form name_date.csv and name_specs_date.csv.
#' Default is NULL.
#'
#' @return A list with two data frames: data and data.specs.
#' 
#' @note The ReadVAD() function reads in a VAD and creates data specs. 
#' 
#' @examples
#' \dontrun{
#' ReadVAD("ars.sas7bdat")  
#' ReadVAD("ars.sas7bdat", ".")  
#' ReadVAD("ars.sas7bdat", "Desktop/R")
#' ReadVAD("Desktop/SAS/ars.sas7bdat", "Desktop/R")
#' }
#'
#' @export

ReadVAD <- function(file, write.csv.to=NULL) {
    
    # Check the input file exists.
    if (file.exists(file)) {
        data <- as.data.frame(read_sas(file))
    } else {
        stop("The VAD does not exist!")
    }
    
    # Parse SAS file and create specs.
    Label <- sapply(data, attr, "label")
    Type <- sapply(data, mode)    
    Variable <- names(Type)
    data.specs <- data.frame(Variable, Type, Label, row.names=NULL)
    
    # Write data and data.specs to csv.
    if (!is.null(write.csv.to)) {
        filename <- strsplit(basename(file), "[.]")[[1]][1]
        data_filepath <- paste(write.csv.to, "/", filename, "_", Sys.Date(), ".csv", sep="")
        data_specs_filepath <- paste(write.csv.to, "/", filename, "_specs_", Sys.Date(), ".csv", sep="")
        write.csv(data, data_filepath, row.names=FALSE)
        write.csv(data.specs, data_specs_filepath, row.names=FALSE)
    }

    return(list(data=data, data.specs=data.specs))
}