#' Show available documents
#'
#' Show available documents.
#'
#' @author Kwame Okrah \email{okrah.kwame@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}
#'
#' @export
AvailableDocs = function() {
    cat("=============================================================\n")
    cat("The following documents are currently available in gClinBiomarker:\n")
    cat("-------------------------------------------------------------\n")
    allowed.docs = c("1. gClinbiomarker_vig_v1", "2. gClinbiomarker_example_use_case")
    allowed.docs = paste(allowed.docs, collapse="\n")
    cat(allowed.docs)
    cat("\n")
    cat("-------------------------------------------------------------\n")
    cat("To open a document, call the function **OpenDoc(doc)**, where\n")
    cat("doc is any of the documents specified above;\n")
    cat('e.g. OpenDoc("gClinbiomarker_vig_v1")\n')
}

#' Access package documents
#'
#' This function allows the user the get access to pdf documents associated
#' with this package.
#'
#' @author Kwame Okrah, \email{okrah.kwame@gene.com}, Alexey Pronin \email{pronin.alexey@gene.com}
#'
#' @param doc a character vector specifying the document to open;
#' see \code{\link{AvailableDocs}}
#' @param verbose a logical vector indicating whether the path of the
#' document should be printed (Default = FALSE)
#'
#' @export
OpenDoc = function(doc, verbose = FALSE) {
    allowed.docs = c("gClinbiomarker_vig_v1", "gClinbiomarker_example_use_case")

    check.doc = doc %in% allowed.docs

    if (!check.doc) {
        m0 = "\n==================================="
        m1 = "\n*doc* must be one of the following:"
        m2 = "\n-----------------------------------\n"
        m3 = paste(allowed.docs, collapse="\n")
        msg = paste0(m0, m1, m2, m3)
        stop(msg, call. = FALSE)
    }

    if (doc == "gClinbiomarker_vig_v1") doc = "gClinbiomarker_vig_v1.pdf"
    if (doc == "gClinbiomarker_example_use_case") doc = "gClinbiomarker_example_use_case.pdf"

    f = system.file("doc", doc, package = "gClinBiomarker")

    if (.Platform$OS.type == "windows") {
        shell.exec(f)
    }else{
        system(paste(Sys.getenv("R_PDFVIEWER"), f, "&"))
    }

    if (verbose) {
        return(f)
    }
}
