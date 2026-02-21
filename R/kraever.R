#'
#'           kraever.R
#' 
#'    Requiring a package namespace
#'
#'    $Revision: 1.9 $  $Date: 2026/02/21 03:20:59 $

# require a namespace and optionally check whether it is attached
kraever <- function(package, fatal=TRUE, loaded=TRUE, ..., preamble=NULL) {
  if(!requireNamespace(package, quietly=TRUE)) {
    if(fatal) {
      p <- paste("package", sQuote(package))
      q <- if(is.null(preamble)) paste("The", p) else paste(preamble, "the", p)
      stop(paste(q, "is required"), call.=FALSE)
    }      
    return(FALSE)
  }
  if(loaded && !isNamespaceLoaded(package)){
    if(fatal) {
      p <- paste("package", sQuote(package))
      q <- if(is.null(preamble)) paste("The", p) else paste(preamble, "the", p)
      stop(paste(q, "must be loaded: please type",
                 sQuote(paste0("library", paren(package)))),
           call.=FALSE)
    }
    return(FALSE)
  }
  return(TRUE)
}

