#'
#'           kraever.R
#' 
#'    Requiring a package namespace
#'
#'    $Revision: 1.8 $  $Date: 2025/12/01 15:16:01 $

# require a namespace and optionally check whether it is attached
kraever <- function(package, fatal=TRUE, loaded=TRUE) {
  if(!requireNamespace(package, quietly=TRUE)) {
    if(fatal)
      stop(paste("The package", sQuote(package), "is required"),
           call.=FALSE)
    return(FALSE)
  }
  if(loaded && !isNamespaceLoaded(package)){
    if(fatal)
      stop(paste("The package", sQuote(package),
                 "must be loaded: please type",
                 sQuote(paste0("library", paren(package)))),
           call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

