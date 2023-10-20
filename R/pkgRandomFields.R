#'
#'           pkgRandomFields.R
#' 
#'    Dealing with the (DEFUNCT) Random Fields package
#'
#'    $Revision: 1.6 $  $Date: 2023/10/20 11:17:19 $

kraeverRandomFields <- function() {
  stop("The package RandomFields is no longer available.")
##  kraever("RandomFieldsUtils")
##  kraever("RandomFields")
# should no longer be needed:  
#  capture.output(RandomFieldsUtils:::.onLoad())
#  capture.output(RandomFields:::.onLoad())
  return(invisible(NULL))
}

# require a namespace and optionally check whether it is attached
kraever <- function(package, fatal=TRUE) {
  if(!requireNamespace(package, quietly=TRUE)) {
    if(fatal)
      stop(paste("The package", sQuote(package), "is required"),
           call.=FALSE)
    return(FALSE)
  }
  if(spatstat.options(paste("check", package, "loaded", sep=".")) &&
    !isNamespaceLoaded(package)){
    if(fatal)
      stop(paste("The package", sQuote(package),
                 "must be loaded: please type",
                 sQuote(paste0("library", paren(package)))),
           call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

# legacy functions

RandomFieldsSafe <- function() { FALSE }

getRandomFieldsModelGen <- function(model) {
  return(NULL)
}


