## clusterfunctions.R
##
## Contains the generic functions:
##  - clusterkernel
##  - clusterfield
##  - clusterradius.
##
##   $Revision: 1.11 $  $Date: 2022/02/21 02:24:34 $
##

clusterkernel <- function(model, ...) {
  UseMethod("clusterkernel")
}
clusterkernel.character <- function(model, ...){
  info <- spatstatClusterModelInfo(model, onlyPCP = TRUE)
  internalkernel <- info$kernel
  dots <- list(...)
  par <- c(kappa = 1, scale = dots$scale)
  par <- info$checkpar(par, native = TRUE)
  nam <- info$shapenames
  margs <- NULL
  if(!is.null(nam))
    margs <- dots[nam]
  f <- function(x, y = 0, ...){
    internalkernel(par = par, rvals = sqrt(x^2+y^2), margs = margs)
  }
  return(f)
}

## The method clusterkernel.kppm is in spatstat.core

clusterfield <- function(model, locations = NULL, ...) {
    UseMethod("clusterfield")
}

clusterfield.character <- function(model, locations = NULL, ...){
    f <- clusterkernel(model, ...)
    clusterfield.function(f, locations, ...)
}

clusterfield.function <- function(model, locations = NULL, ..., mu = NULL) {
  if(is.null(locations)){
    locations <- ppp(.5, .5, window=square(1))
  } else if(!is.ppp(locations))
    stop("Argument ", sQuote("locations"), " must be a point pattern (ppp).")

  if("sigma" %in% names(list(...)) && "sigma" %in% names(formals(model)))
    warning("Currently ", sQuote("sigma"),
            "cannot be passed as an extra argument to the kernel function. ",
            "Please redefine the kernel function to use another argument name.")

  if(requireNamespace("spatstat.core")) {
    rslt <- spatstat.core::density.ppp(locations, kernel=model, ..., edge=FALSE)
  } else {
    message("The package spatstat.core is required.")
    return(NULL)
  }
  
  if(is.null(mu))
    return(rslt)
  
  mu <- as.im(mu, W=rslt)
  if(min(mu)<0)
    stop("Cluster reference intensity ", sQuote("mu"), " is negative.")

  return(rslt*mu)
}


## The method clusterfield.kppm is in spatstat.core

clusterradius <- function(model, ...){
    UseMethod("clusterradius")
}

clusterradius.character <- function(model, ..., thresh = NULL, precision = FALSE){
  info <- spatstatClusterModelInfo(model, onlyPCP=FALSE)
  if(!isTRUE(info$isPCP)) {
    warning("cluster radius is only defined for cluster processes", call.=FALSE)
    return(NA)
  }
  rmax <- info$range(..., thresh = thresh)
  if(precision && is.function(info$ddist)){
    ddist <- function(r) info$ddist(r, ...)
    prec <- integrate(ddist, 0, rmax)
    attr(rmax, "prec") <- prec
  }
  return(rmax)
}

## The method clusterradius.kppm is in spatstat.core

