#'        clusterinfo.R
#' 
#'   Lookup table of information about cluster processes and Cox processes
#'
#'   $Revision: 1.64 $ $Date: 2024/06/09 00:05:33 $
#'
#'   Information is extracted by calling
#'             spatstatClusterModelInfo(<name>)
#'   where <name> is the name of the cluster type (e.g. 'Thomas', 'MatClust', 'Cauchy').
#'         
#'   Information is stored in the named list .Spatstat.ClusterModelInfoTable
#'   with names 'Thomas', 'MatClust' etc.
#' 
#'   Each list entry contains information about a particular cluster mechanism.
#'   It is a list with the following entries:
#'
#'       modelname       (String) The name of the process as it would be stated in an article
#'                       e.g. "Neyman-Scott process with Cauchy kernel"
#' 
#'       descname        (String) abbreviated name of the process for use in fv objects
#'                       e.g. "Cauchy process"
#' 
#'       modelabbrev     (String) short name of the process for use in kppm object
#'                       e.g. "Cauchy process"
#' 
#'       printmodelname  function(obj)
#'                       Invoked by print.kppm to construct the name of the process
#'                       (where the name includes shape parameters of the kernel)
#'                       e.g. produces value "Variance-Gamma process with nu=0.3"
#'
#'       parnames        (Character vector of length 2)
#'                       The names of the "NATIVE" cluster parameters
#'                       e.g. Thomas -> c('kappa', 'sigma2')
#'                            MatClust -> c('kappa', 'R')
#'
#'                       **NOTE** that there are two parametrisations:
#'                             - the 'GENERIC' parameters (e.g. 'kappa' and 'scale' for cluster models)
#'                             - the 'NATIVE' parameters depend on the model.
#' 
#'                       The native parameters are designed to be a more efficient representation
#'                       for internal use when fitting models. Starting values for the parameters
#'                       are expected to be given in the native parametrisation.
#'
#'       shapenames      (Character vector)
#'                       The names of any shape parameters of the kernel
#'                       that could/should be provided by the user
#'                       e.g.   'nu'
#'
#'       clustargsnames  (Character vector)
#'                       DEPRECATED
#'                       Identical to shapenames
#'       
#'       checkpar        function(par, native=TRUE, ...)
#'                       Validates the parameters 'par' (in either format)
#'                       and converts them to the native parameters (if native=TRUE)
#'                       or converts them to the generic parameters (if native=FALSE).
#'                       Transitional syntax: function(par, native=old, ..., old=TRUE)
#'                                            ('old' is a synonym for 'native')
#'
#'       native2generic  function(par)
#'                       Streamlined function to convert 'par' from native to generic format.
#'
#'       outputshape     function(margs, ..)
#'                       Convert 'margs' to the format required for printed output.
#' 
#'       checkclustargs  function(margs, old=TRUE, ...)
#'                       DEPRECATED: equivalent to outputshape(...)
#'                       If old=TRUE, return 'margs' (NEVER USED)
#'                       If old=FALSE, convert 'margs' to the format required for output.
#'
#'       resolveshape    function(...) 
#'                       Extracts any shape parameters that may be present in '...'
#'                       under different aliases or formats (e.g. 'nu.pcf' versus 'nu.ker').
#'                       Returns list(margs, covmodel) where covmodel=list(type, model, margs)
#'                       where 'margs' is in the format required for internal use.
#' 
#'       resolvedots     function(...)
#'                       DEPRECATED - equivalent to resolveshape(...)
#'
#'       parhandler      function(...)
#'                       DEPRECATED - equivalent to function(...) { resolveshape(...)$covmodel }
#'
#'       ddist           function(r, scale, ...)
#'                       One-dimensional probability density of distance from parent to offspring
#'                       ('scale' is in generic format)
#'
#'       range           function(..., par, thresh)
#'                       Computes cluster radius
#'                       'par' is in generic format
#'                       (thresh = probability that parent-to-offspring distance exceeds this radius)
#'
#'       kernel          function(par, rvals, ..., margs)
#'                       [DEFINED ONLY FOR POISSON CLUSTER PROCESSES]
#'                       compute kernel (two-dimensional probability density of offspring)
#'                       as a function of distance from parent to offspring.
#'                       'par' is in native format
#'                       'margs' is a list of shape arguments, if required
#'
#'       isPCP           logical.
#'                       TRUE iff the model is a Poisson cluster process
#'
#'       iscompact       logical.
#'                       TRUE if the kernel has compact support.
#'
#'       roffspring      function(n, par, ..., margs)
#'                       [DEFINED ONLY FOR POISSON CLUSTER PROCESSES]
#'                       Random generator of cluster.
#'                       Generates n offspring of a parent at the origin.
#'
#'       K               function(par, rvals, ..., model, margs)
#'                       Compute K-function
#'                       'par' is in native format
#'                       Arguments 'model', 'margs' are required if there are shape parameters
#'                        
#'       pcf             function(par, rvals, ..., model, margs)
#'                       Compute pair correlation function
#'                       'par' is in native format
#'                       Arguments 'model', 'margs' are required if there are shape parameters
#'                       
#'       Dpcf            function(par, rvals, ..., model, margs)
#'                       Compute vector of partial derivatives of pair correlation function with respect to 'par'
#'                       'par' is in native format
#'                       Arguments 'model', 'margs' are required if there are shape parameters
#'
#'       funaux          DEPRECATED
#'                       List of additional functions used in computation
#'                       (These should now be defined as stand-alone objects)
#'
#'       selfstart       function(X)
#'                       Calculates reasonable default estimates of 'par' from point pattern X
#'                       Returns 'par' in native format
#'
#'       interpret       function(par, lambda)
#'                       Return a full set of model parameters in a meaningful format for printing
#'
#'       roffspring      function(n, par, ..., model, margs)
#'                       Generates random offspring of a parent at the origin
#'

#' ..................................................................................................
#'      This file defines each entry in the table separately, then creates the table
#' ..................................................................................................
#'  


## ................. general helper functions (exported) ....................

## The following function simplifies code maintenance
## (due to changes in subscripting behaviour in recent versions of R)

retrieve.param <- function(desired, aliases, ..., par=NULL) {
  ## Retrieve the generic parameter named <desired> (or one of its <aliases>)
  ## from (...) or from 'par'
  dots <- list(...)
  par  <- as.list(par) # may be empty
  dnames <- names(dots)
  pnames <- names(par)
  for(key in c(desired, aliases)) {
    if(key %in% dnames) return(dots[[key]])
    if(key %in% pnames) return(par[[key]])
  }
  ## failed
  nali <- length(aliases)
  if(nali == 0) {
    explain <- NULL
  } else {
    explain <- paren(paste("also tried", ngettext(nali, "alias", "aliases"), commasep(sQuote(aliases))))
  }
  mess <- paste("Unable to retrieve argument", sQuote(desired), explain)
  stop(mess, call.=FALSE)
}

detect.par.format <- function(par, native, generic) {
  a <- check.named.vector(par, native, onError="null", xtitle="par")
  if(!is.null(a)) return("native")
  a <- check.named.vector(par, generic, onError="null", xtitle="par")
  if(!is.null(a)) return("generic")
  whinge <- paste("'par' should be a named vector with elements",
                  paren(paste(sQuote(native), collapse=" and "), "["),
                  "or",
                  paren(paste(sQuote(generic), collapse=" and "), "["))
  stop(whinge, call.=FALSE)
}


#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#'  >>>>>>>>>>>>>>>>>>>>>> Thomas process <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

.ThomasInfo <- list(
  modelname      = "Thomas process", # In modelname field of mincon fv obj.
  descname       = "Thomas process", # In desc field of mincon fv obj.
  modelabbrev    = "Thomas process", # In fitted kppm obj.
  printmodelname = function(...) "Thomas process", # Used by print.kppm
  parnames       = c("kappa", "sigma2"),
  shapenames     = NULL,
  clustargsnames = NULL,
  checkpar = function(par, native = old, ..., old=TRUE, strict=TRUE){
    ## 'par' is in either format
    if(is.null(par))
      par <- c(kappa=1,scale=1)
    if(strict && any(par<=0))
      stop("par values must be positive.", call.=FALSE)
    fmt <- detect.par.format(par, native=c("kappa", "sigma2"), generic=c("kappa", "scale"))
    if(fmt == "generic" && native) {
      ## convert generic to native
      par[2L] <- par[2L]^2
      names(par)[2L] <- "sigma2"
    } else if(fmt == "native" && !native) {
      ## convert native to generic
      par[2L] <- sqrt(par[2L])
      names(par)[2L] <- "scale"
    }
    return(par)
  },
  native2generic = function(par) {
    par[2L] <- sqrt(par[2L])
    names(par) <- c("kappa", "scale")
    return(par)
  },
  outputshape = function(margs, ...) list(),
  checkclustargs = function(margs, old = TRUE, ...) list(),
  resolveshape = function(...){ return(list(...)) },
  resolvedots = function(...){ return(list(...)) },
  parhandler = NULL,
  ## density function for the distance to offspring
  ddist = function(r, scale, ...) {
    ## 'scale' is generic format
    2 * pi * r * dnorm(r, 0, scale)/sqrt(2*pi*scale^2)
  },
  ## Practical range of clusters
  range = function(..., par=NULL, thresh=NULL){
    ## 'par' is in generic format
    scale <- retrieve.param("scale", "sigma", ..., par=par)
    if(!is.null(thresh)){
      ## The squared length of isotropic Gaussian (sigma)
      ## is exponential with mean 2 sigma^2
      rmax <- scale * sqrt(2 * qexp(thresh, lower.tail=FALSE))
    } else {
      rmax <- 4*scale
    }
    return(rmax)
  },
  kernel = function(par, rvals, ...) {
    ## 'par' is in native format
    scale <- sqrt(par[2L])
    dnorm(rvals, 0, scale)/sqrt(2*pi*scale^2)
  },
  isPCP=TRUE,
  iscompact=FALSE,
  roffspring=function(n, par, ...) {
    ## 'par' is in native format
    sigma <- sqrt(par[2L])
    list(x=rnorm(n, sd=sigma), y=rnorm(n, sd=sigma))
  },
  ## K-function
  K = function(par,rvals, ..., strict=TRUE){
    ## 'par' is in native format
    if(strict && any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    pi*rvals^2+(1-exp(-rvals^2/(4*par[2L])))/par[1L]
  },
  ## pair correlation function
  pcf= function(par,rvals, ..., strict=TRUE){
    ## 'par' is in native format
    if(strict && any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    1 + exp(-rvals^2/(4 * par[2L]))/(4 * pi * par[1L] * par[2L])
  },
  ## gradient of pcf (contributed by Chiara Fend)
  Dpcf= function(par,rvals, ..., strict=TRUE){
    ## 'par' is in native format
    if(strict && any(par <= 0)){
      dsigma2 <- rep.int(Inf, length(rvals))
      dkappa <- rep.int(Inf, length(rvals))
    } else {
      dsigma2 <- exp(-rvals^2/(4 * par[2L])) * (rvals/(4^2 * pi * par[1L] * par[2L]^3) - 1/(4 * pi * par[1L] * par[2L]^2))
      dkappa <- -exp(-rvals^2/(4 * par[2L]))/(4 * pi * par[1L]^2 * par[2L])
    }
    out <- rbind(dkappa, dsigma2)
    rownames(out) <- c("kappa","sigma2")
    return(out)
  },
  ## Convert to/from canonical cluster parameters
  tocanonical = function(par, ...) {
    ## 'par' is in native format
    ## convert to experimental 'canonical' format
    kappa <- par[[1L]]
    sigma2 <- par[[2L]]
    c(strength=1/(4 * pi * kappa * sigma2), scale=sqrt(sigma2))
  },
  tohuman = function(can, ...) {
    ## 'can' is in 'canonical' format
    ## convert to native format
    strength <- can[[1L]]
    scale <- can[[2L]]
    sigma2 <- scale^2
    c(kappa=1/(4 * pi * strength * sigma2), sigma2=sigma2)
  },
  ## sensible starting parameters
  selfstart = function(X) {
    ## return 'par' in native format
    kappa <- intensity(X)
    sigma2 <- 4 * mean(nndist(X))^2
    c(kappa=kappa, sigma2=sigma2)
  },
  ## meaningful model parameters
  interpret = function(par, lambda) {
    kappa <- par[["kappa"]]
    sigma <- sqrt(par[["sigma2"]])
    mu <- if(is.numeric(lambda) && length(lambda) == 1)
            lambda/kappa else NA
    c(kappa=kappa, sigma=sigma, mu=mu)
  },
  roffspring = function(n, par, ...) {
    sd <- sqrt(par[["sigma2"]])
    list(x=rnorm(n, sd=sd), y=rnorm(n, sd=sd))
  }
)

#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#'  >>>>>>>>>>>>>>>>> Matern cluster process <<<<<<<<<<<<<<<<<<<<<<<<<<
#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#' auxiliary functions
.MatClustHfun <- function(zz) {
  ok <- (zz < 1)
  h <- numeric(length(zz))
  h[!ok] <- 1
  z <- zz[ok]
  h[ok] <- 2 + (1/pi) * (
    (8 * z^2 - 4) * acos(z)
    - 2 * asin(z)
    + 4 * z * sqrt((1 - z^2)^3)
    - 6 * z * sqrt(1 - z^2)
  )
  return(h)
}

.MatClustDOH <- function(zz) {
  ok <- (zz < 1)
  h <- numeric(length(zz))
  h[!ok] <- 0
  z <- zz[ok]
  h[ok] <- (16/pi) * (z * acos(z) - (z^2) * sqrt(1 - z^2))
  return(h)
}

## g(z) = DOH(z)/z has a limit at z=0.
.MatClustgfun <- function(zz) {
  ok <- (zz < 1)
  h <- numeric(length(zz))
  h[!ok] <- 0
  z <- zz[ok]
  h[ok] <- (2/pi) * (acos(z) - z * sqrt(1 - z^2))
  return(h)
}

.MatClustgprime <- function(zz) {
  ok <- (zz < 1)
  h <- numeric(length(zz))
  h[!ok] <- 0
  z <- zz[ok]
  h[ok] <- -(2/pi) * 2 * sqrt(1 - z^2)
  return(h)
}

#' main list of info

.MatClustInfo = list(
  modelname      = "Matern cluster process", # In modelname field of mincon fv obj.
  descname       = "Matern cluster process", # In desc field of mincon fv obj.
  modelabbrev    = "Matern cluster process", # In fitted obj.
  printmodelname = function(...) "Matern cluster process", # Used by print.kppm
  parnames       = c("kappa", "R"),
  shapenames     = NULL,
  clustargsnames = NULL,
  checkpar = function(par, native = old, ..., old=TRUE, strict=TRUE){
    ## 'par' is in either format
    if(is.null(par))
      par <- c(kappa=1,scale=1)
    if(strict && any(par<=0))
      stop("par values must be positive.", call.=FALSE)
    detect.par.format(par, native=c("kappa", "R"), generic=c("kappa", "scale"))
    names(par)[2L] <- if(native) "R" else "scale"
    return(par)
  },
  native2generic = function(par) {
    names(par) <- c("kappa", "scale")
    return(par)
  },
  ## density function for the distance to offspring
  ddist = function(r, scale, ...) {
    ## 'scale' is generic format
    ifelse(r>scale, 0, 2 * r / scale^2)
  },
  ## Practical range of clusters
  range = function(..., par=NULL, thresh=NULL){
    ## 'par' is in generic format
    scale <- retrieve.param("scale", "R", ..., par=par)
    return(scale)
  },
  outputshape = function(margs,  ...) list(),
  checkclustargs = function(margs, old = TRUE, ...) list(),
  resolveshape = function(...){ return(list(...)) },
  resolvedots = function(...){ return(list(...)) },
  parhandler = NULL,
  kernel = function(par, rvals, ...) {
    ## 'par' is in native format
    scale <- par[2L]
    ifelse(rvals>scale, 0, 1/(pi*scale^2))
  },
  isPCP=TRUE,
  iscompact=TRUE,
  roffspring = function(n, par, ...) {
    ## 'par' is in native format
    R <- par[2L]
    rad <- R * sqrt(runif(n))
    theta <- runif(n, max=2*pi)
    list(x = rad * cos(theta), y = rad * sin(theta))
  },
  K = function(par,rvals, ...){
    ## 'par' is in native format
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    kappa <- par[1L]
    R <- par[2L]
    y <- pi * rvals^2 + (1/kappa) * .MatClustHfun(rvals/(2 * R))
    return(y)
  },
  pcf= function(par,rvals, ...){
    ## 'par' is in native format
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    kappa <- par[1L]
    R <- par[2L]
    y <- 1 + (1/(pi * kappa * R^2)) * .MatClustgfun(rvals/(2 * R))
    return(y)
  },
  Dpcf= function(par,rvals, ...){
    ## 'par' is in native format
    kappa <- par[1L]
    R <- par[2L]
    if(any(par <= 0)){
      dkappa <- rep.int(Inf, length(rvals))
      dR <- rep.int(Inf, length(rvals))
    } else {
      dkappa <- -.MatClustgfun(rvals/(2 * R)) / (pi * kappa^2 * R^2)
      dR <- -2*.MatClustgfun(rvals/(2 * R))/(pi * kappa * R^3) - (1/(pi * kappa * R^2)) * .MatClustgprime(rvals/(2 * R))*rvals/(2*R^2)
    }
    out <- rbind(dkappa, dR)
    rownames(out) <- c("kappa","R")
    return(out)
  },         
  ## Convert to/from canonical cluster parameters
  tocanonical = function(par, ...) {
    ## 'par' is in native format
    ## convert to experimental 'canonical' format
    kappa <- par[[1L]]
    R <- par[[2L]]
    c(strength=1/(pi * kappa * R^2), scale=R)
  },
  tohuman = function(can, ...) {
    ## 'can' is in 'canonical' format
    ## convert to native format
    strength <- can[[1L]]
    scale <- can[[2L]]
    c(kappa=1/(pi * strength * scale^2), R=scale)
  },
  ## sensible starting paramters
  selfstart = function(X) {
    ## return 'par' in native format
    kappa <- intensity(X)
    R <- 2 * mean(nndist(X)) 
    c(kappa=kappa, R=R)
  },
  ## meaningful model parameters
  interpret = function(par, lambda) {
    kappa <- par[["kappa"]]
    R     <- par[["R"]]
    mu    <- if(is.numeric(lambda) && length(lambda) == 1)
               lambda/kappa else NA           
    c(kappa=kappa, R=R, mu=mu)
  },
  roffspring = function(n, par, ...) {
    R <- par[["R"]]
    rad <- R * sqrt(runif(n))
    theta <- runif(n, max=2*pi)
    list(x = rad * cos(theta), y = rad * sin(theta))
  }
)

#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#'  >>>>>>>>>>>>>>>>>>>>> Cauchy kernel cluster process <<<<<<<<<<<<<<<
#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


.CauchyInfo <- list(
  modelname      = "Neyman-Scott process with Cauchy kernel", # In modelname field of mincon fv obj.
  descname       = "Neyman-Scott process with Cauchy kernel", # In desc field of mincon fv obj.
  modelabbrev    = "Cauchy process", # In fitted obj.
  printmodelname = function(...) "Cauchy process", # Used by print.kppm
  parnames       = c("kappa", "eta2"),
  shapenames     = NULL,
  clustargsnames = NULL,
  checkpar = function(par, native = old, ..., old=TRUE, strict=TRUE){
    ## 'par' is in either format
    if(is.null(par))
      par <- c(kappa=1,scale=1)
    if(strict && any(par<=0))
      stop("par values must be positive.", call.=FALSE)
    fmt <- detect.par.format(par, native=c("kappa", "eta2"), generic=c("kappa", "scale"))
    if(fmt == "generic" && native) {
      ## convert generic to native
      ## eta2 = 4 * scale^2
      par[2L] <- (2*par[2L])^2
      names(par)[2L] <- "eta2"
    } else if(fmt == "native" && !native) {
      ## convert native to generic
      ## scale = sqrt(eta2/4)
      par[2L] <- sqrt(par[2L])/2
      names(par)[2L] <- "scale"
    }
    return(par)
  },
  native2generic = function(par) {
    par[2L] <- sqrt(par[2L])/2
    names(par) <- c("kappa", "scale")
    return(par)
  },
  outputshape = function(margs, ...) list(),
  checkclustargs = function(margs, old = TRUE, ...) list(),
  resolveshape = function(...){ return(list(...)) },
  resolvedots = function(...){ return(list(...)) },
  parhandler = NULL,
  ## density function for the distance to offspring
  ddist = function(r, scale, ...) {
    ## 'scale' is generic format
    r/(scale^2) *  (1 + (r / scale)^2)^(-3/2)
  },
  ## Practical range of clusters
  range = function(..., par=NULL, thresh=0.01){
    ## 'par' is in generic format
    thresh <- as.numeric(thresh %orifnull% 0.01)
    scale <- retrieve.param("scale", character(0), ..., par=par)
    ## integral of ddist(r) dr is 1 - (1+(r/scale)^2)^(-1/2)
    ## solve for integral = 1-thresh:
    rmax <- scale * sqrt(1/thresh^2 - 1)
    return(rmax)
  },
  kernel = function(par, rvals, ...) {
    ## 'par' is in native format
    scale <- sqrt(par[2L])/2
    1/(2*pi*scale^2)*((1 + (rvals/scale)^2)^(-3/2))
  },
  isPCP=TRUE,
  iscompact=FALSE,
  roffspring=function(n, par, ...) {
    ## 'par' is in native format
    rate <- par[["eta2"]]/8 
    b <- 1/sqrt(rgamma(n, shape=1/2, rate=rate))
    list(x = b * rnorm(n), y = b * rnorm(n))
  },
  K = function(par,rvals, ...){
    ## 'par' is in native format
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    pi*rvals^2 + (1 - 1/sqrt(1 + rvals^2/par[2L]))/par[1L]
  },
  pcf= function(par,rvals, ...){
    ## 'par' is in native format
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    1 + ((1 + rvals^2/par[2L])^(-1.5))/(2 * pi * par[2L] * par[1L])
  },
  Dpcf= function(par,rvals, ...){
    ## 'par' is in native format
    if(any(par <= 0)){
      dkappa <- rep.int(Inf, length(rvals))
      deta2 <- rep.int(Inf, length(rvals))
    } else {
      dkappa <- -(1 + rvals^2/par[2L])^(-1.5)/(2 * pi * par[2L] * par[1L]^2)
      deta2 <- 1.5 * rvals^2 * (1 + rvals^2/par[2L])^(-2.5)/(2 * par[2L]^3 * par[1L] * pi) - (1 + rvals^2/par[2L])^(-1.5)/(2*pi*par[1L]*par[2L]^2)
    }
    out <- rbind(dkappa, deta2)
    rownames(out) <- c("kappa","eta2")
    return(out)
  },
  ## Convert to/from canonical cluster parameters
  tocanonical = function(par, ...) {
    ## 'par' is in native format
    ## convert to experimental 'canonical' format
    kappa <- par[[1L]]
    eta2 <- par[[2L]]
    c(strength=1/(2 * pi * kappa * eta2), scale=sqrt(eta2)/2)
  },
  tohuman = function(can, ...) {
    ## 'can' is in 'canonical' format
    ## convert to native format
    strength <- can[[1L]]
    scale <- can[[2L]]
    eta2 <- 4 * scale^2
    c(kappa=1/(2 * pi * strength * eta2), eta2=eta2)
  },
  selfstart = function(X) {
    ## return 'par' in native format
    kappa <- intensity(X)
    eta2 <- 4 * mean(nndist(X))^2
    c(kappa = kappa, eta2 = eta2)
  },
  ## meaningful model parameters
  interpret = function(par, lambda) {
    #' par is in native format
    kappa <- par[["kappa"]]
    omega <- sqrt(par[["eta2"]])/2
    mu <- if(is.numeric(lambda) && length(lambda) == 1)
            lambda/kappa else NA
    c(kappa=kappa, omega=omega, mu=mu)
  },
  roffspring = function(n, par, ...) {
    #' par is in native format
    rate <- par[["eta2"]]/8 
    b <- 1/sqrt(rgamma(n, shape=1/2, rate=rate))
    list(x = b * rnorm(n), y = b * rnorm(n))
  }
)


#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#'  >>>>>>>>>>>>>>>>> Variance Gamma kernel cluster process <<<<<<<<<<<
#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#' helper functions

resolve.vargamma.shape <- function(...,
                                   nu.ker=NULL,  nu.pcf=NULL,
                                   nu = NULL, allow.nu = FALSE,
                                   allow.default = FALSE) {
  ## ingest any kind of 'nu' argument by name
  if(is.null(nu.ker) && is.null(nu.pcf)) {
    if(allow.nu && !is.null(nu)) {
      nu.ker <- nu
    } else if(allow.default) {
      nu.ker <- -1/4
    } else stop("Must specify nu.ker or nu.pcf", call.=FALSE)
  }
  if(!is.null(nu.ker) && !is.null(nu.pcf))
    stop("Only one of nu.ker and nu.pcf should be specified",
         call.=FALSE)
  if(!is.null(nu.ker)) {
    check.1.real(nu.ker)
    stopifnot(nu.ker > -1/2)
    nu.pcf <- 2 * nu.ker + 1
  } else {
    check.1.real(nu.pcf)
    stopifnot(nu.pcf > 0)
    nu.ker <- (nu.pcf - 1)/2
  }
  return(list(nu.ker=nu.ker, nu.pcf=nu.pcf))
}

.VarGammaResolveShape <- function(...){
  nu.ker <- resolve.vargamma.shape(...,
                                   allow.default=TRUE, allow.nu=TRUE)$nu.ker
  check.1.real(nu.ker)
  stopifnot(nu.ker > -1/2)
  margs <- list(nu.ker=nu.ker, nu.pcf=2*nu.ker+1)
  cmodel <- list(type="Kernel", model="VarGamma", margs=margs)
  out <- list(margs = margs, covmodel = cmodel)
  return(out)
}

.VarGammaOutputShape <- function(margs, ...) { list(nu=margs$nu.ker) }

.VarGammaDdist <- function(r, scale, nu, ...) {
  ## one-dimensional probability density function for the distance to offspring
  ## 'scale' is generic format
  numer <- ((r/scale)^(nu+1)) * besselK(r/scale, nu)
  numer[r==0] <- 0
  denom <- (2^nu) * scale * gamma(nu + 1)
  y <- numer/denom
  y[!is.finite(y)] <- 0
  return(y)
}

.VarGammaPdist <- function(r, scale, nu, ...) {
  n <- length(r)
  z <- numeric(n)
  for(i in seq_len(n)) 
    z[i] <- integrate(.VarGammaDdist, 0, r[i], scale=scale, nu=nu)$value
  return(z)
}

.VarGammaPdistContrast <- function(r, scale, nu, p, ...) {
  .VarGammaPdist(r, scale, nu, ...) - p
}

.VarGammaKintegrand <- function(x, par, nu.pcf) {
  ## x * pcf(x) without check on par values
  numer <- (x/par[2L])^nu.pcf * besselK(x/par[2L], nu.pcf)
  denom <- 2^(nu.pcf+1) * pi * par[2L]^2 * par[1L] * gamma(nu.pcf + 1)
  return(x * (1 + numer/denom))
}

#'    main list of info
#'   NOTE: 'nu' represents 'nu.ker',
#'         the parameter \nu^\prime in equation (12)
#'         of Jalilian, Guan & Waagepetersen 2013

.VarGammaInfo <- list(
  modelname      = "Neyman-Scott process with Variance Gamma kernel", # In modelname field of mincon fv obj.
  descname       = "Neyman-Scott process with Variance Gamma kernel", # In desc field of mincon fv obj.
  modelabbrev    = "Variance Gamma process", # In fitted obj.
  printmodelname = function(obj){ # Used by print.kppm
    paste0("Variance Gamma process (nu=",
           signif(obj$clustargs[["nu"]], 2), ")")
  },
  parnames = c("kappa", "eta"),
  shapenames     = "nu",
  clustargsnames = "nu",
  checkpar = function(par, native = old, ..., old = TRUE, strict=TRUE){
    ## 'par' is in either format
    if(is.null(par))
      par <- c(kappa=1,scale=1)
    if(strict && any(par<=0))
      stop("par values must be positive.", call.=FALSE)
    detect.par.format(par, native=c("kappa", "eta"), generic=c("kappa", "scale"))
    names(par)[2L] <- if(native) "eta" else "scale"
    return(par)
  },
  native2generic = function(par) {
    names(par) <- c("kappa", "scale")
    return(par)
  },
  outputshape = .VarGammaOutputShape,
  checkclustargs = function(margs, old = TRUE, ...){
    if(!old)
      margs <- list(nu=margs$nu.ker)
    return(margs)
  },
  resolveshape = .VarGammaResolveShape,
  resolvedots  = .VarGammaResolveShape,
  parhandler = function(...){ .VarGammaResolveShape(...)$covmodel },
  ddist = .VarGammaDdist,
  ## Practical range of clusters
  range = function(..., par=NULL, thresh=0.001){
    ## 'par' is in generic format
    thresh <- as.numeric(thresh %orifnull% 0.001)
    scale <- retrieve.param("scale", character(0), ..., par=par)
    ## Find value of nu:
    margs <- .VarGammaResolveShape(...)$margs
    nu <-    .VarGammaOutputShape(margs)$nu
    if(is.null(nu))
      stop(paste("Argument ", sQuote("nu"), " must be given."),
           call.=FALSE)
    rmax <- uniroot(.VarGammaPdistContrast, lower = scale, upper = 1000 * scale,
                    scale=scale, nu=nu, p=1-thresh)$root
    return(rmax)
  },
  ## kernel function in polar coordinates (no angular argument).
  kernel = function(par, rvals, ..., margs) {
    ## 'par' is in native format
    scale <- as.numeric(par[2L])
    nu.ker <- margs$nu.ker %orifnull% margs$nu
    ## evaluate
    numer <- ((rvals/scale)^nu.ker) * besselK(rvals/scale, nu.ker)
    numer[rvals==0] <- if(nu.ker > 0) 2^(nu.ker-1)*gamma(nu.ker) else Inf
    denom <- pi * (2^(nu.ker+1)) * scale^2 * gamma(nu.ker + 1)
    numer/denom
  },
  isPCP=TRUE,
  iscompact=FALSE,
  roffspring=function(n, par, ..., margs) {
    ## 'par' is in native format
    scale <- par[["eta"]] ## eta = omega = scale
    nu.ker <- margs$nu.ker
    alpha <- 2 * (nu.ker + 1)
    beta  <- 1/(2 * scale^2)
    sdee <- sqrt(rgamma(n, shape=alpha/2, rate=beta))
    list(x = sdee * rnorm(n),
         y = sdee * rnorm(n))
  },
  K = function(par,rvals, ..., margs){
    ## 'par' is in native format
    ## margs = list(.. nu.pcf.. )
    ## K function requires integration of pair correlation
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    nu.pcf <- margs$nu.pcf
    out <- numeric(length(rvals))
    for (i in which(rvals > 0))
      out[i] <- 2 * pi * integrate(.VarGammaKintegrand,
                                   lower=0, upper=rvals[i],
                                   par=par, nu.pcf=nu.pcf)$value
    return(out)
  },
  pcf= function(par,rvals, ..., margs){
    ## 'par' is in native format
    ## margs = list(..nu.pcf..)
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    nu.pcf <- margs$nu.pcf
    sig2 <- 1 / (4 * pi * (par[2L]^2) * nu.pcf * par[1L])
    denom <- 2^(nu.pcf - 1) * gamma(nu.pcf)
    rr <- rvals / par[2L]
    ## Matern correlation function
    fr <- ifelseXB(rr > 0,
                  (rr^nu.pcf) * besselK(rr, nu.pcf) / denom,
                  1)
    return(as.numeric(1 + sig2 * fr))
  },
  Dpcf = NULL,
  ## Convert to/from canonical cluster parameters
  tocanonical = function(par, ..., margs) {
    ## 'par' is in native format
    ## convert to experimental 'canonical' format
    kappa <- par[[1L]]
    eta <- par[[2L]]
    nu.pcf <- margs$nu.pcf
    c(strength=1/(4 * pi * nu.pcf * kappa * eta^2), scale=eta)
  },
  tohuman = function(can, ..., margs) {
    ## 'can' is in 'canonical' format
    ## convert to native format
    strength <- can[[1L]]
    eta <- scale <- can[[2L]]
    nu.pcf <- margs$nu.pcf
    c(kappa=1/(4 * pi * nu.pcf * strength * eta^2), eta=scale)
  },
  ## sensible starting values
  selfstart = function(X) {
    ## return 'par' in native format
    kappa <- intensity(X)
    eta <- 2 * mean(nndist(X))
    c(kappa=kappa, eta=eta)
  },
  ## meaningful model parameters
  interpret = function(par, lambda) {
    #' par is in native format
    kappa <- par[["kappa"]]
    omega <- par[["eta"]]
    mu <- if(is.numeric(lambda) && length(lambda) == 1)
            lambda/kappa else NA
    c(kappa=kappa, omega=omega, mu=mu)
  },
  roffspring = function(n, par, ..., margs) {
    #' par is in native format
    shape <- margs$nu.ker + 1
    scale <- par[[2L]]
    rate <- 1/(2 * scale^2)
    b <- sqrt(rgamma(n, shape=shape, rate=rate))
    list(x= b * rnorm(n), y = b * rnorm(n))
  }
)

#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#'  >>>>>>>>>>>>>>>>>>  Log-Gaussian Cox process <<<<<<<<<<<<<<<<<<<<<<
#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#' helper functions

.LGCPResolveShape <- function(...){
  ## resolve shape parameters for kppm and friends allowing for native/generic par syntax
  dots <- list(...)
  nam <- names(dots)
  out <- list()
  cmod <- dots$covmodel
  model <- cmod$model %orifnull% dots$model %orifnull% "exponential"
  margs <- NULL
  ## extract shape parameters and validate them
  switch(model,
         exponential = ,
         gauss = {
           ## no shape parameters
         },
         stable = ,
         fastStable = {
           stuff <- cmod %orifnull% dots
           ok <- "alpha" %in% names(stuff)
           if(!ok) stop("Parameter 'alpha' is required")
           margs <- stuff["alpha"]
           with(margs, {
             check.1.real(alpha)
             stopifnot(0 < alpha && alpha <= 2)
           })
         },
         gencauchy = ,
         fastGencauchy = {
           stuff <- cmod %orifnull% dots
           ok <- c("alpha", "beta") %in% names(stuff)
           if(!ok[1]) stop("Parameter 'alpha' is required")
           if(!ok[2]) stop("Parameter 'beta' is required")
           margs <- stuff[c("alpha", "beta")]
           with(margs, {
             check.1.real(alpha)
             check.1.real(beta)
             stopifnot(0 < alpha && alpha <= 2)
             stopifnot(beta > 0)
           })
         },
         matern = {
           stuff <- cmod %orifnull% dots
           ok <- "nu" %in% names(stuff)
           if(!ok) stop("Parameter 'nu' is required")
           margs <- stuff["nu"]
           with(margs, {
             check.1.real(nu)
           })
         },
         stop(paste("Covariance model", sQuote(model),
                    "is not yet supported"), call.=FALSE)
         )
  if(length(margs)==0) {
    margs <- NULL
  } else {
    ## detect anisotropic model
    if("Aniso" %in% names(margs))
      stop("Anisotropic covariance models cannot be used",
           call.=FALSE)
  }
  out$margs <- margs
  out$model <- model
  out$covmodel <- list(type="Covariance", model=model, margs=margs)
  return(out)
}

#'      main list of info

.LGCPInfo <- list(
  ## Log Gaussian Cox process: native par = (sigma2, alpha) 
  ## Log Gaussian Cox process: generic par = (var, scale) 
  modelname      = "Log-Gaussian Cox process", # In modelname field of mincon fv obj.
  descname       = "LGCP", # In desc field of mincon fv obj.
  modelabbrev    = "log-Gaussian Cox process", # In fitted obj.
  printmodelname = function(...) "log-Gaussian Cox process", # Used by print.kppm
  parnames       = c("sigma2", "alpha"),
  checkpar = function(par, native = old, ..., old=TRUE, strict=TRUE){
    ## 'par' is in either format
    if(is.null(par))
      par <- c(var=1,scale=1)
    if(strict && any(par<=0))
      stop("par values must be positive.", call.=FALSE)
    detect.par.format(par, native=c("sigma2", "alpha"), generic=c("var", "scale"))
    names(par) <- if(native) c("sigma2", "alpha") else c("var","scale")
    return(par)
  },
  native2generic = function(par) {
    names(par) <- c("var","scale")
    return(par)
  },
  outputshape = function(margs, ...) return(margs),
  checkclustargs = function(margs, old = TRUE, ...) return(margs),
  resolveshape = .LGCPResolveShape,
  resolvedots  = .LGCPResolveShape,
  parhandler   = function(...) { .LGCPResolveShape(...)$covmodel },
  isPCP=FALSE,
  iscompact=FALSE,
  roffspring=NULL,
  K = function(par, rvals, ..., model, margs) {
    ## 'par' is in native format
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    ## Determine covariance function 
    switch(model,
           exponential = {
             Cfun <- function(x) exp(-x)
           },
           gauss = ,
           fastGauss = {
             Cfun <- function(x) exp(-x^2)
           },
           stable = ,
           fastStable = {
             alpha <- margs[["alpha"]]
             Cfun <- function(x) exp(-x^alpha)
           },
           gencauchy = ,
           fastGencauchy = {
             alpha <- margs[["alpha"]]
             beta  <- margs[["beta"]]
             Cfun <- function(x) { (1 + x^alpha)^(-beta/alpha) }
           },
           matern = {
             nu <- margs[["nu"]]
             Cfun <- function(x) {
               z <- x * sqrt(2 * nu)
               ifelse(x == 0,
                      1, 
                      (z^nu) * besselK(z, nu) * (2^(1-nu))/gamma(nu))
             }
           },
           stop(paste("Model", sQuote(model), "is not recognised"))
           )
    ## hence determine integrand for K function
    integrand <- function(r) 2 * pi * r * exp(par[1L] * Cfun(r/par[2L]))
    ## compute indefinite integral
    imethod <- if(spatstat.options("fastK.lgcp")) "trapezoid" else "quadrature"
    th <- indefinteg(integrand, rvals, lower=0, method=imethod)
    return(th)
  },
  pcf= function(par, rvals, ..., model, margs) {
    ## 'par' is in native format
    if(any(par <= 0))
      return(rep.int(Inf, length(rvals)))
    ## Determine covariance function 
    switch(model,
           exponential = {
             Cfun <- function(x) exp(-x)
           },
           gauss = ,
           fastGauss = {
             Cfun <- function(x) exp(-x^2)
           },
           stable = ,
           fastStable = {
             alpha <- margs[["alpha"]]
             Cfun <- function(x) exp(-x^alpha)
           },
           gencauchy = ,
           fastGencauchy = {
             alpha <- margs[["alpha"]]
             beta  <- margs[["beta"]]
             Cfun <- function(x) { (1 + x^alpha)^(-beta/alpha) }
           },
           matern = {
             nu <- margs[["nu"]]
             Cfun <- function(x) {
               z <- x * sqrt(2 * nu)
               ifelse(x == 0,
                      1, 
                      (z^nu) * besselK(z, nu) * (2^(1-nu))/gamma(nu))
             }
           },
           stop(paste("Model", sQuote(model), "is not recognised"))
           )
    ## Hence evaluate pcf
    gtheo <- exp(par[1L] * Cfun(rvals/par[2L]))
    return(gtheo)
  },
  Dpcf= function(par,rvals, ..., model){
    ## 'par' is in native format
    if(!(model %in% c("exponential", "stable")))
      stop("Gradient of the pcf is not available for this model")
    if(model=="exponential") {
      dsigma2 <- exp(-rvals/par[2L]) * exp(par[1L]*exp(-rvals/par[2L]))
      dalpha <- rvals * par[1L] * exp(-rvals/par[2L]) * exp(par[1L]*exp(-rvals/par[2L]))/par[2L]^2
    } else if(model=="stable"){
      dsigma2 <- exp(-sqrt(rvals/par[2L])) * exp(par[1L]*exp(-sqrt(rvals/par[2L])))
      dalpha <- sqrt(rvals/par[2L]^3) * par[1L] * exp(-sqrt(rvals/par[2L])) * exp(par[1L] * exp(-sqrt(rvals/par[2L])))
    }
    out <- rbind(dsigma2, dalpha)
    rownames(out) <- c("sigma2","alpha")
    return(out)
  },
  ## sensible starting values
  selfstart = function(X) {
    ## return 'par' in native format
    alpha <- 2 * mean(nndist(X))
    c(sigma2=1, alpha=alpha)
  },
  ## meaningful model parameters
  interpret = function(par, lambda) {
    sigma2 <- par[["sigma2"]]
    alpha  <- par[["alpha"]]
    mu <- if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0)
            log(lambda) - sigma2/2 else NA
    c(sigma2=sigma2, alpha=alpha, mu=mu)
  }
)

#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#'              C O N S T R U C T     T A B L E
#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#'  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

.Spatstat.ClusterModelInfoTable <- 
  list(
    Thomas   = .ThomasInfo,
    MatClust = .MatClustInfo,
    Cauchy   = .CauchyInfo,
    VarGamma = .VarGammaInfo,
    LGCP     = .LGCPInfo)


spatstatClusterModelInfo <- function(name, onlyPCP = FALSE) {
  if(inherits(name, "detpointprocfamily")) {
    if(requireNamespace("spatstat.model")) {
      return(spatstat.model::spatstatDPPModelInfo(name))
    } else {
      message("The package 'spatstat.model' is required")
      return(NULL)
    }
  }
  if(!is.character(name) || length(name) != 1)
    stop("Argument must be a single character string", call.=FALSE)
  nama2 <- names(.Spatstat.ClusterModelInfoTable)
  if(onlyPCP){
    ok <- sapply(.Spatstat.ClusterModelInfoTable, getElement, name="isPCP")
    nama2 <- nama2[ok]
  } 
  if(!(name %in% nama2))
    stop(paste(sQuote(name), "is not recognised;",
               "valid names are", commasep(sQuote(nama2))),
         call.=FALSE)
  out <- .Spatstat.ClusterModelInfoTable[[name]]
  return(out)
}

