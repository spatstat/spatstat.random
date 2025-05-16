#'
#'   clusterprocess.R
#'
#'   $Revision: 1.2 $ $Date: 2025/05/16 07:14:33 $
#'


clusterprocess <- function(name="Thomas", ..., mu, kappa, scale) {
  if(missing(kappa)) stop("The parent intensity kappa must be given")
  if(missing(mu)) stop("The mean cluster size mu must be given")
  if(missing(scale)) stop("The cluster scale must be given")
  rules <- spatstatClusterModelInfo(name)
  if(!(rules$isPCP)) stop("clusterprocess only supports Neyman-Scott processes")
  par.std <- c(kappa=kappa, scale=scale)
  par.std <- rules$checkpar(par.std, native=FALSE)
  par.idio <- rules$checkpar(par.std, native=TRUE)
  other <- rules$resolveshape(...)
  clustargs <- rules$outputshape(other$margs)
  out <- list(name=name, rules=rules,
              par.std=par.std, par.idio=par.idio,
              mu=mu, clustargs=clustargs,
              other=other)
  class(out) <- "clusterprocess"
  return(out)
}

print.clusterprocess <- local({

  print.clusterprocess <- function(x, ...) {
    with(x, {
      splat(rules$printmodelname(list(clustargs=clustargs)))
      newpar <- rules$checkpar(par.std, native=FALSE)
      splat("Parent intensity kappa =", blurb("kappa", newpar["kappa"]))
      splat("Cluster scale = ", newpar["scale"])
      splat("Mean cluster size mu =", blurb("mu", mu))
      if(length(clustargs) > 0) {
        hdr <- paste("Cluster shape",
                     ngettext(length(clustargs), "parameter:", "parameters:"))
        if(is.list(clustargs) &&
           all(sapply(clustargs, is.numeric)) &&
           all(lengths(clustargs) == 1)) {
          splat(hdr,
                paste(names(clustargs), as.numeric(clustargs),
                      sep="=",
                      collapse=", "))
        } else {
          splat(hdr)
          print(clustargs)
        }
      }
    })
    return(invisible(NULL))
  }

  blurb <- function(name, value) {
    if(is.numeric(value)) as.character(value) else
    if(is.im(value)) "[image]" else "[unrecognized format]"
  }

  print.clusterprocess
})


                             
intensity.clusterprocess <- function(X, ...) {
  X$par.std[["kappa"]] * X$mu
}
  
predict.clusterprocess <- function(object, ...,
                                 locations,
                                 type="intensity",
                                 ngrid=NULL) {
  ## limited use!!!
  if(!identical(type, "intensity"))
    stop("Sorry, only type='intensity' is implemented")
  lambda <- object$par.std[["kappa"]] * object$mu
  if(is.numeric(lambda)) {
    ## stationary
    if(is.ppp(locations))
      return(rep(lambda, npoints(locations)))
    W <- as.owin(locations)
    if(!is.mask(W))
      W <- as.mask(W, dimyx=ngrid, ...)
    return(as.im(lambda, W=W))
  } else {
    ## nonstationary; lambda is an image
    if(missing(locations) || is.null(locations)) 
      return(lambda)
    return(lambda[locations, drop=FALSE])
  }
}

clusterradius.clusterprocess <- function(model, ...,
                                        thresh = NULL, precision = FALSE) {
  do.call(clusterradius.character,
          resolve.defaults(
            list(model = model$name, thresh = thresh, precision = precision),
            list(...),
            as.list(model$par.std), # sic
            model$clustargs)
          )
}

reach.clusterprocess <- function(x, ..., epsilon) {
  thresh <- if(missing(epsilon)) NULL else epsilon
  2 * clusterradius(x, ..., thresh=thresh)
}

simulate.clusterprocess <- function(object, nsim=1, ...,
                                    win=unit.square(), window=win) {
  with(object, {
    switch(name,
         Thomas = {
           rThomas(kappa=par.std[["kappa"]],
                   scale=par.std[["scale"]],
                   mu=mu,
                   win=window,
                   nsim=nsim,
                   ...)
         },
         MatClust = {
           rMatClust(kappa=par.std[["kappa"]],
                     scale=par.std[["scale"]],
                     mu=mu,
                     win=window,
                     nsim=nsim,
                     ...)
         },
         Cauchy = {
           rCauchy(kappa=par.std[["kappa"]],
                   scale=par.std[["scale"]],
                   mu=mu,
                   win=window,
                   nsim=nsim,
                   ...)
         },
         VarGamma = {
           do.call(rVarGamma,
                   resolve.defaults(
                     list(kappa=par.std[["kappa"]],
                          scale=par.std[["scale"]],
                          mu=mu,
                          win=window,
                          nsim=nsim,
                          ...),
                     clustargs))
         },
         stop(paste("Unrecognised cluster process model name",
                    sQuote(object$name)),
              call.=FALSE)
         )})
}

clusterkernel.clusterprocess <- function(model, ...) {
  info <- spatstatClusterModelInfo(model$name)
  k <- info$kernel
  par <- model$par.std
  clustargs <- model$clustargs
  f <- function(x, y=0, ...) { k(par, sqrt(x^2+y^2), margs=clustargs) }
  return(f)
}


#' for other methods, see also spatstat.model
