#'
#'   rLGCP.R
#'
#'   simulation of log-Gaussian Cox process
#'
#'  original code by Abdollah Jalilian
#'
#'  modifications by Adrian Baddeley, Ege Rubak and Tilman Davies
#' 
#'  $Revision: 1.30 $    $Date: 2023/10/20 15:36:57 $
#'

rLGCP <- local({

  rLGCP <- function(model=c("exponential", "gauss", "stable",
                            "gencauchy", "matern"),
                    mu = 0, param = NULL, ...,
                    win=NULL, saveLambda=TRUE, nsim=1, drop=TRUE) {
    ## validate
    model <- match.arg(model)
    if (is.numeric(mu)) {
      check.1.real(mu, paste("if", sQuote("mu"), "is numeric,"))
    } else if(!is.function(mu) && !is.im(mu))
      stop(paste(sQuote("mu"), "must be a constant, a function or an image"))
    check.1.integer(nsim)
    stopifnot(nsim >= 0)
    ## check for outdated usage
    if(!all(nzchar(names(param))))
      stop("Outdated syntax of argument 'param' to rLGCP", call.=FALSE)
    ##
    do.call(do.rLGCP,
            append(list(model=model,
                        mu=mu,
                        win=win,
                        saveLambda=saveLambda,
                        nsim=nsim,
                        drop=drop),
                   resolve.defaults(list(...), param)))
  }

  do.rLGCP <- function(model=c("exponential", "gauss", "stable",
                            "gencauchy", "matern"),
                       mu = 0, ...,
                       win=NULL, saveLambda=TRUE,
                       Lambdaonly=FALSE,
                       nsim=1, drop=TRUE) {

    ## empty list
    if(nsim == 0) return(simulationresult(list()))
    
    model <- match.arg(model)

    ## simulation window
    win.given <- !is.null(win)
    mu.image <- is.im(mu)
    win <- if(win.given) as.owin(win) else if(mu.image) as.owin(mu) else owin()
  
    if(win.given && mu.image && !is.subset.owin(win, as.owin(mu)))
      stop(paste("The spatial domain of the pixel image", sQuote("mu"),
                 "does not cover the simulation window", sQuote("win")))

    ## get shape parameters
    needed <- switch(model,
                     exponential = ,
                     gauss       = character(0),
                     stable      = "alpha",
                     gencauchy   = c("alpha", "beta"),
                     matern      = "nu")
    if(length(needed)) {
      stuff <- list(...)
      missed <- is.na(match(needed, names(stuff)))
      if(any(missed)) {
        nbad <- sum(missed)
        stop(paste(ngettext(nbad, "Parameter", "Parameters"),
                   commasep(sQuote(needed[missed])),
                   ngettext(nbad, "is", "are"), "required"),
             call.=FALSE)
      }
    }

    ## generate Gaussian Random Field
    Zlist <- switch(model,
                    exponential = {
                      rGRFexpo(W=win, mu=mu, 
                               ..., nsim=nsim, drop=FALSE)
                    },
                    gauss       = {
                      rGRFgauss(W=win, mu=mu,
                                ..., nsim=nsim, drop=FALSE)
                    },
                    stable = {
                      rGRFstable(W=win, mu=mu, 
                                 ..., nsim=nsim, drop=FALSE)
                    },
                    gencauchy = {
                      rGRFgencauchy(W=win, mu=mu, 
                                    ..., nsim=nsim, drop=FALSE)
                    },
                    matern = {
                      rGRFmatern(W=win, mu=mu, 
                                 ..., nsim=nsim, drop=FALSE)
                    },
                    stop(paste("Model", sQuote(model), "not matched")))

    if(length(Zlist) != nsim)
      stop("Internal error in generating realisations")
    
    ## exponentiate
    Lambdalist <- solapply(Zlist, exp)

    if(Lambdaonly) {
      ## undocumented exit - return Lambda only
      return(simulationresult(Lambdalist, nsim, drop))
    }
    
    ## generate realisations of LGCP
    result <- vector(mode="list", length=nsim)
    for(isim in seq_len(nsim)) {
      Lambda <- Lambdalist[[isim]]
      ## generate Poisson points
      X <- rpoispp(Lambda)[win]
      ## 
      if(saveLambda)
        attr(X, "Lambda") <- Lambda
      result[[isim]] <- X
    }
    return(simulationresult(result, nsim, drop))
  }

  rLGCP
})
