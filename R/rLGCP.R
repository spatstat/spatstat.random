#'
#'   rLGCP.R
#'
#'   simulation of log-Gaussian Cox process
#'
#'  original code by Abdollah Jalilian
#'
#'  modifications by Adrian Baddeley, Ege Rubak and Tilman Davies
#' 
#'  $Revision: 1.37 $    $Date: 2025/04/19 03:41:26 $
#'

rLGCP <- local({

  rLGCP <- function(model=c("exponential", "gauss", "stable",
                            "gencauchy", "matern"),
                    mu = 0, param = NULL, ...,
                    win=NULL, saveLambda=TRUE, nsim=1, drop=TRUE,
                    n.cond=NULL, w.cond=NULL) {
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
    do.call(dispatch.rLGCP,
            append(list(model=model,
                        mu=mu,
                        win=win,
                        saveLambda=saveLambda,
                        nsim=nsim,
                        drop=drop,
                        n.cond=n.cond),
                   resolve.defaults(list(...), param)))
  }

  dispatch.rLGCP <- function(..., n.cond=NULL, w.cond=NULL) {
    if(is.null(n.cond)) {
      do.rLGCP(...)
    } else {
      do.cond.rLGCP(..., n.cond=n.cond, w.cond=w.cond)
    }
  }
  
  do.rLGCP <- function(model=c("exponential", "gauss", "stable",
                            "gencauchy", "matern"),
                       mu = 0, ...,
                       win=NULL, saveLambda=TRUE,
                       LambdaOnly=FALSE, Lambdaonly=FALSE,
                       nsim=1, drop=TRUE,
                       n.cond=NULL) {

    ## empty list
    if(nsim == 0) return(simulationresult(list()))
    
    model <- match.arg(model)

    LambdaOnly <- LambdaOnly || Lambdaonly  # old spelling

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

    if(LambdaOnly) {
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

  do.cond.rLGCP <- function(model=c("exponential", "gauss", "stable",
                            "gencauchy", "matern"),
                            mu = 0, ...,
                            n.cond=NULL, w.cond=NULL, 
                            giveup=1000, maxchunk=100,
                            win=NULL, saveLambda=TRUE,
                            LambdaOnly=FALSE, Lambdaonly=FALSE,
                            nsim=1,
                            verbose=FALSE, drop=FALSE) {
    ## simulation window
    win.given <- !is.null(win)
    mu.image <- is.im(mu)
    win <- if(win.given) as.owin(win) else if(mu.image) as.owin(mu) else owin()

    LambdaOnly <- LambdaOnly || Lambdaonly  # old spelling
    
    ## type of simulation
    w.sim <- as.owin(win)
    fullwindow <- is.null(w.cond)
    if(fullwindow) {
      w.cond <- w.sim
      w.free <- NULL
    } else {
      #' not yet public
      stopifnot(is.owin(w.cond))
      w.free <- setminus.owin(w.sim, w.cond)
    }
  
    nremaining <- nsim
    ntried <- 0
    accept <- FALSE
    nchunk <- 1
    phistory <- Lamhistory <- numeric(0)
    results <- list()
    while(nremaining > 0) {
      ## increase chunk length
      nchunk <- min(maxchunk, giveup - ntried, 2 * nchunk)
      ## bite off next chunk of unconditional simulations
      lamlist <- do.rLGCP(model=model, mu=mu, win=win, ..., 
                          nsim=nchunk,
                          LambdaOnly=TRUE,
                          drop=FALSE)
      ## compute acceptance probabilities
      lamlist <- lapply(lamlist, "[", i=w.sim, drop=FALSE, tight=TRUE)
      if(fullwindow) {
        Lam <- sapply(lamlist, integral)
      } else {
        Lam <- sapply(lamlist, integral, domain=w.cond)
      }
      p <- exp(n.cond * log(Lam/n.cond) + n.cond - Lam)
      phistory <- c(phistory, p)
      Lamhistory <- c(Lamhistory, Lam)
      ## accept/reject
      accept <- (runif(length(p)) < p)
      if(any(accept)) {
        jaccept <- which(accept)
        if(length(jaccept) > nremaining)
          jaccept <- jaccept[seq_len(nremaining)]
        naccepted <- length(jaccept)
        if(verbose)
          splat("Accepted the",
                commasep(ordinal(ntried + jaccept)),
                ngettext(naccepted, "proposal", "proposals"))
        nremaining <- nremaining - naccepted
        for(j in jaccept) {
          lamj <- lamlist[[j]]
          if(min(lamj) < 0)
            lamj <- eval.im(pmax(lamj, 0))
          if(LambdaOnly) {
            #' undocumented: return the driving intensities only
            results <- append(results, list(lamj))
          } else {
            #' normal: return the simulated patterns
            if(fullwindow) {
              Y <- rpoint(n.cond, lamj, win=w.sim, forcewin=TRUE)
            } else {
              lamj.cond <- lamj[w.cond, drop=FALSE, tight=TRUE]
              lamj.free <- lamj[w.free, drop=FALSE, tight=TRUE]
              Ycond <- rpoint(n.cond, lamj.cond, win=w.cond)
              Yfree <- rpoispp(lamj.free)
              Y <- superimpose(Ycond, Yfree, W=w.sim)
            }
            if(saveLambda) attr(Y, "Lambda") <- lamj
            results <- append(results, list(Y))
          }
        }
      }
      ntried <- ntried + nchunk
      if(ntried >= giveup && nremaining > 0) {
        message(paste("Gave up after", ntried,
                      "proposals with", nsim - nremaining, "accepted"))
        message(paste("Mean acceptance probability =",
                      signif(mean(phistory), 3)))
        break
      }
    }
    nresults <- length(results)
    results <- simulationresult(results, nresults, drop)
    attr(results, "history") <- data.frame(Lam=Lamhistory, p=phistory)
    if(verbose && nresults == nsim)
      splat("Mean acceptance probability", signif(mean(phistory), 3))
    return(results)
  }
  
  rLGCP
})
