#'
#'   rLGCP.R
#'
#'   simulation of log-Gaussian Cox process
#'
#'  original code by Abdollah Jalilian
#'
#'  modifications by Adrian Baddeley, Ege Rubak and Tilman Davies
#' 
#'  $Revision: 1.26 $    $Date: 2023/05/02 06:19:29 $
#'

rLGCP <- local({

  rLGCP <- function(model="exp", mu = 0, param = NULL, ...,
                    win=NULL, saveLambda=TRUE, nsim=1, drop=TRUE) {
    ## validate
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
    do.rLGCP(model=model, mu=mu, param=param, ...,
             win=win, saveLambda=saveLambda, nsim=nsim, drop=drop)
  }

  do.rLGCP <- function(model="exp", mu = 0, param = NULL, ...,
                       win=NULL, saveLambda=TRUE,
                       eps = NULL, dimyx = NULL, xy = NULL,
                       rule.eps = c("adjust.eps", "grow.frame", "shrink.frame"),
                       modelonly=FALSE, Lambdaonly=FALSE,
                       nsim=1, drop=TRUE) {
    ## make RF model object from RandomFields package
    ## get the 'model generator'
    modgen <- getRandomFieldsModelGen(model)
    ## now create a RandomFields 'model' object
    rfmodel <- do.call(modgen, append(as.list(param), list(...)))
    if(!inherits(rfmodel, "RMmodel"))
      stop("Unable to create RandomFields model object", call.=FALSE)

    ## undocumented exit - return the RandomFields model object only
    if(modelonly)
      return(rfmodel)

    ## empty list
    if(nsim == 0) return(simulationresult(list()))
    
    ## simulation window
    win.given <- !is.null(win)
    mu.image <- is.im(mu)
    win <- if(win.given) as.owin(win) else if(mu.image) as.owin(mu) else owin()
  
    if(win.given && mu.image && !is.subset.owin(win, as.owin(mu)))
      stop(paste("The spatial domain of the pixel image", sQuote("mu"),
                 "does not cover the simulation window", sQuote("win")))

    ## convert win to a mask
    rule.eps <- match.arg(rule.eps)
    w <- as.mask(w=win, eps=eps, dimyx=dimyx, xy=xy, rule.eps=rule.eps)
    xcol <- w$xcol
    yrow <- w$yrow
    dimw <- w$dim

    ## evaluate 'mu' at pixels of mask
    if(is.numeric(mu)) {
      muxy <- mu
    } else {
      xy <- rasterxy.mask(w, drop=FALSE)
      xx <- xy$x
      yy <- xy$y
      muxy <- if (is.function(mu)) mu(xx,yy) else
              lookup.im(mu, xx, yy, naok=TRUE, strict=TRUE)
      muxy[is.na(muxy)] <- -Inf
    }
    ## corresponding image template
    Lambda <- as.im(w)

    ## generate 'nsim' realisations of a zero-mean Gaussian random field Z
    spc <- RandomFields::RFoptions()$general$spConform
    if(spc) RandomFields::RFoptions(spConform=FALSE)
    z <- RandomFields::RFsimulate(rfmodel, xcol, yrow, grid = TRUE, n=nsim)
    if(spc) RandomFields::RFoptions(spConform=TRUE)
    if(is.null(dim(z)))
      stop("RFsimulate did not return a matrix or array", call.=FALSE)
    ## ensure 3D array
    if(length(dim(z)) == 2) z <- array(z, dim=c(dim(z), 1))
    ## transform to spatstat convention 
    z <- aperm(z, c(2,1,3))
    ## safety checks
    if(!all(dim(z)[1:2] == dim(Lambda))) 
      stop("Internal error: wrong matrix dimensions in rLGCP", call.=FALSE)

    if(Lambdaonly) {
      ## undocumented exit - return Lambda only
      Lambdalist <- vector(mode="list", length=nsim)
      for(i in seq_len(nsim)) {
        ## Extract i-th realisation of Z; convert to log-Gaussian image
        Lambda$v[] <- exp(muxy + z[,,i])
        ## save as i-th realisation of Lambda
        Lambdalist[[i]] <- Lambda
      }
      return(simulationresult(Lambdalist, nsim, drop))
    }
    
    ## generate realisations of LGCP
    result <- vector(mode="list", length=nsim)
    for(i in seq_len(nsim)) {
      ## Extract i-th realisation of Z; convert to log-Gaussian image
      Lambda$v[] <- exp(muxy + z[,,i])
      ## generate Poisson points
      X <- rpoispp(Lambda)[win]
      ## 
      if(saveLambda)
        attr(X, "Lambda") <- Lambda
      result[[i]] <- X
    }
    return(simulationresult(result, nsim, drop))
  }

  rLGCP
})
