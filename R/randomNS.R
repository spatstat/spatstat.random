##
##   randomNS.R
##
##   Simulating from Neyman-Scott process
##   'Naive' algorithm
##
##   $Revision: 1.39 $  $Date: 2024/06/09 00:11:29 $
##
##    Original code for naive simulation of Neyman-Scott by Adrian Baddeley
##    Original code for rCauchy and rVarGamma offspring by Abdollah Jalilian
##    Other code and modifications by Adrian Baddeley
##    Bug fixes by Abdollah, Adrian, and Rolf Turner
##
##   Copyright (c) 2000-2023 Adrian Baddeley and Abdollah Jalilian
##   GNU Public Licence >= 2.0

rNeymanScott <- 
  function(kappa, expand, rcluster, win = unit.square(),
           ..., nsim=1, drop=TRUE, nonempty=TRUE, saveparents=TRUE,
           kappamax=NULL, mumax=NULL)
{
  ## Generic Neyman-Scott process
  ## Implementation for bounded cluster radius
  ##
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  if(nsim == 0) return(simulationresult(list()))

  ## Catch old argument name 'rmax'
  if(missing(expand) && !is.null(rmax <- list(...)$rmax)) {
    warning("outdated usage in rNeymanScott: 'rmax' should be 'expand'")
    expand <- rmax
  }
    
  ## Catch old argument name 'lmax'
  if(is.null(kappamax) && !is.null(lmax <- list(...)$lmax)) {
    warning("outdated usage in rNeymanScott: 'lmax' should be 'kappamax'")
    kappamax <- lmax
  }
    
  ## 'rcluster' may be
  ##
  ##     (1) a function(x,y, ...) that takes the coordinates
  ##         (x,y) of the parent point and generates a list(x,y) of offspring
  ##

  if(is.function(rcluster))
    return(rPoissonCluster(kappa, expand, rcluster, win, ...,
                           kappamax=kappamax, nsim=nsim, drop=drop,
                           saveparents=saveparents))

  ##     (2) a list(mu, f) where mu is a numeric value, function, or pixel image
  ##         and f is a function(n, ...) generating n i.i.d. offspring at 0,0
  
  if(!(is.list(rcluster) && length(rcluster) == 2))
    stop("rcluster should be either a function, or a list of two elements")
  win <- as.owin(win)
  mu <- rcluster[[1]]
  rdisplace <- rcluster[[2]]
  if(!(is.numeric(mu) || is.im(mu) || is.function(mu))) 
    stop("rcluster[[1]] should be a number, a function or a pixel image")  
  if(is.numeric(mu) && !(length(mu) == 1 && mu >= 0))
    stop("rcluster[[1]] should be a single nonnegative number")
  if(!is.function(rdisplace))
    stop("rcluster[[2]] should be a function")
  if(is.null(mumax))
    mumax <- if(is.numeric(mu)) mu else
             if(is.im(mu)) max(mu) else
             (1.05 * max(as.im(mu, W=win, ..., strict=TRUE)))

  ## Generate parents in dilated window
  frame <- boundingbox(win)
  dilated <- grow.rectangle(frame, expand)
  if(is.im(kappa) && !is.subset.owin(dilated, as.owin(kappa)))
    stop(paste("The window in which the image",
               sQuote("kappa"),
               "is defined\n",
               "is not large enough to contain the dilation of the window",
               sQuote("win")))
  if(nonempty) {
    if(is.function(kappa)) {
      kappa <- as.im(kappa, W=dilated, ..., strict=TRUE)
      kappamax <- NULL
    }
    ## intensity of parents with at least one offspring point
    kappa <- kappa * (1 - exp(-mumax))
  }
  ## generate
  parentlist <- rpoispp(kappa, lmax=kappamax, win=dilated, nsim=nsim,
                        drop=FALSE, warnwin=FALSE)

  resultlist <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {

    ## if(i > 1) gc(FALSE)

    parents <- parentlist[[i]]
    np <- npoints(parents)
    
    ## generate cluster sizes
    if(np == 0) {
      ## no parents - empty pattern
      result <- ppp(numeric(0), numeric(0), window=win)
      parentid <- integer(0)
      noff <- 0
    } else {
      if(!nonempty) {
        ## cluster sizes are Poisson
        csize <- rpois(np, mumax)
      } else {
        ## cluster sizes are Poisson conditional on > 0
        csize <- qpois(runif(np, min=dpois(0, mumax)), mumax)
      }
      noff <- sum(csize)
      xparent <- parents$x
      yparent <- parents$y
      x0 <- rep.int(xparent, csize)
      y0 <- rep.int(yparent, csize)
      ## invoke random generator
      dd <- rdisplace(noff, ...)
      mm <- if(is.ppp(dd)) marks(dd) else NULL
      ## validate
      xy <- xy.coords(dd)
      dx <- xy$x
      dy <- xy$y
      if(!(length(dx) == noff))
        stop("rcluster returned the wrong number of points")
      ## create offspring and offspring-to-parent map
      xoff <- x0 + dx
      yoff <- y0 + dy
      parentid <- rep.int(1:np, csize)
      ## trim to window
      retain <- inside.owin(xoff, yoff, win)
      if(is.im(mu))
        retain[retain] <- inside.owin(xoff[retain], yoff[retain], as.owin(mu))
      xoff <- xoff[retain]
      yoff <- yoff[retain]
      parentid <- parentid[retain]
      if(!is.null(mm)) mm <- marksubset(mm, retain)
      ## done
      result <- ppp(xoff, yoff, window=win, check=FALSE, marks=mm)
    }

    if(is.im(mu)) {
      ## inhomogeneously modulated clusters a la Waagepetersen
      result <- rthin(result, P=mu, Pmax=mumax)
    }

    if(saveparents) {
      attr(result, "parents") <- parents
      attr(result, "parentid") <- parentid
      attr(result, "expand") <- expand
      attr(result, "cost") <- np + noff
    }
    
    resultlist[[i]] <- result
  }

  result <- simulationresult(resultlist, nsim, drop)
  return(result)
}  

fakeNeyScot <- function(Y, lambda, win, saveLambda, saveparents) {
  ## Y is a ppp or ppplist obtained from rpoispp
  ## which will be returned as the realisation of a Neyman-Scott process
  ## when the process is degenerately close to Poisson.
  if(saveLambda || saveparents) {
    if(saveLambda && !is.im(lambda)) lambda <- as.im(lambda, W=win)
    if(saveparents) emptyparents <- ppp(window=win) # empty pattern
    if(isSingle <- is.ppp(Y)) Y <- solist(Y)
    for(i in seq_along(Y)) {
      Yi <- Y[[i]]
      if(saveLambda) attr(Yi, "Lambda") <- lambda
      if(saveparents) {
        attr(Yi, "parents") <- emptyparents
        attr(Yi, "parentid") <- integer(0)
        attr(Yi, "cost") <- npoints(Yi)
      }
      Y[[i]] <- Yi
    }
    if(isSingle) Y <- Y[[1L]]
  }
  return(Y)
}

thinParents <- function(X, P, Pmax=1) {
  ## Thin the parents and remove orphans
  Offspring <- X
  Parents <- attr(Offspring, "parents")
  Surname <- attr(Offspring, "parentid")
  retainparents <- rthin(Parents,
                         P=P, Pmax=Pmax,
                         what="fate",
                         na.zero=TRUE, fatal=FALSE)
  retainoffspring <- retainparents[Surname]
  Offspring <- Offspring[retainoffspring]
  attr(Offspring, "parents") <- list(x=Parents$x[retainparents],
                             y=Parents$y[retainparents])
  newserial <- cumsum(retainparents)
  attr(Offspring, "parentid") <- newserial[Surname[retainoffspring]]
  return(Offspring)
}


validate.kappa.mu <- function(kappa, mu, kappamax=NULL, mumax=NULL,
                              win, expand, ..., context="") {
  #' validate 'kappa' and 'mu' arguments
  if(!(is.numeric(mu) || is.im(mu) || is.function(mu)))
    stop(paste(context,
               "mu should be a number, a function or a pixel image"),
         call.=FALSE)
  if(is.numeric(mu)) {
    check.1.real(mu)
    check.finite(mu, xname="mu")
    stopifnot(mu >= 0)
  }
  if(!(is.numeric(kappa) || is.im(kappa) || is.function(kappa)))
    stop(paste(context,
               "kappa should be a number, a function or a pixel image"),
         call.=FALSE)
  if(is.numeric(kappa)) {
    check.1.real(kappa)
    check.finite(kappa, xname="kappa")
    stopifnot(kappa >= 0)
  }

  if(is.im(kappa)) {
    #' check domain 
    frame <- boundingbox(win)
    dilated <- grow.rectangle(frame, expand)
    if(!is.subset.owin(dilated, as.owin(kappa)))
      stop(paste("The window in which the image",
                 sQuote("kappa"),
                 "is defined\n",
                 "is not large enough to contain",
                 "the dilation of the window",
                 sQuote("win")),
             call.=FALSE)
  }

  ## get upper bounds on kappa and mu
  if(is.null(kappamax)) {
    if(is.numeric(kappa)) {
      kappamax <- kappa
    } else if(is.im(kappa)) {
      kappamax <- max(kappa)
    } else if(is.function(kappa)) {
      ## rough upper bound
      frame <- boundingbox(win)
      dilated <- grow.rectangle(frame, expand)
      dotargs <- list(...)
      kargs <- dotargs[names(dotargs) %in% names(args(kappa))]
      kim <- do.call(as.im,
                     append(list(kappa, W=dilated, strict=TRUE),
                            kargs))
      kra <- range(kim)
      kappamax <- kra[2] + 0.05 * diff(kra)
    }
  }
  if(is.null(mumax)) {
    if(is.numeric(mu)) {
      mumax <- mu
    } else if(is.im(mu)) {
      mumax <- max(mu)
    } else if(is.function(mu)) {
      ## rough upper bound
      mim <- as.im(mu, W=win, ..., strict=TRUE)
      mra <- range(mim)
      mumax <- mra[2] + 0.05 * diff(mra)
    } 
  }
  return(list(kappamax=kappamax, mumax=mumax))
}

