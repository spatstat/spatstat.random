#'
#'    rthomas.R
#'
#'   $Revision: 1.6 $ $Date: 2023/01/25 00:41:02 $
#' 
#'   Simulation of modified Thomas cluster process
#'   using either naive algorithm or BKBC algorithm
#'
#'   rThomasHom    Interface to C code for stationary case (BKBC)
#'   rThomas       General case (naive or BKBC)
#'
#'   Copyright (C) Adrian Baddeley, Rolf Turner and Ya-Mei Chang 2000-2023
#'   Licence: GNU Public Licence >= 2

rThomasHom <-function(kappa, mu, sigma, W=unit.square(), ..., nsim=1, drop=TRUE,
                      inflate=NULL, saveparents=FALSE, maxinflate=10) {
  check.1.real(kappa) && check.finite(kappa)
  check.1.real(mu) && check.finite(mu)
  check.1.real(sigma) && check.finite(sigma)
  check.1.integer(nsim)
  stopifnot(kappa >= 0)
  stopifnot(mu >= 0)
  stopifnot(sigma > 0)
  if(!is.null(inflate)) {
    check.1.real(inflate) && check.finite(inflate)
    stopifnot(inflate >= 1)
  }
  ## trivial cases
  if(nsim == 0) return(simulationresult(list()))
  if(kappa == 0 || mu == 0) {
    ## intensity is zero - patterns are empty
    empt <- ppp(window=W)
    if(saveparents) {
      attr(empt, "parents") <- list(x=numeric(0), y=numeric(0))
      attr(empt, "parentid") <- integer(0)
      attr(empt, "cost") <- 0
    }
    result <- rep(list(empt), nsim)
    return(simulationresult(result, nsim=nsim, drop=drop))
  }
  ## shift window to convenient origin
  oldW <- W
  oldcentre <- as.numeric(centroid.owin(Frame(oldW)))
  W <- shift(oldW, -oldcentre)
  ## enclose it in a disc
  rD <- with(vertices(Frame(W)), sqrt(max(x^2+y^2)))
  ## optimal inflation
  if(is.null(inflate)) {
    a <- if(mu == 0) 1 else (1 + (1-exp(-mu))/mu)
    b <- (rD^2)/(2*a*sigma^2)
    if(b <= 1) {
      inflate <- 1
    } else {
      delta <- 2 * sigma * sqrt(log(b)/2)
      inflate <- 1 + delta/rD
      inflate <- min(inflate, maxinflate)
    }
  }
  ## Prepare for C code
  storage.mode(kappa) <- "double"
  storage.mode(mu) <- "double"
  storage.mode(sigma) <- "double"
  storage.mode(rD) <- "double"
  storage.mode(inflate) <- "double"
  ##
  resultlist <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    ## call C code
    if(saveparents) {
      z <- .Call(SR_rthomasAll,
                 kappa, mu, sigma, rD, inflate,
                 PACKAGE="spatstat.random")
    } else {
      z <- .Call(SR_rthomasOff,
                 kappa, mu, sigma, rD, inflate,
                 PACKAGE="spatstat.random")
    }
    ## unpack
    xo <- z[[1]]
    yo <- z[[2]]
    if(saveparents) {
      xp <- z[[3]]
      yp <- z[[4]]
      parentid <- z[[5]]
    }
    ## shift back to original window
    xo <- xo + oldcentre[1L]
    yo <- yo + oldcentre[2L]
    if(saveparents) {
      xp <- xp + oldcentre[1L]
      yp <- yp + oldcentre[2L]
    }
    ## restrict to original window
    retain <- inside.owin(xo, yo, oldW)
    if(!all(retain)) {
      xo <- xo[retain]
      yo <- yo[retain]
      if(saveparents) {
        parentid <- parentid[retain]
        retainedparents <- sort(unique(parentid))
        parentid <- match(parentid, retainedparents)
        xp <- xp[retainedparents]
        yp <- yp[retainedparents]
      }
    }
    ## save as point pattern
    Y <- ppp(xo, yo, window=oldW, check=FALSE)
    if(saveparents) {
      attr(Y, "parents") <- list(x = xp, y = yp)
      attr(Y, "parentid") <- parentid
      attr(Y, "cost") <- length(xo) + length(xp)
    }
    resultlist[[isim]] <- Y
  }
  result <- simulationresult(resultlist, nsim, drop=drop)
  return(result)
}



rThomas <- local({

  ## random displacements
  gaus <- function(n, sigma) {
    matrix(rnorm(2 * n, mean=0, sd=sigma), ncol=2)
  }

  ## main function
  rThomas <-
    function(kappa, scale, mu, win = square(1),
             nsim=1, drop=TRUE,
             ...,
             algorithm=c("BKBC", "naive"),
             nonempty=TRUE, 
             poisthresh=1e-6,
             expand = 4*scale,
             saveparents=FALSE, saveLambda=FALSE,
             kappamax=NULL, mumax=NULL, sigma) {
      ## modified Thomas process
      ## Poisson(mu) number of offspring
      ## at isotropic Normal(0,sigma^2) displacements from parent
      
      check.1.integer(nsim)
      stopifnot(nsim >= 0)
      if(nsim == 0) return(simulationresult(list()))

      ## Catch old scale syntax (sigma)
      if((missing(scale) || is.null(scale)) && !missing(sigma)) {
        ## message("Argument 'sigma' is deprecated; it has been replaced by 'scale'")
        scale <- sigma
      }
      check.1.real(scale)
      stopifnot(scale > 0)

      ## determine the effective maximum radius of clusters
      ## (for the naive algorithm, or when kappa is not constant)
      if(missing(expand))
          expand <- clusterradius("Thomas", scale = scale, ...)

      #' validate 'kappa' and 'mu'
      km <- validate.kappa.mu(kappa, mu, kappamax, mumax,
                              win, expand, ...,
                              context="In rThomas")
      kappamax <- km[["kappamax"]]
      mumax    <- km[["mumax"]]

      ## detect trivial case where patterns are empty
      if(kappamax == 0 || mumax == 0) {
        empt <- ppp(window=win)
        if(saveparents) {
          attr(empt, "parents") <- list(x=numeric(0), y=numeric(0))
          attr(empt, "parentid") <- integer(0)
          attr(empt, "cost") <- 0
        }
        if(saveLambda) 
          attr(empt, "Lambda") <- as.im(0, W=win)
        result <- rep(list(empt), nsim)
        return(simulationresult(result, nsim=nsim, drop=drop))
      }

      #' determine algorithm
      algorithm <- match.arg(algorithm)
      do.parents <- saveparents || saveLambda || !is.numeric(kappa)
      do.hybrid <- (algorithm == "BKBC") && nonempty 

      if(do.hybrid) {
        ## ........ Fast algorithm (BKBC) .................................
        ## run BKBC algorithm for stationary model
        result <- rThomasHom(kappamax, mumax, sigma=scale,
                             W=win, ..., nsim=nsim, drop=FALSE,
                             saveparents=do.parents)

        ## thin 
        if(!is.numeric(kappa))
          result <- solapply(result, thinParents,
                             P=kappa, Pmax=kappamax)

        if(!is.numeric(mu)) 
          result <- solapply(result, rthin,
                             P=mu, Pmax=mumax,
                             na.zero=TRUE, fatal=FALSE)
      } else {
        ## .......... Slower algorithm ('naive') ..........................
        ## trap case of large clusters, close to Poisson
        if(is.numeric(kappa) && 1/(4*pi * kappa * scale^2) < poisthresh) {
          if(is.function(mu)) mu <- as.im(mu, W=win, ...)
          kapmu <- kappa * mu
          result <- rpoispp(kapmu, win=win, nsim=nsim, drop=drop, warnwin=FALSE)
          result <- fakeNeyScot(result, kapmu, win, saveLambda, saveparents)
          return(result)
        }

        result <- rNeymanScott(kappa=kappa,
                               expand=expand,
                               rcluster=list(mu, gaus),
                               win=win,
                               sigma=scale, # formal argument of 'gaus'
                               nsim=nsim, drop=FALSE,
                               nonempty=nonempty,
                               saveparents = do.parents,
                               kappamax=kappamax, mumax=mumax)
      }
      
      if(saveLambda){
        BW <- Frame(win)
        for(i in 1:nsim) {
          parents <- attr(result[[i]], "parents")
          BX <- boundingbox(BW, bounding.box.xy(parents))
          parents <- as.ppp(parents, W=BX, check=FALSE)
          Lambda <- clusterfield("Thomas", parents, scale=scale, mu=mu, ...)
          attr(result[[i]], "Lambda") <- Lambda[win, drop=FALSE]
        }
      }
      return(simulationresult(result, nsim, drop))
    }

  rThomas
})

