#'
#'    rmatclust.R
#'
#'   $Revision: 1.3 $ $Date: 2023/01/06 11:37:11 $
#'
#'   Simulation of Matern cluster process
#'   naive algorithm or BKBC algorithm
#'
#'   rMatClustHom    Interface to C code for stationary case (BKBC)
#'   rMatClust       General case (naive or BKBC)
#' 
#'   Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022-2023
#'   Licence: GNU Public Licence >= 2

rMatClustHom <- function(kappa, mu, R, W=unit.square(), ...,
                         nsim=1, drop=TRUE, inflate=NULL,
                         saveparents=FALSE) {
  check.1.real(kappa) && check.finite(kappa)
  check.1.real(mu) && check.finite(mu)
  check.1.real(R) && check.finite(R)
  if(!is.null(inflate)) {
    check.1.real(inflate) && check.finite(inflate)
    stopifnot(inflate >= 1)
  }
  check.1.integer(nsim)
  stopifnot(kappa >= 0)
  stopifnot(mu >= 0)
  stopifnot(R > 0)
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
    rE <- if(R < rD) (rD + R) else rD
    inflate <- rE/rD
  }
  ## Prepare for C code
  storage.mode(kappa) <- "double"
  storage.mode(mu) <- "double"
  storage.mode(R) <- "double"
  storage.mode(rD) <- "double"
  storage.mode(inflate) <- "double"
  ##
  resultlist <- vector(mode="list", length=nsim)
  for(isim in seq_len(nsim)) {
    ## call C code
    if(saveparents) {
      z <- .Call(SR_rmatclusAll,
                 kappa, mu, R, rD, inflate,
                 PACKAGE="spatstat.random")
    } else {
      z <- .Call(SR_rmatclusOff,
                 kappa, mu, R, rD, inflate,
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


rMatClust <- local({
  
  ## like runifdisc but returns only the coordinates
  rundisk <- function(n, radius) {
    R <- radius * sqrt(runif(n, min=0, max=1))
    Theta <- runif(n, min=0, max=2*pi)
    cbind(R * cos(Theta), R * sin(Theta))
  }

  rMatClust <- 
  function(kappa, scale, mu, win = square(1),
           nsim=1, drop=TRUE, 
           ...,
           algorithm=c("BKBC", "naive"),
           nonempty=TRUE, 
           poisthresh=1e-6,
           saveparents=FALSE, saveLambda=FALSE,
           kappamax=NULL, mumax=NULL) {
    ## Matern Cluster Process
    ## Poisson (mu) number of offspring, uniform inside disc

    check.1.integer(nsim) && check.finite(nsim)
    stopifnot(nsim >= 0)
    if(nsim == 0) return(simulationresult(list()))

    ## Catch old scale syntax (r)
    if(missing(scale)) scale <- list(...)$r
    check.1.real(scale)
    stopifnot(scale > 0)

    #' validate 'kappa' and 'mu'
    km <- validate.kappa.mu(kappa, mu, kappamax, mumax,
                            win, scale, ..., 
                            context="In rMatClust")
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
      result <- rMatClustHom(kappamax, mumax, scale, W=win, ...,
                             nsim=nsim, drop=FALSE,
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
      if(is.numeric(kappa) && 1/(pi * kappa * scale^2) < poisthresh) {
        if(is.function(mu)) mu <- as.im(mu, W=win, ...)
        kapmu <- kappa * mu
        result <- rpoispp(kapmu, win=win, nsim=nsim, drop=drop, warnwin=FALSE)
        result <- fakeNeyScot(result, kapmu, win, saveLambda, saveparents)
        return(result)
      }

      result <- rNeymanScott(kappa=kappa,
                             expand=scale,
                             rcluster=list(mu, rundisk),
                             win=win,
                             radius=scale, # formal argument of 'rundisk'
                             nsim=nsim, drop=FALSE,
                             nonempty=nonempty,
                             saveparents = do.parents,
                             kappamax=kappamax, mumax=mumax)
    }
    
    if(saveLambda){
      B <- grow.rectangle(Frame(win), scale)
      for(i in 1:nsim) {
        parents <- attr(result[[i]], "parents")
        parents <- as.ppp(parents, W=B, check=FALSE)
        Lambda <- clusterfield("MatClust", parents, scale=scale, mu=mu, ...)
        attr(result[[i]], "Lambda") <- Lambda[win, drop=FALSE]
      }
    }
    return(simulationresult(result, nsim, drop))
  }

  rMatClust
})