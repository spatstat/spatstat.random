#'
#'   randompp3.R
#'
#'   Random point patterns in 3 dimensions
#'
#'   runifpoint3    i.i.d. uniformly distributed
#'   rpoint3        i.i.d.
#'   rpoispp3       Poisson (homogeneous or inhomogeneous)
#' 
#'   $Revision: 1.8 $ $Date: 2026/03/19 06:40:37 $
#'
#'   Copyright (c) Adrian Baddeley, Ege Ruba and Rolf Turner 1993-2026
#'   GNU Public Licence (>= 2.0)

runifpoint3 <- function(n, domain=box3(), ..., nsim=1, drop=TRUE, ex=NULL) {
  if(!is.null(ex)) {
    ## make a pattern like 'ex'
    stopifnot(is.pp3(ex))
    if(missing(n) || is.null(n)) n <- npoints(ex)
    if(missing(domain) || is.null(domain)) domain <- domain(ex)
  }
  domain <- as.box3(domain)
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  result <- vector(mode="list", length=nsim)
  xr <- domain$xrange
  yr <- domain$yrange
  zr <- domain$zrange
  for(i in seq_len(nsim)) {
    x <- runif(n, min=xr[1], max=xr[2])
    y <- runif(n, min=yr[1], max=yr[2])
    z <- runif(n, min=zr[1], max=zr[2])
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  if(nsim > 0) names(result) <- paste("Simulation", seq_len(nsim))
  return(result)
}

rpoint3 <- function(n, f, fmax=NULL, domain=box3(), ...,
                    nsim=1, drop=TRUE, ex=NULL, giveup=1000) {
  if(missing(f) || is.null(f) || (is.numeric(f) && length(f) == 1)) {
    return(runifpoint3(n, domain=domain, nsim=nsim, drop=drop, ex=ex))
  }

  stopifnot(is.function(f))
  if(is.null(fmax)) stop("Argument fmax is required", call.=FALSE)
  check.1.real(fmax)
  stopifnot(fmax >= 0)
  
  if(!is.null(ex)) {
    ## make a pattern like 'ex'
    stopifnot(is.pp3(ex))
    if(missing(n) || is.null(n)) n <- npoints(ex)
    if(missing(domain) || is.null(domain)) domain <- domain(ex)
  }
  domain <- as.box3(domain)
  check.1.integer(nsim)
  stopifnot(nsim >= 0)

  if(fmax == 0) {
    ## empty patterns
    result <- rep(list(NAobject("pp3")), nsim)
    result <- simulationresult(result, drop=drop)
    return(result)
  }

  ## prepare to simulate
  ## first generate all point coordinates
  ntotal <- nsim * n
  xx <- yy <- zz <- numeric(0)
  nn <- 0
  ntries <- 0
  xr <- domain$xrange
  yr <- domain$yrange
  zr <- domain$zrange
  prange <- c(0,1)
  while(nn < ntotal && ntries < giveup) {
    ntries <- ntries + 1
    xprop <- runif(ntotal, min=xr[1], max=xr[2])
    yprop <- runif(ntotal, min=yr[1], max=yr[2])
    zprop <- runif(ntotal, min=zr[1], max=zr[2])
    ## rejection method
    p <- f(xprop,yprop,zprop, ...)/fmax
    if(!all(is.finite(p)))
      stop("Infinite, NA or NaN values in f(x,y,z)", call.=FALSE)
    prange <- range(p, prange)
    retain <- (runif(ntotal) <= p)
    xx <- c(xx, xprop[retain])
    yy <- c(yy, yprop[retain])
    zz <- c(zz, zprop[retain])
    nn <- length(xx)
  }
  if(prange[1L] < 0)
    warning("Negative values of f were encountered", call.=FALSE)
  if(prange[2L] > 1)
    warning("Some values of f exceeded fmax", call.=FALSE)
  ## now divide coordinate vectors into separate point patterns
  result <- vector(mode="list", length=nsim)
  sn <- seq_len(n)
  for(i in seq_len(nsim)) {
    if(nn >= n) {
      result[[i]] <- pp3(xx[sn], yy[sn], zz[sn], domain)
      xx <- xx[-sn]
      yy <- yy[-sn]
      zz <- zz[-sn]
      nn <- length(xx)
    } else {
      warning("Exited after", giveup, "tries with only",
              i - 1L, "patterns generated")
      result[i:nsim] <- rep(list(NAobject("pp3")), nsim-i+1L)
      break
    } 
  }
  result <- simulationresult(result, drop=drop)
  return(result)
}

rpoispp3 <- function(lambda, domain=box3(), ...,
                     nsim=1, drop=TRUE, ex=NULL, lmax=NULL) {
  if(!is.null(ex)) {
    stopifnot(is.pp3(ex))
    if(missing(lambda) || is.null(lambda)) lambda <- intensity(unmark(ex))
    if(missing(domain) || is.null(domain)) domain <- domain(ex)
  }
  domain <- as.box3(domain)
  v <- volume(domain)
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  if(is.numeric(lambda) && length(as.numeric(lambda)) == 1) {
    lmax <- as.numeric(lambda)
  } else if(is.function(lambda)) {
    if(is.null(lmax))
      stop("lmax is required when lambda is a function", call.=FALSE)
    check.1.real(lmax)
  } else {    
    stop("lambda must be a single numeric value, or a function(x,y,z")
  }
  inhom <- is.function(lambda)
  if(lmax <= 0) {
    ## zero intensity => empty patterns
    none <- numeric(0)
    emptypattern <- pp3(none, none, none, domain)
    result <- rep(list(emptypattern), nsim)
  } else {
    ## generate homogeneous Poisson
    np <- rpois(nsim, lmax * v)
    dom <- as.list(domain)[c("xrange", "yrange", "zrange")]
    result <- vector(mode="list", length=nsim)
    probrange <- numeric(0)
    for(i in seq_len(nsim)) {
      ni <- np[i]
      x <- with(dom, runif(ni, min=xrange[1], max=xrange[2]))
      y <- with(dom, runif(ni, min=yrange[1], max=yrange[2]))
      z <- with(dom, runif(ni, min=zrange[1], max=zrange[2]))
      X.i <- pp3(x,y,z,domain)
      if(inhom && ni > 0) {
        ## Lewis-Shedler thinning 
        p <- lambda(x,y,z)/lmax
        probrange <- range(p, probrange)
        retain <- (runif(ni) <= p)
        X.i <- X.i[retain]
      }
      result[[i]] <- X.i
    }
    if(inhom && length(probrange)) {
      if(probrange[1L] < 0)
        warning("Negative values of lambda(x,y,z) occurred", call.=FALSE)
      if(probrange[2L] > 1)
        warning("Some values of lambda(x,y,z) exceeded lmax", call.=FALSE)
    }
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  if(nsim > 0) names(result) <- paste("Simulation", seq_len(nsim))
  return(result)
}

