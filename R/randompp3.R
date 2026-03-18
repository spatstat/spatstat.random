#'
#'   randompp3.R
#'
#'   $Revision: 1.4 $ $Date: 2026/03/18 07:46:31 $
#'

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
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for(i in seq_len(nsim)) {
    x <- with(dd, runif(n, min=xrange[1], max=xrange[2]))
    y <- with(dd, runif(n, min=yrange[1], max=yrange[2]))
    z <- with(dd, runif(n, min=zrange[1], max=zrange[2]))
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  if(nsim > 0) names(result) <- paste("Simulation", seq_len(nsim))
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
    dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
    result <- vector(mode="list", length=nsim)
    probrange <- numeric(0)
    for(i in seq_len(nsim)) {
      ni <- np[i]
      x <- with(dd, runif(ni, min=xrange[1], max=xrange[2]))
      y <- with(dd, runif(ni, min=yrange[1], max=yrange[2]))
      z <- with(dd, runif(ni, min=zrange[1], max=zrange[2]))
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

