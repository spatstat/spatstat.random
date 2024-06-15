#
#  Perfect Simulation 
#
#  $Revision: 1.27 $ $Date: 2024/06/09 00:22:17 $
#
#  rStrauss
#  rHardcore
#  rStraussHard
#  rDiggleGratton
#  rDGS
#  rPenttinen

rStrauss <- function(beta, gamma=1, R=0, W=owin(), expand=TRUE,
                     nsim=1, drop=TRUE) {

  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(gamma)
  check.1.real(R)

  check.finite(beta, xname="beta")
  check.finite(gamma, xname="gamma")
  check.finite(R, xname="R")
  
  stopifnot(beta > 0)
  stopifnot(gamma >= 0)
  stopifnot(gamma <= 1)
  stopifnot(R >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in seq_len(nsim)) {
    storage.mode(beta) <- storage.mode(gamma) <- storage.mode(R) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call(SR_PerfectStrauss,
               beta,
               gamma,
               R,
               xrange,
               yrange,
               PACKAGE="spatstat.random")

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    times <- c(start=z[[4]], end=z[[5]])
    
    if(nout<0)
      stop("internal error: copying failed in PerfectStrauss")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]
    attr(P, "times") <- times

    if(nsim == 1 && drop) return(P)
    
    result[[i]] <- P
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

#  Perfect Simulation of Hardcore process

rHardcore <- function(beta, R=0, W=owin(), expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(R)

  check.finite(beta, xname="beta")
  check.finite(R, xname="R")

  stopifnot(beta > 0)
  stopifnot(R    >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in seq_len(nsim)) {
    storage.mode(beta) <- storage.mode(R) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call(SR_PerfectHardcore,
               beta,
               R,
               xrange,
               yrange,
               PACKAGE="spatstat.random")

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    
    if(nout<0)
      stop("internal error: copying failed in PerfectHardcore")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

#
#  Perfect simulation of hybrid Strauss-Hardcore
#        provided gamma <= 1
#

rStraussHard <- function(beta, gamma=1, R=0, H=0, W=owin(),
                         expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(gamma)
  check.1.real(R)
  check.1.real(H)

  check.finite(beta, xname="beta")
  check.finite(gamma, xname="gamma")
  check.finite(R, xname="R")
  check.finite(H, xname="H")
  
  stopifnot(beta > 0)
  stopifnot(gamma >= 0)
  if(gamma > 1)
    stop("Sorry, perfect simulation is only implemented for gamma <= 1")
  stopifnot(R >= 0)
  stopifnot(H >= 0)
  stopifnot(H <= R)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in seq_len(nsim)) {
    storage.mode(beta) <- storage.mode(gamma) <-
      storage.mode(R) <- storage.mode(H) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call(SR_PerfectStraussHard,
               beta,
               gamma,
               R,
               H,
               xrange,
               yrange,
               PACKAGE="spatstat.random")

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]

    if(nout<0)
      stop("internal error: copying failed in PerfectStraussHard")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

#
#  Perfect Simulation of Diggle-Gratton process
#

rDiggleGratton <- function(beta, delta, rho, kappa=1, W=owin(),
                           expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(delta)
  check.1.real(rho)
  check.1.real(kappa)

  check.finite(beta, xname="beta")
  check.finite(delta, xname="delta")
  check.finite(rho, xname="rho")
  check.finite(kappa, xname="kappa")

  stopifnot(beta > 0)
  stopifnot(delta >= 0)
  stopifnot(rho   >= 0)
  stopifnot(delta <= rho)
  stopifnot(kappa >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*rho))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in seq_len(nsim)) {
    storage.mode(beta) <- "double"
    storage.mode(delta) <- storage.mode(rho) <- storage.mode(kappa) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call(SR_PerfectDiggleGratton,
               beta,
               delta,
               rho,
               kappa,
               xrange,
               yrange,
               PACKAGE="spatstat.random")

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]

    if(nout<0)
      stop("internal error: copying failed in PerfectDiggleGratton")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}


#
#  Perfect Simulation of Diggle-Gates-Stibbard process
#

rDGS <- function(beta, rho, W=owin(), expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(rho)

  check.finite(beta, xname="beta")
  check.finite(rho, xname="rho")

  stopifnot(beta > 0)
  stopifnot(rho  >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*rho))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in seq_len(nsim)) {
    storage.mode(beta) <- "double"
    storage.mode(rho) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call(SR_PerfectDGS,
               beta,
               rho,
               xrange,
               yrange,
               PACKAGE="spatstat.random")

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    
    if(nout<0)
      stop("internal error: copying failed in PerfectDGS")
    
    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}


#
#  Perfect Simulation of Penttinen process
#

rPenttinen <- function(beta, gamma=1, R, W=owin(),
                       expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(gamma)
  check.1.real(R)

  check.finite(beta, xname="beta")
  check.finite(gamma, xname="gamma")
  check.finite(R, xname="R")

  stopifnot(beta > 0)
  stopifnot(gamma >= 0)
  stopifnot(gamma <= 1)
  stopifnot(R >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in seq_len(nsim)) {
    storage.mode(beta) <- storage.mode(gamma) <- storage.mode(R) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call(SR_PerfectPenttinen,
               beta,
               gamma,
               R,
               xrange,
               yrange,
               PACKAGE="spatstat.random")

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    
    if(nout<0)
      stop("internal error: copying failed in PerfectPenttinen")
    
    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}


## .......  utilities .................................

expandwinPerfect <- function(W, expand, amount) {
  ## expand 'W' if expand=TRUE according to default 'amount'
  ## or expand 'W' using rmhexpand(expand)
  if(!is.logical(expand)) {
    amount <- rmhexpand(expand)
    expand <- TRUE
  }
  changed <- FALSE
  if(expand) {
    W <- expand.owin(W, amount)
    changed <- TRUE
  }
  if(!is.rectangle(W)) {
    W <- as.rectangle(W)
    changed <- TRUE
    warning(paste("Simulation will be performed in the containing rectangle",
                  "and clipped to the original window."),
            call.=FALSE)
  }
  attr(W, "changed") <- changed
  return(W)
}
