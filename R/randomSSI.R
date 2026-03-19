#'
#'   randomSSI.R
#'
#'   Simple sequential inhibition
#'
#'   $Revision: 1.4 $ $Date: 2026/03/19 06:11:55 $
#'
#'    Copyright (c) Adrian Baddeley, Ege Rubak and Rolf Turner 1994-2026
#'    GNU Public Licence (>= 2.0)


rSSI <- function(r, n=Inf, win = square(1), 
                 giveup = 1000, x.init=NULL, ...,
                 f=NULL, fmax=NULL,
                 nsim=1, drop=TRUE, verbose=TRUE)
{
  win.given <- !missing(win) && !is.null(win)
  stopifnot(is.numeric(r) && length(r) == 1 && r >= 0)
  stopifnot(is.numeric(n) && length(n) == 1 && n >= 0)
  must.reach.n <- is.finite(n)
  if(!missing(nsim)) {
    check.1.integer(nsim)
    stopifnot(nsim >= 1)
  }
  ## 'verbose' applies only when nsim > 1
  verbose <- verbose && (nsim > 1)
  #' resolve spatial domain
  if(win.given) {
    ## ensure 'win' is a window, a box, or a linear network
    if(!inherits(win, c("owin", "box3", "boxx", "linnet"))) {
      winwin <- try(as.owin(win), silent = TRUE)
      if(inherits(winwin, "try-error"))
        winwin <- try(as.boxx(win), silent = TRUE)
      if(inherits(winwin, "try-error") && requireNamespace("spatstat.linnet"))
        winwin <- try(spatstat.linnet::as.linnet(win), silent=TRUE)
      if(inherits(winwin, "try-error"))
        stop(paste("Could not coerce argument", sQuote("win"), 
                   "to a window, a box, or a linear network"),
             call.=FALSE)
      win <- winwin
    }
  } else if(!is.null(x.init)) {
    if(!inherits(x.init, c("ppp", "lpp", "pp3", "ppx")))
      stop(paste("x.init should be a point pattern of class",
                 commasep(dQuote(c("ppp", "pp3", "ppx", "lpp"), "or"))),
           call.=FALSE)
    win <- domain(x.init)
  } else if(inherits(f, c("im", "funxy", "linfun"))) {
    win <- domain(f)
  } else {
    win <- square(1)
  }

  rSSIengine(win, r=r, n=n, giveup = giveup,
             x.init=x.init, f=f, fmax=fmax,
             nsim=nsim, drop=drop, verbose=verbose, ...)
}

#' Algorithm depends on type of domain

rSSIengine <- function(win, ...) {
  UseMethod("rSSIengine")
}

rSSIengine.owin <- function(win, ..., r,  n=Inf, 
                            giveup = 1000, x.init=NULL,
                            f=NULL, fmax=NULL,
                            nsim=1, drop=TRUE, verbose=TRUE, warn=TRUE) {
  must.reach.n <- is.finite(n)
  if(!is.null(f)) {
    stopifnot(is.numeric(f) || is.im(f) || is.function(f))
    if(is.null(fmax) && !is.numeric(f))
      fmax <- if(is.im(f)) max(f) else max(as.im(f, win))
  }
  ## determine initial state
  if(is.null(x.init)) {
    ## start with empty pattern in specified window
    win <- as.owin(win)
    x.init <- ppp(numeric(0),numeric(0), window=win)
  } else {
    ## start with specified pattern
    stopifnot(is.ppp(x.init))
    ## check compatibility of windows
    if(!identical(win, as.owin(x.init)))
      warning(paste("Argument", sQuote("win"),
                    "is not the same as the window of", sQuote("x.init")))
    x.init.new <- x.init[win]
    if(npoints(x.init.new) == 0)
      stop(paste("No points of x.init lie inside the specified window",
                 sQuote("win")))
    nlost <- npoints(x.init) - npoints(x.init.new)
    if(nlost > 0) 
      warning(paste(nlost, "out of",
                    npoints(x.init), "points of the pattern x.init",
                    "lay outside the specified window",
                    sQuote("win")))
    x.init <- x.init.new
    if(is.finite(n)) {
      if(n < npoints(x.init))
        stop(paste("x.init contains", npoints(x.init), "points",
                   "but a pattern containing only n =", n, "points", 
                   "is required"))
      if(n == npoints(x.init)) {
        warning(paste("Initial state x.init already contains", n, "points;",
                      "no further points were added"))
        result <- rep(list(x.init), nsim)
        result <- simulationresult(result, nsim, drop)
        return(result)
      }
    }
  }
  #' assess packing density
  r2 <- r^2
  winArea <- area(win)
  discarea <- pi * r2/4   # area of disc of diameter r
  arearatio <- winArea/discarea
  nmelt       <- floor( arearatio )
  nclosepack  <- floor( arearatio * pi * sqrt(3)/6 )
  nrandompack <- floor( arearatio * 0.88 )
  nsquarepack <- floor( arearatio * pi/4 )
  if(warn && is.finite(n)) {
    if(n > nmelt) {
      warning(paste("Window is far too small to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("absolute maximum number is", nmelt))))
    } else if(n > nclosepack) {
      warning(paste("Window is too small to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("deterministic close packing limit is",
                                nclosepack))))
    } else if(n > nrandompack) {
      warning(paste("Window is very unlikely to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("random close packing limit is",
                                nrandompack))))
    } else if(n > nsquarepack) {
      warning(paste("Window is unlikely to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("deterministic square-packing limit is",
                                nsquarepack))))
    }
  }

  #' start simulation 		    
  pstate <- list()
  if(verbose) splat("Generating", nsim, "realisations...")
  result <- vector(mode="list", length=nsim)
  for(isim in seq_len(nsim)) {
    if(verbose) pstate <- progressreport(isim, nsim, state=pstate)
    ## Simple Sequential Inhibition process
    ## fixed number of points
    xx <- coords(x.init)$x
    yy <- coords(x.init)$y
    nn <- npoints(x.init)
    ## Naive implementation, proposals are uniform
    xprop <- yprop <- numeric(0)
    nblock <- if(is.finite(n)) n else min(1024, nmelt)
    ntries <- 0
    while(ntries < giveup) {
      ntries <- ntries + 1
      if(length(xprop) == 0) {
        ## generate some more proposal points
        prop <- if(is.null(f)) runifpoint(nblock, win) else
                               rpoint(nblock, f, fmax, win)
        xprop <- coords(prop)$x
        yprop <- coords(prop)$y
      }
      ## extract next proposal
      xnew <- xprop[1L]
      ynew <- yprop[1L]
      xprop <- xprop[-1L]
      yprop <- yprop[-1L]
      ## check hard core constraint
      dx <- xnew - xx
      dy <- ynew - yy
      if(!any(dx^2 + dy^2 <= r2)) {
        xx <- c(xx, xnew)
        yy <- c(yy, ynew)
        nn <- nn + 1L
        ntries <- 0
      }
      if(nn >= n)
        break
    }
    if(must.reach.n && nn < n)
      warning(paste("Gave up after", giveup,
                    "attempts with only", nn, "points placed out of", n))
    X <- ppp(xx, yy, window=win, check=FALSE)
    result[[isim]] <- X
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}


rSSIengine.box3 <- function(win, ..., r,  n=Inf, 
                            giveup = 1000, x.init=NULL,
                            f=NULL, fmax=NULL,
                            nsim=1, drop=TRUE, verbose=TRUE, warn = TRUE) {
  must.reach.n <- is.finite(n)
  if(!is.null(f)) {
    stopifnot(is.numeric(f) || is.function(f))
    if(is.null(fmax) && is.function(f))
      stop("Argument fmax is required when f is a function", call.=FALSE)
  }
  ## determine initial state
  if(is.null(x.init)) {
    ## start with empty pattern in specified window
    x.init <- pp3(numeric(0),numeric(0), numeric(0), win)
  } else {
    ## start with specified pattern
    stopifnot(is.pp3(x.init))
    ## check compatibility of domains
    if(!identical(win, domain(x.init)))
      warning(paste("Argument", sQuote("win"),
                    "is not the same as the domain of", sQuote("x.init")))
    x.init.new <- x.init[win]
    if(npoints(x.init.new) == 0)
      stop(paste("No points of x.init lie inside the specified window",
                 sQuote("win")))
    nlost <- npoints(x.init) - npoints(x.init.new)
    if(nlost > 0) 
      warning(paste(nlost, "out of",
                    npoints(x.init), "points of the pattern x.init",
                    "lay outside the specified window",
                    sQuote("win")))
    x.init <- x.init.new
    if(is.finite(n)) {
      if(n < npoints(x.init))
        stop(paste("x.init contains", npoints(x.init), "points",
                   "but a pattern containing only n =", n, "points", 
                   "is required"))
      if(n == npoints(x.init)) {
        warning(paste("Initial state x.init already contains", n, "points;",
                      "no further points were added"))
        result <- rep(list(x.init), nsim)
        result <- simulationresult(result, nsim, drop)
        return(result)
      }
    }
  }
  #' validate radius and 'n' 
  r2 <- r^2
  winVol <- volume(win)
  sphVol <- (pi/6) * r^3
  volratio <- winVol/sphVol
  nmelt       <- floor( volratio )
  nclosepack  <- floor( volratio * pi * sqrt(2)/6 )
  nrandompack <- floor( volratio * 0.64 )
  ncubepack   <- floor( volratio * pi/8 )
  if(warn && is.finite(n) && n > ncubepack) {
    if(n > nmelt) {
      warning(paste("Window is far too small to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("absolute maximum number is", nmelt))))
    } else if(n > nclosepack) {
      warning(paste("Window is probably too small to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("deterministic close packing limit is",
                                nclosepack))))
    } else if(n > nrandompack) {
      warning(paste("Algorithm is very unlikely to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("random close packing limit is",
                                nclosepack))))
    } else if(n > ncubepack) {
      warning(paste("Algorithm may struggle to fit", n, "points",
                    "at minimum separation", r,
                    paren(paste("deterministic cubic packing limit is",
                                ncubepack))))
    }
  }

  #' start simulation 		    
  pstate <- list()
  if(verbose) splat("Generating", nsim, "realisations...")
  result <- vector(mode="list", length=nsim)
  for(isim in seq_len(nsim)) {
    if(verbose) pstate <- progressreport(isim, nsim, state=pstate)
    ## Simple Sequential Inhibition process
    ## fixed number of points
    co <- coords(x.init)
    xx <- co$x
    yy <- co$y
    zz <- co$z
    nn <- npoints(x.init)
    ## Naive implementation, proposals are uniform
    xprop <- yprop <- numeric(0)
    nblock <- if(is.finite(n)) n else min(1024, nmelt)
    ntries <- 0
    while(ntries < giveup) {
      ntries <- ntries + 1
      if(length(xprop) == 0) {
        ## generate some more proposal points
        prop <- runifpoint3(nblock, win)
        if(is.function(f)) {
          ## randomly thin
          p <- with(coords(prop), f(x, y, z, ...)/fmax)
          retain <- runif(npoints(prop)) < p
          prop <- prop[retain]
        }
        cop <- coords(prop)
        xprop <- cop$x
        yprop <- cop$y
        zprop <- cop$z
      }
      if(length(xprop) > 0) {
        ## extract next proposal
        xnew <- xprop[1L]
        ynew <- yprop[1L]
        znew <- zprop[1L]
        xprop <- xprop[-1L]
        yprop <- yprop[-1L]
        zprop <- zprop[-1L]
        ## check hard core constraint
        dx <- xnew - xx
        dy <- ynew - yy
        dz <- znew - zz
        if(!any(dx^2 + dy^2 + dz^2 <= r2)) {
          xx <- c(xx, xnew)
          yy <- c(yy, ynew)
          zz <- c(zz, znew)
          nn <- nn + 1L
          ntries <- 0
        }
        if(nn >= n)
          break
      }
    }
    ## end of rejection method loop for simulation 'isim'
    if(must.reach.n && nn < n) {
      warning(paste("Gave up after", giveup,
                    "attempts with only", nn, "points placed out of", n))
      X <- NAobject("pp3")
    } else {
      X <- pp3(xx, yy, zz, win)
    }
    result[[isim]] <- X
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

