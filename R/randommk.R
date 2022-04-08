#
#
#   randommk.R
#
#   Random generators for MULTITYPE point processes
#
#   $Revision: 1.41 $   $Date: 2022/04/08 06:25:53 $
#
#   rmpoispp()   random marked Poisson pp
#   rmpoint()    n independent random marked points
#   rmpoint.I.allim()  ... internal
#   rpoint.multi()   temporary wrapper 
#
rmpoispp <- local({

  ## Argument checking
  is.numvector <- function(x) {is.numeric(x) && is.vector(x)}
  is.constant <- function(x) {is.numvector(x) && length(x) == 1}
  checkone <- function(x) {
    if(is.constant(x)) {
      if(x >= 0) return(TRUE) else stop("Intensity is negative!")
    }
    return(is.function(x) || is.im(x))
  }

  ## Ensure that m can be passed as a single value to function(x,y,m,...)
  slice.fun <- function(x,y,fun,mvalue, ...) {
    m <- if(length(mvalue) == 1) rep.int(mvalue, length(x)) else mvalue
    result <- fun(x,y,m, ...)
    return(result)
  }

  ## Main function
  rmpoispp <- 
    function(lambda, lmax=NULL, win = owin(c(0,1),c(0,1)),
             types, ..., nsim=1, drop=TRUE, warnwin=!missing(win)) {
      ## arguments:
      ##     lambda  intensity:
      ##                constant, function(x,y,m,...), image,
      ##                vector, list of function(x,y,...) or list of images
      ##
      ##     lmax     maximum possible value of lambda
      ##                constant, vector, or list
      ##
      ##     win     default observation window (of class 'owin')
      ##
      ##     types    possible types for multitype pattern
      ##    
      ##     ...     extra arguments passed to lambda()
      ##

      check.1.integer(nsim)
      stopifnot(nsim >= 0)
      if(nsim == 0) return(simulationresult(list()))

      if(missing(types)) types <- NULL
      force(warnwin)
      
      ## Validate arguments
      single.arg <- checkone(lambda)
      vector.arg <- !single.arg && is.numvector(lambda) 
      list.arg <- !single.arg && is.list(lambda)
      if(! (single.arg || vector.arg || list.arg))
        stop(paste("argument", sQuote("lambda"), "not understood"))
    
      if(list.arg && !all(unlist(lapply(lambda, checkone))))
        stop(paste("Each entry in the list",
                   sQuote("lambda"),
                   "must be either a constant, a function or an image"))
      if(vector.arg && any(lambda < 0))
        stop(paste("Some entries in the vector",
                   sQuote("lambda"), "are negative"))

      ## Determine & validate the set of possible types
      if(is.null(types)) {
        if(single.arg) {
          stop(paste(sQuote("types"), "must be given explicitly when",
                     sQuote("lambda"), "is a constant, a function or an image"))
        } else if(!is.null(nama <- names(lambda)) &&
                  sum(nzchar(nama)) == length(lambda)) {
          types <- nama
        } else {
          types <- seq_along(lambda)
        }
      } 
      ntypes <- length(types)
      if(!single.arg && (length(lambda) != ntypes))
        stop(paste("The lengths of", sQuote("lambda"),
                   "and", sQuote("types"), "do not match"))

      factortype <- factor(types, levels=types)

      ## Validate `lmax'
      if(! (is.null(lmax) || is.numvector(lmax) || is.list(lmax) ))
        stop(paste(sQuote("lmax"),
                   "should be a constant, a vector, a list or NULL"))
       
      ## coerce lmax to a vector, to save confusion
      if(is.null(lmax))
        maxes <- rep(NULL, ntypes)
      else if(is.numvector(lmax) && length(lmax) == 1)
        maxes <- rep.int(lmax, ntypes)
      else if(length(lmax) != ntypes)
        stop(paste("The length of",
                   sQuote("lmax"),
                   "does not match the number of possible types"))
      else if(is.list(lmax))
        maxes <- unlist(lmax)
      else maxes <- lmax

      ## coerce lambda to a list, to save confusion
      lam <- if(single.arg) rep(list(lambda), ntypes) else
             if(vector.arg) as.list(lambda) else lambda

      ## ------------- SIMULATE ----------------------------
      resultlist <- vector(mode="list", length=nsim)
      for(isim in seq_len(nsim)) {
        for(i in 1:ntypes) {
          if(single.arg && is.function(lambda)) {
            ## call f(x,y,m, ...)
            Y <- rpoispp(slice.fun, lmax=maxes[i], win=win,
                         fun=lambda, mvalue=types[i], ..., warnwin=warnwin)
          } else {
            ## call f(x,y, ...) or use other formats
            Y <- rpoispp(lam[[i]], lmax=maxes[i], win=win, ..., warnwin=warnwin)
          }
          Y <- Y %mark% factortype[i]
          X <- if(i == 1) Y else superimpose(X, Y, W=X$window, check=FALSE)
        }
        ## Randomly permute, just in case the order is important
        permu <- sample(X$n)
        X <- X[permu]
        ## save 
        resultlist[[isim]] <- X
      }
      return(simulationresult(resultlist, nsim, drop))
    }

  rmpoispp
})

## ------------------------------------------------------------------------

rmpoint <- local({

  ## argument validation
  is.numvector <- function(x) {is.numeric(x) && is.vector(x)}
  is.constant <- function(x) {is.numvector(x) && length(x) == 1}
  checkone <- function(x) {
    if(is.constant(x)) {
      if(x >= 0) return(TRUE) else stop("Intensity is negative!")
    }
    return(is.function(x) || is.im(x))
  }

  # integration..
  integratexy <- function(f, win, ...) {
    imag <- as.im(f, W=win, ...)
    integral.im(imag)
  }
  ## create a counterpart of f(x,y,m) that works when m is a single value
  funwithfixedmark <- function(xx, yy, ..., m, fun) {
    mm <- rep.int(m, length(xx))
    fun(xx, yy, mm, ...)
  }
  integratewithfixedmark <- function(m, fun, win, ...) {
    integratexy(funwithfixedmark, win=win, m=m, fun=fun, ...)
  }

  # Main function
  rmpoint <- function(n, f=1, fmax=NULL, 
                      win = unit.square(), 
                      types, ptypes, ...,
                      giveup = 1000, verbose = FALSE,
                      nsim = 1, drop=TRUE) {
    if(!is.numeric(n))
      stop("n must be a scalar or vector")
    if(any(ceiling(n) != floor(n)))
      stop("n must be an integer or integers")
    if(any(n < 0))
      stop("n must be non-negative")
    if(missing(types)) types <- NULL
    if(missing(ptypes)) ptypes <- NULL

    check.1.integer(nsim)
    stopifnot(nsim >= 0)

    ######  Trivial cases

    if(nsim == 0) return(simulationresult(list()))
    
    if(sum(n) == 0) {
      ## results are empty patterns
      nopoints <- ppp(x=numeric(0), y=numeric(0), window=win, check=FALSE)
      if(!is.null(types)) {
        nomarks <- factor(types[integer(0)], levels=types)
        nopoints <- nopoints %mark% nomarks
      }
      resultlist <- rep(list(nopoints), nsim)
      return(simulationresult(resultlist, nsim, drop))
    }

    ####### Validate arguments and determine type of simulation ##############
  
    Model <- if(length(n) == 1) {
      if(is.null(ptypes)) "I" else "II"
    } else "III"
  
    ##############  Validate f argument
    single.arg <- checkone(f)
    vector.arg <- !single.arg && is.numvector(f) 
    list.arg <- !single.arg && is.list(f)
    if(! (single.arg || vector.arg || list.arg))
      stop(paste("argument", sQuote("f"), "not understood"))
    
    if(list.arg && !all(unlist(lapply(f, checkone))))
      stop(paste("Each entry in the list", sQuote("f"),
                 "must be either a constant, a function or an image"))
    if(vector.arg && any(f < 0))
      stop(paste("Some entries in the vector",
                 sQuote("f"), "are negative"))
    
    ## cases where it's known that all types of points 
    ## have the same conditional density of location (x,y)
    const.density <- vector.arg ||
                     (list.arg && all(unlist(lapply(f, is.constant))))
    same.density <- const.density || (single.arg && !is.function(f))

    ################   Determine & validate the set of possible types
    if(is.null(types)) {
      if(single.arg && length(n) == 1)
        stop(paste(sQuote("types"), "must be given explicitly when",
                   sQuote("f"),
                   "is a single number, a function or an image and",
                   sQuote("n"), "is a single number"))
      else {
        basis <- if(single.arg) n else f
        if(!is.null(nama <- names(basis)) &&
           sum(nzchar(nama)) == length(basis)) {
          types <- nama
        } else {
          types <- seq_along(basis)
        }
      }
    }

    ntypes <- length(types)
    if(!single.arg && (length(f) != ntypes))
      stop(paste("The lengths of",
                 sQuote("f"), "and", sQuote("types"),
                 "do not match"))
    if(length(n) > 1 && ntypes != length(n))
      stop(paste("The lengths of",
                 sQuote("n"), "and", sQuote("types"),
                 "do not match"))

    factortype <- factor(types, levels=types)
    seqtypes <- seq_len(ntypes)
  
    #######################  Validate `fmax'
    if(! (is.null(fmax) || is.numvector(fmax) || is.list(fmax) ))
      stop(paste(sQuote("fmax"),
                 "should be a constant, a vector, a list or NULL"))
       
    ## coerce fmax to a vector, to save confusion
    if(is.null(fmax))
      maxes <- rep(NULL, ntypes)
    else if(is.constant(fmax))
      maxes <- rep.int(fmax, ntypes)
    else if(length(fmax) != ntypes)
      stop(paste("The length of", sQuote("fmax"),
                 "does not match the number of possible types"))
    else if(is.list(fmax))
      maxes <- unlist(fmax)
    else maxes <- fmax

    ## coerce f to a list, to save confusion
    flist <- if(single.arg) rep(list(f), ntypes) else
             if(vector.arg) as.list(f) else f

    #################### Finished validating arguments ########################

    #################### Handle special case ########################
    
    if(Model == "I" && !same.density && all(sapply(flist, is.im))) {
      resultlist <- rmpoint.I.allim(n, flist, types, nsim)
      return(simulationresult(resultlist, nsim, drop))
    }
    
    #################### Prepare for simulations ########################

    ## Set up some data for repeated use in the simulations

    if(Model == "I") {
      ## Compute approximate marginal distribution of type
      if(vector.arg)
        ptypes <- f/sum(f)
      else if(list.arg) {
        fintegrals <- unlist(lapply(flist, integratexy, win=win, ...))
        ptypes <- fintegrals/sum(fintegrals)
      } else {
        ## single argument
        if(is.constant(f) || is.im(f)) {
          ptypes <- rep.int(1/ntypes, ntypes)
        } else {
          ## f is a function (x,y,m)
          ## convert to images and integrate
          fintegrals <- unlist(lapply(types,
                                      integratewithfixedmark,
                                      win=win, fun=f, ...))
          ## normalise
          ptypes <- fintegrals/sum(fintegrals)
        }
      }
    }
    if(Model == "III") {
      ## multinomial: fixed number n[i] of types[i]
      repmarks <- factor(rep.int(types, n), levels=types)
    }
      
    #################### RANDOM GENERATION, general case  ########################

    resultlist <- vector(mode="list", length=nsim)

    for(isim in seq_len(nsim)) {
      ## First generate types, then locations given types
      switch(Model,
             I = ,
             II = {
               ## i.i.d.: n marks with distribution 'ptypes'
               marques <- sample(factortype, n, prob=ptypes, replace=TRUE)
               nn <- table(marques)
             },
             III = {
               ## multinomial: fixed number n[i] of types[i]
               marques <- sample(repmarks)
               nn <- n
             })
      ntot <- sum(nn)

      if(same.density) {
        ## All types have the same conditional density of location;
        ## generate the locations using rpoint
        X <- rpoint(ntot, flist[[1]], maxes[[1]], win=win, ...,
                    giveup=giveup, verbose=verbose)
        resultlist[[isim]] <- X %mark% marques
      } else {
        ## Invoke rpoint() for each type separately
        ## Set up 'blank' pattern
        X <- ppp(numeric(ntot), numeric(ntot), window=win, marks=marques,
                 check=FALSE)
        for(i in seqtypes) {
          if(verbose) cat(paste("Type", i, "\n"))
          if(single.arg && is.function(f)) {
            ## want to call f(x,y,m, ...)
            Y <- rpoint(nn[i], funwithfixedmark, fmax=maxes[i], win=win,
                        ..., m=factortype[i], fun=f, giveup=giveup, verbose=verbose)
          } else {
            ## call f(x,y, ...) or use other formats
            Y <- rpoint(nn[i], flist[[i]], fmax=maxes[i], win=win,
                        ..., giveup=giveup, verbose=verbose)
          }
          Y <- Y %mark% factortype[i]
          X[marques == factortype[i]] <- Y
        }
        resultlist[[isim]] <- X
      }
    }
    return(simulationresult(resultlist, nsim, drop))
  }

  rmpoint
})

rmpoint.I.allim <- local({

  ## Extract pixel coordinates and probabilities
  get.stuff <- function(imag) {
    w <- as.mask(as.owin(imag))
    dx <- w$xstep
    dy <- w$ystep
    rxy <- rasterxy.mask(w, drop=TRUE)
    xpix <- rxy$x
    ypix <- rxy$y
    ppix <- as.vector(imag$v[w$m]) ## not normalised - OK
    npix <- length(xpix)
    return(list(xpix=xpix, ypix=ypix, ppix=ppix,
                dx=rep.int(dx,npix), dy=rep.int(dy, npix),
                npix=npix))
  }

  rmpoint.I.allim <- function(n, f, types, nsim=1) {
    ## Internal use only!
    ## Generates random marked points (Model I *only*)
    ## when all f[[i]] are pixel images.
    ##
    stuff <- lapply(f, get.stuff)
    ## Concatenate into loooong vectors
    xpix <- unlist(lapply(stuff, getElement, name="xpix"))
    ypix <- unlist(lapply(stuff, getElement, name="ypix"))
    ppix <- unlist(lapply(stuff, getElement, name="ppix"))
    dx   <- unlist(lapply(stuff, getElement, name="dx"))
    dy   <- unlist(lapply(stuff, getElement, name="dy"))
    ## replicate types
    numpix <- unlist(lapply(stuff, getElement, name="npix"))
    tpix <- rep.int(seq_along(types), numpix)
    npix <- sum(numpix)
    ##
    resultlist <- vector(mode="list", length=nsim)
    for(isim in seq_len(nsim)) {
      ## sample pixels from union of all images
      id <- sample(npix, n, replace=TRUE, prob=ppix)
      ## get pixel centre coordinates and randomise within pixel
      x <- xpix[id] + (runif(n) - 1/2) * dx[id]
      y <- ypix[id] + (runif(n) - 1/2) * dy[id]
      ## compute types
      marx <- factor(types[tpix[id]],levels=types)
      ## et voila!
      resultlist[[isim]] <- ppp(x, y, window=as.owin(f[[1]]), marks=marx, check=FALSE)
    }
    return(resultlist)
  }

  rmpoint.I.allim
})

##
##     wrapper for Rolf's function
##
rpoint.multi <- function (n, f, fmax=NULL, marks = NULL,
                          win = unit.square(),
                          giveup = 1000, verbose = FALSE,
                          warn=TRUE, nsim=1, drop=TRUE) {
  check.1.integer(nsim)
  stopifnot(nsim >= 0)

  if(is.null(marks) || (is.factor(marks) && length(levels(marks)) == 1)) {
    ## Unmarked
    if(is.function(f)) {
      return(rpoint(n, f, fmax, win, nsim=nsim, warn=warn, giveup=giveup, verbose=verbose))
    } else {
      return(rpoint(n, f, fmax,      nsim=nsim, warn=warn, giveup=giveup, verbose=verbose))
    }
  }

  ## multitype case: validate arguments
  if(length(marks) != n)
    stop("length of marks vector != n")
  if(!is.factor(marks))
    stop("marks should be a factor")
  types <- levels(marks)
  types <- factor(types, levels=types)
  nums <- table(marks)
  
  if(warn) {
    nhuge <- spatstat.options("huge.npoints")
    if(n > nhuge) {
      whinge <- paste("Attempting to generate", n, "random points")
      if(nsim > 1) whinge <- paste(whinge, paren(paste(nsim, "times")))
      warning(whinge, call.=FALSE)
    }
  }

  ## simulate
  resultlist <- vector(mode="list", length=nsim)
  ## generate required number of points of each type
  Xlist <- rmpoint(nums, f, fmax, win=win, types=types,
                   nsim=nsim, drop=FALSE,
                   giveup=giveup, verbose=verbose)
  ## Reorder them to correspond to the desired 'marks' vector
  for(isim in seq_len(nsim)) {
    X <- Xlist[[isim]]
    if(any(table(marks(X)) != nums))
      stop("Internal error: output of rmpoint illegal")
    Y <- X
    Xmarks <- marks(X)
    for(ty in types) {
      to   <- (marks == ty)
      from <- (Xmarks == ty)
      if(sum(to) != sum(from))
        stop(paste("Internal error: mismatch for mark =", ty))
      if(any(to)) {
        Y$x[to] <- X$x[from]
        Y$y[to] <- X$y[from]
        Y$marks[to] <- ty
      }
    }
    resultlist[[isim]] <- Y
  }
  return(simulationresult(resultlist, nsim, drop))
}


  
  

    
