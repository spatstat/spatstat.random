#'   rclusterBKBC.R
#'
#'   $Revision: 1.9 $ $Date: 2025/11/21 01:24:47 $
#' 
#'   Simulation of stationary cluster process
#'   using Brix-Kendall algorithm and Baddeley-Chang modification
#'
#'   Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022-2023
#'   GNU Public Licence >= 2

rclusterBKBC <- function(clusters="Thomas",
               kappa, mu, scale, ..., 
               W=unit.square(), nsim=1, drop=TRUE,
               best=FALSE,
               external=c("BK", "superBK", "border"),
               internal=c("dominating", "naive"),
               inflate=1,
               psmall=1e-4,
               use.inverse=TRUE,
               use.special=TRUE,
               integralmethod=c("quadrature", "trapezoid"),
               verbose=TRUE, warn=TRUE) {
  external.given <- !missing(external)
  integralmethod <- match.arg(integralmethod)
  if(!missing(inflate)) {
    if(is.numeric(inflate)) {
      check.1.real(inflate)
      stopifnot(inflate >= 1)
      if(verbose && inflate > 1) splat("Inflating disc by factor", inflate)
    } else if(is.function(inflate)) {
      if(verbose) splat("Using a user-supplied inflation rule")
    } else if(identical(inflate, "optimal")) {
      if(verbose) splat("Using optimal inflation rule")
    } else stop("Argument 'inflate' should be a number or a function")
  }
  ## ............. Information about the model .............................
  cinfo <- spatstatClusterModelInfo(clusters) # error if unrecognised
  sinfo <- spatstatClusterSimInfo(clusters)      # NULL if unrecognised
  iscompact <- isTRUE(sinfo$iscompact)
  ## Validate parameters and convert to native parametrisation
  par.generic <- c(kappa=kappa, scale=scale)
  par <- cinfo$checkpar(par.generic, old=TRUE)
  if(length(cinfo$shapenames)) { ## there are parameters controlling shape of kernel
    shapestuff <- cinfo$resolveshape(...)$covmodel ## resolve the parameter values
    shapemodel <- shapestuff$model ## optional: name of covariance model e.g. 'exponential'
    shapeargs <- shapestuff$margs ## values of the shape parameters
  } else shapemodel <- shapeargs <- NULL
  ## Assemble model information in format required by cluster simulation functions in bkinfo.R
  mod <- list(par=par, mu=mu, shapemodel=shapemodel, shapeargs=shapeargs)
  
  ## ................ Geometry .............................
  ## shift window to convenient origin
  oldW <- W
  oldcentre <- as.numeric(centroid.owin(Frame(oldW)))
  W <- shift(oldW, -oldcentre)
  ## enclose it in a disc
  rD <- with(vertices(Frame(W)), sqrt(max(x^2+y^2)))
  if(verbose) splat("Disc radius rD =", rD)
  ## D <- disc(radius=rD)
  ## inflated disc
  if(identical(inflate, "optimal")) {
    rE <- optimalinflation(clusters, mod, rD)
  } else if(is.function(inflate)) {
    rE <- inflate(mod, rD)
    if(rE < rD) stop("Function 'inflate' yielded rE < rD")
  } else if(inflate == 1) {
    ## default
    rE <- rD
  } else {
    rE <- inflate * rD
    if(iscompact) {
      rmax <- cinfo$range(par=par.generic, margs=shapeargs)
      rEmax <- rD + rmax
      if(rE > rEmax) {
        if(verbose) splat("Kernel has compact support;",
                          "reducing inflated radius from", rE, "to", rEmax)
        rE <- rEmax
      }
    }
  }
  if(rE == rD) {
    if(verbose) splat("Disc will not be inflated")
    discEname <- "disc"
  } else {
    if(verbose) splat("Inflated radius rE =", rE)
    discEname <- "inflated disc"
  }
  
  ## ............. Decide on computation policy ..........................................
  use.special <- isTRUE(use.special) && !is.null(sinfo) # TRUE if using specialised code in 'sinfo'
  use.inverse <- isTRUE(use.inverse)
  sizecode <- 1L + (scale > rD/10) + (scale > rD/2)
  external <- if(best) switch(sizecode, "border", "BK", "superBK")           else match.arg(external) 
  internal <- if(best) switch(sizecode, "naive", "dominating", "dominating") else match.arg(internal)
  if(use.special && external == "superBK") {
    ## check whether this is supported
    if(!all(c("rhoplusplus", "Mplusplus", "MplusplusInf") %in% names(sinfo))) {
      external <- "BK"
      if(!best && external.given)
        message("Superdominating strategy (external='superBK') is not supported for this cluster process; reverting to 'BK'")
    }
  }
  ## create empty pattern in required format
  emptypattern <- ppp(window=oldW)
  attr(emptypattern, "parents") <- list(x=numeric(0), y=numeric(0))
  attr(emptypattern, "parentid") <- integer(0)
  ## Threshold for warnings about large number of points
  if(warn) 
    nhuge <- spatstat.options("huge.npoints")
  ## ...........................................
  ## Define components of simulation algorithm
  ## ...........................................
  ## intensity of parents conditioned on having any offspring anywhere
  kappadag <- kappa * (1 - exp(-mu))
  ## integral over D of dominating kernel, given distance to origin from parent
  if(use.special) {
    ## function to compute expected number of dominating offspring of a parent at given distance from origin
    Eplus <- sinfo$Eplus
    hdomfun <- sinfo$hdom
    ## random generator of offspring
    roffspring <- sinfo$roffspring
  } else {
    ## generic code
    kernel <- cinfo$kernel # function(par, r, ..., model, margs)
    hdomfun <- function(r, mod, rD) {
      z <- numeric(length(r))
      high <- (r > rD)
      z[!high] <- with(mod, kernel(par, 0, model=shapemodel, margs=shapeargs))
      z[high] <- with(mod,  kernel(par, (r[high]-rD), model=shapemodel, margs=shapeargs))
      return(z)
    }
    Eplus <- function(r, ...) { mu * pi * rD^2 * hdomfun(r, mod, rD) }
    ## random generator of offspring
    cinfo.roff <- cinfo$roffspring
    if(!is.function(cinfo.roff))
      stop("I don't know how to generate clusters of this kind")
    roffspring <- function(n, mod) {
      with(mod, cinfo.roff(n, par, model=shapemodel, margs=shapeargs))
    }
  }
  ## For Brix-Kendall and super-dominating algorithms,
  ## kernel must be bounded at r=0
  LowerLimit <- 0
  if(is.infinite(hdomfun(0, mod, rD))) {
    offence <- remedy <- NULL
    if(internal != "naive") {
      offence <- "Cannot use Brix-Kendall type algorithm"
      remedy  <- "Switching to hybrid algorithm"
      internal <- "naive"
    }
    if(inflate == 1) {
      offence <- c(offence, "Cannot use inflate = 1")
      remedy <- c(remedy, "Setting inflate = 2")
      inflate <- 2
      rE <- inflate * rD
    }
    if(length(offence)) 
      message(paste("The kernel is infinite at distance zero.",
                    paste0(paste(offence, collapse="; "), "."),
                    paste0(paste(remedy, collapse="; "), ".")))
    if(!use.special) 
      LowerLimit <- rE
  }
  ## intensity of Brix-Kendall dominating parent process, given distance to origin
  switch(external,
         border = {
           ## determine border width 
           if(iscompact) {
             ## compact support: determine radius of disc containing support
             rbord <- cinfo$range(par=par.generic, model=shapemodel, margs=shapeargs)
           } else {
             ## non-compact support: determine radius of disc which has probability (1 - psmall)
             rbord <- cinfo$range(par=par.generic, model=shapemodel, margs=shapeargs, thresh=psmall)
           }
           Rmax <- rD + rbord
           if(Rmax > rE) {
             Eborder <- kappadag * pi * (Rmax^2 - rE^2)
             if(verbose) splat("External parents: border method",
                               "\n\tBorder width:", signif(Rmax - rE, 3),
                               "\n\tMaximum distance of parent from origin:", signif(Rmax, 3),
                               "\n\tExpected number of parents in border:", round(Eborder, 3))
           } else {
             Eborder <- 0
             if(verbose)
               splat("External parents: border method\n\tNot required; covered by inflated disc")
           }
         },
         BK = {
           ## .......................
           ## original B-K algorithm 
           ## .......................
           if(verbose) splat("External parents: Brix-Kendall dominating process")
           ## Radial cumulative integral of dominating parent intensity
           if(use.special) {
             rhoplus <- sinfo$rhoplus
             sinfoMplus <- sinfo$Mplus
             Mplus <- function(r) {
               sinfoMplus(r, mod, rD, method=integralmethod)
             }
             ## Integral of dominating parent intensity = Mplus(Inf)
             MplusInf <- sinfo$MplusInf(mod, rD)
           } else {
             ## generic code
             rhoplus <- function(r, ...) {
               kappa <- mod$par[["kappa"]]
               z <- numeric(length(r))
               above <- (r > rD)
               z[!above] <- Eplus(0)
               z[above]  <- Eplus(r[above])
               return(kappa * (1 - exp(-z)))
             }
             MplusIntegrand <- function(r) {
               2 * pi * r * rhoplus(r)
             }
             ## Mplus(r) is the integral from LowerLimit to r
             Mplus <- function(r) {
               rstart <- LowerLimit
               if(rstart == 0) {
                 z <- pi * pmin(r, rD)^2 * rhoplus(0)
                 high <- (r > rD)
                 if(any(high)) 
                   z[high] <- z[high] + indefinteg(MplusIntegrand,
                                                   r[high],
                                                   lower=rD,
                                                   method=integralmethod)
               } else {
                 z <- numeric(length(r))
                 high <- (r > rstart)
                 if(any(high))
                   z[high] <- indefinteg(MplusIntegrand,
                                         r[high],
                                         lower=rstart,
                                         method=integralmethod)
               }
               return(z)
             }
             MplusInf <- Mplus(Inf)
           }
           if(verbose) splat("\tTotal integral of intensity of dominating parents:", round(MplusInf, 2))
           ## We will generate values uniformly distributed between Mmin and MplusInf
           Mmin <- switch(internal,
                          dominating=0,
                          naive=Mplus(rE))
           if(verbose && internal == "naive") {
             splat("\tIntegral inside", paste0(discEname, ":"), round(Mmin, 2))
             splat("\tDifference:", MplusInf - Mmin)
           }
           ## use a slightly lower maximum, to ensure finite radii
           Margin <- MplusInf - Mmin
           eps <- sqrt(.Machine$double.eps)
           if(Margin >= 0) {
             delta <- 1e-4 * Margin
             delta <- max(min(delta, 0.1), eps)
             Mmax <- MplusInf - delta
           } else {
             if(Margin < -eps)
               warning(paste0("Internal numerical problem: MplusInf < Mmin ",
                              paren(paste("difference", Margin)), 
                              "; reset MplusInf = Mmin"))
             Mmax <- MplusInf <- Mmin 
           }
           if(Mmin >= Mmax) {
             Mmax <- MplusInf <- Mmin
             if(verbose) splat("\tNo dominating parents required (Mmax <= Mmin)")
           } else {
             if(verbose) splat("\tExpected number of dominating parents to be generated:", round(Mmax-Mmin, 2))
             ## Inverse function of Mplus available?
             inverseMplus <- sinfo$invMplus
             use.inverse <- use.special && isTRUE(use.inverse) && is.function(inverseMplus)
             if(!use.inverse) {
               ## Numerical root-finding is required
               ## Mplus(r) is proportional to r^2 on [0, rD]
               MplusrD <- Mplus(rD)
               ## Find upper bound on distance to origin corresponding to 'Mmax'
               Rmax <- rD + scale
               doublings <- 0L
               while(Mplus(Rmax) < Mmax) {
                 Rmax <- 2 * Rmax
                 doublings <- doublings + 1L
                 if(doublings > 1e5) stop("Internal error: no solution for Mplus(Rmax) >= Mmax")
               }
               if(verbose) splat("\tUpper bound on distance from parent to centre of window:", signif(Rmax, 4))
               Contrast <- function(x, m) { Mplus(x) - m  }
               SolveRadius <- function(m, rD, rmax, MrD) {
                 if(m <= MrD) return(rD * sqrt(m/MrD))
                 if(m >= Mmax) return(Rmax)
                 uniroot(Contrast, c(rD, 1.1*rmax), m=m)[["root"]]
               }
             }
           }
         },
         superBK = {
           ## ..............................................................................................
           ## Use a process that dominates the dominating process
           ## ..............................................................................................
           if(verbose) splat("External parents: super-dominating process")
           if(use.special) {
             rhoplus <- sinfo$rhoplus
             ## superdominating parent intensity, given distance to origin
             rhoplusplus <- sinfo$rhoplusplus
             ## radial cumulative integral of intensity for superdominating process
             ## i.e. Mplusplus(r) = int_0^r { s * rhoplusplus(s) } ds
             sinfoMplusplus <- sinfo$Mplusplus
             Mplusplus <- function(r) {
               sinfoMplusplus(r, mod, rD, method=integralmethod)
             }
             ## integral of superparent intensity = Mplusplus(Inf)
             MplusplusInf <- sinfo$MplusplusInf(mod, rD)
           } else {
             ## generic code
             rhoplus <- function(r, ...) {
               kappa <- mod$par[["kappa"]]
               z <- numeric(length(r))
               above <- (r > rD)
               z[!above] <- Eplus(0)
               z[above]  <- Eplus(r[above])
               return(kappa * (1 - exp(-z)))
             }
             rhoplusplus <- function(r, ...) {
               kappa <- mod$par[["kappa"]]
               z <- numeric(length(r))
               above <- (r > rD)
               z[!above] <- Eplus(0)
               z[above]  <- Eplus(r[above])
               return(kappa * z)
             }
             MplusplusIntegrand <- function(r) {
               2 * pi * r * rhoplusplus(r)
             }
             Mplusplus <- function(r) {
               rstart <- LowerLimit
               if(rstart == 0) {
                 z <- pi * pmin(r, rD)^2 * rhoplusplus(0)
                 high <- (r > rD)
                 if(any(high)) 
                   z[high] <- z[high] + indefinteg(MplusplusIntegrand,
                                                   r[high],
                                                   lower=rD,
                                                   method=integralmethod)
               } else {
                 z <- numeric(length(r))
                 high <- (r > rstart)
                 if(any(high)) 
                   z[high] <- indefinteg(MplusplusIntegrand,
                                         r[high],
                                         lower=rstart,
                                         method=integralmethod)
               }
               return(z)
             }
             MplusplusInf <- Mplusplus(Inf)
           }
           if(verbose) splat("\tTotal integral of intensity of superparents:", round(MplusplusInf, 2))
           ## generate values uniformly distributed between Mmin and MplusplusInf
           Mmin <- switch(internal,
                          dominating=0,
                          naive=Mplusplus(rE))
           if(verbose && internal == "naive") {
             splat("\tIntegral inside", paste0(discEname, ":"), round(Mmin, 2))
             splat("\tDifference:", MplusplusInf - Mmin)
           }
           ## use a slightly lower maximum, to ensure finite radii
           Margin <- MplusplusInf - Mmin
           eps <- sqrt(.Machine$double.eps)
           if(Margin >= 0) {
             delta <- 1e-4 * Margin
             delta <- max(min(delta, 0.1), eps)
             Mmax <- MplusplusInf - delta
           } else {
             if(Margin < -eps)
               warning(paste0("Internal numerical problem: ",
                              "MplusplusInf < Mmin ",
                              paren(paste("difference", Margin)),
                              "; reset MplusplusInf = Mmin"))
             Mmax <- MplusplusInf <- Mmin 
           }
           if(Mmin >= Mmax) {
             Mmax <- MplusplusInf <- Mmin 
             if(verbose) splat("\tNo superparents required (Mmax = Mmin)")
           } else {
             if(verbose) splat("\tExpected number of superparents to be generated:", round(Mmax-Mmin, 2))
             ## Inverse function of Mplusplus available?
             inverseMplusplus <- sinfo$invMplusplus
             use.inverse <- use.special && isTRUE(use.inverse) && is.function(inverseMplusplus)
             if(!use.inverse) {
               ## Numerical root-finding required
               ## Mplusplus(r) is proportional to r^2 on [0, rD]
               MplusplusrD <- Mplusplus(rD)
               ## Find upper bound on solution corresponding to 'Mmax'
               Rmax <- rD + scale
               doublings <- 0L
               while(Mplusplus(Rmax) < Mmax) {
                 Rmax <- 2 * Rmax
                 doublings <- doublings + 1L
                 if(doublings > 1e5) stop("Internal error: no solution for Mplusplus(Rmax) >= Mmax")
               }
               if(verbose)
                 splat("Upper bound on distance from superparent to centre of window:", signif(Rmax, 6))
               Contrast <- function(x, m) { Mplusplus(x) - m }
               SolveRadius <- function(m, rD, rmax, MrD) {
                 if(m <= MrD) return(rD * sqrt(m/MrD))
                 if(m >= Mmax) return(Rmax)
                 uniroot(Contrast, c(rD, 1.1*rmax), m=m)[["root"]]
               }
             }
           }
         })
  
  ## >>>>>>>>>>>>>>>  s t a r t    s i m u l a t i o n    l o o p   <<<<<<<<<<<<<<<<<<<<<<<
  resultlist <- vector(mode="list", length=nsim)
  for(isim in seq_len(nsim)) {
    if(verbose && nsim > 1) splat("Generating realisation", isim)
    cost <- 0
    ## .......................
    ##   PARENTS OUTSIDE DISC
    ## .......................
    switch(external,
           border = {
             ## .................................................................
             ##  naive: generate parents uniformly in border region
             ## .................................................................
             externalparentname <- "dominating parents in border region"
             if(Eborder == 0) {
               ndom <- 0
             } else {
               ndom <- rpois(1, Eborder)
               cost <- cost + ndom
               if(verbose) splat("Generated", externalparentname, "\n\tNumber:", ndom)
               if(warn && ndom > nhuge) {
                 whinge <- paste("Generating", ndom, externalparentname)
                 message(whinge)
                 warning(whinge, call.=FALSE)
               }
             } 
             if(ndom > 0)
               rdom <- sqrt(runif(ndom, min=rE^2, max=Rmax^2))
           },
           BK = {
             ## .................................................................
             ## original Brix-Kendall (possibly restricted to parents outside E)
             ## .................................................................      
             ## generate Poisson number of parents in dominating process
             ##                  (possibly restricted to locations outside E)
             externalparentname <- if(internal != "naive") "dominating parents" else
                                   paste("dominating parents outside", discEname)
             if(MplusInf > Mmin) {
               ndom <- rpois(1, MplusInf - Mmin)
               cost <- cost + ndom
               if(verbose) splat("Generated", externalparentname, "\n\tNumber:", ndom)
               if(warn && ndom > nhuge) {
                 whinge <- paste("Generating", ndom, externalparentname)
                 message(whinge)
                 warning(whinge, call.=FALSE)
               }
             } else {
               ndom <- 0
             }
             ## generate parent locations
             if(ndom > 0) {
               mdom <- runif(ndom, min=Mmin, max=Mmax)
               if(use.inverse) {
                 ## inverse function is analytic
                 rdom <- inverseMplus(mdom, mod, rD, MplusInf)
               } else {
                 ## use root-finding
                 mdom <- sort(mdom, decreasing=TRUE)
                 rdom <- numeric(ndom)
                 rmax <- Rmax
                 ## Mplus(r) is proportional to r^2 on [0, rD]
                 below <- (mdom <= MplusrD)
                 rdom[below] <- rD * sqrt(mdom[below]/MplusrD)
                 for(i in which(!below)) {
                   ## solve for radius, and use this as the upper bound for the next problem
                   rmax <- rdom[i] <- SolveRadius(mdom[i], rD, rmax, MplusrD)
                 }
               }
               if(verbose) splat("\tRange of distances from origin:", prange(signif(range(rdom), 3)))
             }
           },
           superBK = {
             ## .................................................................
             ##    Brix-Kendall with super-dominating process
             ##                     (possibly restricted to parents outside E)
             ## .................................................................
             externalparentname <- if(internal != "naive") "super-parents" else
                                   paste("super-parents outside", discEname)
             ## generate Poisson number of superparents
             if(MplusplusInf > Mmin) {
               nsuper <- rpois(1, MplusplusInf - Mmin)
               cost <- cost + nsuper
               if(verbose) splat("Generated", externalparentname, "\n\tNumber:", nsuper)
               if(warn && nsuper > nhuge) {
                 whinge <- paste("Generating", nsuper, externalparentname)
                 message(whinge)
                 warning(whinge, call.=FALSE)
               }
             } else {
               nsuper <- 0
             }
             ## generate superparent locations
             if(nsuper == 0) {
               ndom <- 0
             } else {
               ## generate radii by inverting Mplusplus
               msuper <- runif(nsuper, min=Mmin, max=Mmax)
               rsuper <- numeric(nsuper)
               if(use.inverse) {
                 ## inverse function is analytic
                 rsuper <- inverseMplusplus(msuper, mod, rD, MplusplusInf)
               } else {
                 ## use numerical root-finding
                 msuper <- sort(msuper, decreasing=TRUE)
                 rmax <- Rmax
                 ## Mplusplus(r) is proportional to r^2 on [0, rD]
                 below <- (msuper < MplusplusrD)
                 rsuper[below] <- rD * sqrt(msuper[below]/MplusplusrD)
                 for(i in which(!below))
                   rmax <- rsuper[i] <- SolveRadius(msuper[i], rD, rmax, MplusplusrD)
               }
               if(verbose) splat("\tRange of distances from origin:", prange(signif(range(rsuper), 3)))
               ## thin the superdominating parents to dominating parents
               rdom <- rsuper[rhoplus(rsuper, mod, rD) > runif(nsuper) * rhoplusplus(rsuper, mod, rD)]
               ndom <- length(rdom)
               if(verbose) {
                 splat("Thinned superparents to obtain dominating parents\n",
                       "\tNumber of dominating parents:", ndom,
                       "\n\t", "Acceptance rate:", signif(as.numeric(ndom)/nsuper, 4))
                 if(ndom > 0)
                   splat("\tRange of distances from origin:", prange(signif(range(rdom), 3)))
               }
             }
           })

    ## Now generate offspring 
    if(ndom == 0) {
      noff <- 0L
    } else {
      ## generate locations of [dominating/proposed] parents
      theta <- runif(ndom, max=2*pi)
      xp <- rdom * cos(theta)
      yp <- rdom * sin(theta)
      if(verbose) {
        splat("Generated spatial locations of", ndom, externalparentname)
        splat("\n\nGenerating offspring of", externalparentname)
      }
      ## generate offspring 
      switch(external,
             border = {
               ##  for each proposed parent, generate nonzero number of offspring,
               ##                            according to cluster mechanism
               offspringname <- paste("offspring of", externalparentname)
               noffeach <- rpoisnonzero(ndom, mu)
               noff <- sum(noffeach)
               cost <- cost + noff
               if(verbose) splat("Generated", offspringname,
                                 "\n\tNumber of offspring:", noff)
               if(warn && noff > nhuge) {
                 whinge <- paste("Generating", noff, offspringname)
                 message(whinge)
                 warning(whinge, call.=FALSE)
               }
               ## assign offspring to parents
               parentid <- rep.int(1:ndom, noffeach)
               ## offspring are random displacements of parents
               displace <- roffspring(noff, mod)
               xoff <- xp[parentid] + displace$x
               yoff <- yp[parentid] + displace$y
               ## clip to window
               retain <- inside.owin(xoff, yoff, W)
               noff <- sum(retain)
               if(verbose) splat("Clipped to window\n\tNumber of offspring in window:", noff)
             },
             BK = ,
             superBK = {
               ## generate nonzero number of offspring in disc, for each dominating parent
               Edom <- Eplus(rdom, mod, rD)
               noffeach <- rpoisnonzero(ndom, Edom)
               noff <- sum(noffeach)
               cost <- cost + noff
               offspringname <- paste("offspring (inside disc) of",
                                      externalparentname)
               if(verbose) splat("Generated", offspringname,
                                 "\n\tNumber of offspring:", noff)
               if(warn && noff > nhuge) {
                 whinge <- paste("Generating", noff, offspringname)
                 message(whinge)
                 warning(whinge, call.=FALSE)
               }
               ## offspring are uniform in disc
               roff <- sqrt(runif(noff, max=rD^2))
               toff <- runif(noff, max=2*pi)
               xoff <- roff * cos(toff)
               yoff <- roff * sin(toff)
               ## assign offspring to parents
               parentid <- rep.int(1:ndom, noffeach)
               ## thin using correct kernel
               dx <- xoff - xp[parentid]
               dy <- yoff - yp[parentid]
               dr <- sqrt(dx^2 + dy^2)
               hcorrect <- cinfo$kernel(par=par, rvals=dr, model=shapemodel, margs=shapeargs)
               hdom     <- hdomfun(rdom, mod, rD)[parentid]
               if(any(bad <- (hcorrect > hdom)))
                 stop(paste("Internal error:", sum(bad),
                   "values of the dominating kernel do not dominate the true kernel"),
                   call.=FALSE)
                 retain <- (hcorrect >= runif(noff) * hdom)
                 if(verbose) splat("Thinned offspring to correct process\n",
                                   "\tNumber of retained offspring in disc:", sum(retain),
                                   "\n\tAcceptance rate:", signif(mean(retain), 4))
                 if(any(retain)) {
                   ## finally clip to window
                   retain[retain] <- inside.owin(xoff[retain], yoff[retain], W)
                 }
                 noff <- sum(retain)
                 if(verbose) splat("Clipped to window\n\tNumber of offspring in window:", noff)
             })
    }
    ## ...................................
    ##    Apply the thinning
    ## ...................................
    if(noff > 0) {
      xoff <- xoff[retain]
      yoff <- yoff[retain]
      parentid <- parentid[retain]
      retainedparents <- sort(unique(parentid))
      parentid <- match(parentid, retainedparents)
      xp <- xp[retainedparents]
      yp <- yp[retainedparents]
      np <- length(retainedparents)
    } else {
      xoff <- yoff <- xp <- yp <- numeric(0)
      parentid <- integer(0)
      np <- 0L
    }
    if(verbose)
      splat("Total number of offspring points",
            if(internal == "naive") "(of parents outside disc):" else ":",
            noff)

    ## ................................
    ##   PARENTS INSIDE (inflated) DISC
    ## ................................

    switch(internal,
           dominating = {
             ## Parents were already generated above
             if(verbose) splat("Internal parents:",
                               switch(external,
                                      BK="Brix-Kendall dominating process",
                                      superBK = "super-dominating process",
                                      border = "border method"),
                               "(already generated)")
           },
           naive = {
             ## Generate parents inside (inflated) disc 
             ##  conditional on having at least 1 offspring **anywhere**
             internalparentname <- paste("additional parents inside", discEname)
             if(verbose) {
               splat("Internal parents: naive method",
                     "\n\tGenerating", internalparentname, 
                     "with intensity", signif(kappadag, 4))
             }
             npin <- rpois(1, kappadag * pi * rE^2)
             rpin <- sqrt(runif(npin, max=rE^2))
             tpin <- runif(npin, max=2*pi)
             xpin <- rpin * cos(tpin)
             ypin <- rpin * sin(tpin)
             cost <- cost + npin
             if(verbose) splat("\tNumber of",
                               paste0(internalparentname, ":"), npin)
             if(warn && npin > nhuge) {
               whinge <- paste("Generating", npin, internalparentname)
               message(whinge)
               warning(whinge, call.=FALSE)
             }
             if(npin > 0) {
               offspringname <- paste("offspring of", internalparentname)
               if(verbose) splat("\nGenerating", offspringname)
               noff.in.each <- rpoisnonzero(npin, mu)
               noff.in <- sum(noff.in.each)
               cost <- cost + noff.in
               if(verbose) splat("\tNumber of",
                                 paste0(offspringname, ":"), noff.in)
               if(warn && noff.in > nhuge) {
                 whinge <- paste("Generating", noff.in, offspringname)
                 message(whinge)
                 warning(whinge, call.=FALSE)
               }
               pid.in <- rep.int(1:npin, noff.in.each)
               displace <- roffspring(noff.in, mod)
               xoff.in <- xpin[pid.in] + displace[["x"]]
               yoff.in <- ypin[pid.in] + displace[["y"]]
               ## restrict to window
               retain.in <- inside.owin(xoff.in, yoff.in, W)
               if(verbose) splat("Clipping", offspringname, "to window",
                                 "\n\tNumber retained:", sum(retain.in),
                                 paren(paste("acceptance rate", signif(mean(retain.in), 4))))
               if(any(retain.in)) {
                 ## retain these offspring
                 xoff.in <- xoff.in[retain.in]
                 yoff.in <- yoff.in[retain.in]
                 pid.in  <- pid.in[retain.in]
                 retainedpid.in <- sort(unique(pid.in))
                 pid.in <- match(pid.in, retainedpid.in)
                 xpin <- xpin[retainedpid.in]
                 ypin <- ypin[retainedpid.in]
                 npin <- length(retainedpid.in)
                 ## append to existing
                 xoff <- c(xoff, xoff.in)
                 yoff <- c(yoff, yoff.in)
                 parentid <- c(parentid, np + pid.in)
                 noff <- noff + sum(retain.in)
                 xp <- c(xp, xpin)
                 yp <- c(yp, ypin)
                 np   <- np   + npin
               }
             }
           })

    ## Pack up simulation result [[isim]]
    if(verbose) splat("Final total number of offspring points:", noff)
    if(noff == 0) {
      Y <- emptypattern
    } else {
      Y <- ppp(xoff + oldcentre[1],
               yoff + oldcentre[2],
               window=oldW)
      attr(Y, "parents") <- list(x = xp + oldcentre[1],
                                 y = yp + oldcentre[2])
      attr(Y, "parentid") <- parentid
    }
    attr(Y, "cost") <- cost
    resultlist[[isim]] <- Y
  }
  ## >>>>>>>>>>>>>>>   e n d       l o o p   <<<<<<<<<<<<<<<<<<<<<<<
  result <- simulationresult(resultlist, nsim, drop)
  return(result)
}
  
optimalinflation <- function(clusters, mod, rD) {
  ## Optimal inflated disc radius rE, determined by root-finding
  ## as in section 6.6 of Baddeley and Chang (2023)
  cinfo <- spatstatClusterModelInfo(clusters) # error if unrecognised
  par.native <- cinfo$checkpar(mod$par, native=TRUE)
  mu    <- mod$mu
  margs <- mod$shapeargs
  a <- if(mu == 0) 1 else (1 + (1-exp(-mu))/mu)
  b <- a/(pi * rD^2)
  ## solve the equation h(r) = b where h is the 2D kernel
  if(is.function(kerinv <- cinfo$kerinverse)) {
    ## use specialised code to solve h(r) = b
    delta <- kerinv(par.native, level=b, margs=margs)
    rE <- rD + delta
  } else if(is.function(kernel <- cinfo$kernel)) {
    ## apply root-finding to solve h(r) = b
    h0 <- kernel(par.native, 0, margs=margs)
    if(h0 <= b) {
      rE <- rD
    } else {
      f <- function(x) { b - kernel(par.native, x, margs=margs) }
      u <- try(uniroot(f, c(0, 100 * scale)), silent=TRUE)
      if(inherits(u, "try-error")) {
        rE <- 2 * rD
      } else {
        rE <- rD + u$root
      }
    }
  } else {
    rE <- 2 * rD
  }
  return(rE)
}
