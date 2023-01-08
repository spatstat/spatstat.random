#'
#'    rvargamma.R
#'
#'   $Revision: 1.6 $ $Date: 2023/01/08 02:26:28 $
#'
#'   Simulation of Variance-Gamma cluster process
#'   using either naive algorithm or BKBC algorithm
#'   (R code)
#'
#'   rVarGamma
#'
#'    Original code for naive simulation of Neyman-Scott by Adrian Baddeley
#'    Original code for simulation of rVarGamma offspring by Abdollah Jalilian
#'    Bug fixes by Abdollah, Adrian Baddeley, and Rolf Turner
#'
#'   Implementation of BKBC algorithm by Adrian Baddeley and Ya-Mei Chang
#' 
#'   Copyright (c) 2000-2023 Adrian Baddeley, Abdollah Jalilian, Ya-Mei Chang 
#'   GNU Public Licence >= 2

##    
## =================================================================
## Neyman-Scott process with Variance Gamma (Bessel) kernel function
## =================================================================

## nu.ker: smoothness parameter of Variance Gamma kernel function
## omega: scale parameter of kernel function
## nu.pcf: smoothness parameter of Variance Gamma pair correlation function
## eta: scale parameter of Variance Gamma pair correlation function
## nu.pcf = 2 * nu.ker + 1    and    eta = omega

rVarGamma <- local({
  
  ## simulates mixture of isotropic Normal points in 2D with gamma variances
  rnmix.gamma <- function(n = 1, shape, rate) {
    V <- matrix(rnorm(2 * n, 0, 1), nrow = n, ncol = 2)
    s <- rgamma(n, shape=shape, rate=rate)
    return(sqrt(s) * V)
  }

  ## main function
  rVarGamma <- function(kappa, scale, mu,
                           nu,
                           win = square(1),
                           nsim=1, drop=TRUE, 
                           ...,
                          algorithm=c("BKBC", "naive"),
                          nonempty=TRUE, 
                          thresh = 0.001,
                          poisthresh=1e-6,
                          expand = NULL,
                          saveparents=FALSE, saveLambda=FALSE,
                          kappamax=NULL, mumax=NULL) {
    ## Variance-Gamma cluster process
    
    ## nu / nu.ker: smoothness parameter of Variance Gamma kernel function
    ## scale / omega: scale parameter of kernel function

    check.1.integer(nsim)
    stopifnot(nsim >= 0)
    if(nsim == 0) return(simulationresult(list()))

    ## Catch old nu.ker/nu.pcf syntax and resolve nu-value.
    dots <- list(...)
    if(missing(nu)){
      nu <- resolve.vargamma.shape(nu.ker=dots$nu.ker,
                                   nu.pcf=dots$nu.pcf,
                                   allow.default=TRUE)$nu.ker
    } else {
      check.1.real(nu)
      stopifnot(nu > -1)
    }
    nu.ker <- nu
    ## Catch old scale syntax (omega)
    if(missing(scale)) scale <- dots$omega
    
    ## Catch old name 'eps' for 'thresh':
    if(missthresh <- missing(thresh))
      thresh <- dots$eps %orifnull% 0.001

    ## determine the effective maximum radius of clusters
    ## (for the naive algorithm, or when kappa is not constant)    
    if(missing(expand)){
      expand <- clusterradius("VarGamma", scale = scale, nu = nu.ker,
                              thresh = thresh, ...)
    } else if(!missthresh) {
      warning("Argument ", sQuote("thresh"), " is ignored when ",
              sQuote("expand"), " is given")
    }

    #' validate 'kappa' and 'mu'
    km <- validate.kappa.mu(kappa, mu, kappamax, mumax,
                            win, expand, ...,
                            context="In rCauchy")
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
      ## (generic R implementation using information from cluster model table)
      result <- do.call(rclusterBKBC,
                        resolve.defaults(
                          list(clusters = "VarGamma",
                               kappa    = kappamax,
                               scale    = scale,
                               mu       = mumax,
                               nu.ker   = nu.ker,
                               W        = quote(win),
                               nsim     = nsim,
                               drop     = FALSE),
                          list(...),
                          list(internal = "naive",
                               external = "super",
                               inflate  = "optimal",
                               verbose=FALSE)))
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
      if(is.numeric(kappa) && 1/(4 * pi * kappamax * scale^2) < poisthresh) {
        if(is.function(mu)) mu <- as.im(mu, W=win, ...)
        kapmu <- kappa * mu
        result <- rpoispp(kapmu, win=win, nsim=nsim, drop=drop, warnwin=FALSE)
        result <- fakeNeyScot(result, kapmu, win, saveLambda, saveparents)
        return(result)
      }
    
      ## gamma mixture of normals
      alpha <- 2 * (nu.ker + 1)
      beta  <- 1/(2 * scale^2)
    
      ## simulate
      result <- rNeymanScott(kappa=kappa,
                             expand=expand,
                             rcluster=list(mu, rnmix.gamma),
                             win=win,
                             shape = alpha/2, # formal argument of rnmix.gamma
                             rate = beta, # formal argument of rnmix.gamma
                             nsim=nsim, drop=FALSE,
                             nonempty = nonempty,
                             saveparents = do.parents,
                             kappamax=kappamax, mumax=mumax)

    }
    if(saveLambda){
      BW <- Frame(win)
      for(i in 1:nsim) {
        parents <- attr(result[[i]], "parents")
        BX <- boundingbox(BW, bounding.box.xy(parents))
        parents <- as.ppp(parents, W=BX, check=FALSE)
        Lambda <- clusterfield("VarGamma", parents, scale=scale,
                               nu=nu.ker, mu=mu, ...)
        attr(result[[i]], "Lambda") <- Lambda[win, drop=FALSE]
      }
    }
    return(simulationresult(result, nsim, drop))
  }

  inflateVarGamma <- function(mod, rD) {
    optimalinflation("VarGamma", mod, rD) 
    }
  rVarGamma
})
