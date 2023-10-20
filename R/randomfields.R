#'
#'   R/randomfields.R
#'
#'   Random generators of Gaussian random fields
#'
#'   $Revision: 1.7 $ $Date: 2023/10/20 14:50:06 $
#'
#'   Copyright (c) 2023 Adrian Baddeley
#'   GNU Public Licence (>= 2.0)

rGRFgauss <- function(W=owin(), mu=0, var=1, scale, ..., nsim=1, drop=TRUE) {
  ## Gaussian random field with Gaussian covariance function
  check.1.real(scale)
  stopifnot(scale > 0)
  gfun <- function(x) { exp(-(x/scale)^2) }
  rGRFcircembed(W=W, mu=mu, var=var, corrfun = gfun,
             ...,
             nsim=nsim, drop=drop)
}

rGRFexpo <- function(W=owin(), mu=0, var=1, scale, ...,
                     nsim=1, drop=TRUE) {
  ## Gaussian random field with exponential covariance function
  check.1.real(scale)
  stopifnot(scale > 0)
  efun <- function(x) { exp(-x/scale) }
  rGRFcircembed(W=W, mu=mu, var=var, corrfun = efun,
             ...,
             nsim=nsim, drop=drop)
}

rGRFstable <- function(W=owin(), mu=0, var=1, scale,
                       alpha,
                       ...,
                       nsim=1, drop=TRUE) {
  ## Gaussian random field with stable covariance
  check.1.real(alpha)
  sfun <- function(x) { exp(-(x/scale)^alpha) }
  rGRFcircembed(W=W, mu=mu, var=var, corrfun = sfun,
             ...,
             nsim=nsim, drop=drop)
}

rGRFgencauchy <- function(W=owin(), mu=0, var=1, scale,
                          alpha, beta,
                          ...,
                       nsim=1, drop=TRUE) {
  ## Gaussian random field with generalised Cauchy covariance
  check.1.real(alpha)
  check.1.real(beta)
  cfun <- function(x) { (1 + (x/scale)^alpha)^(-beta/alpha) }
  rGRFcircembed(W=W, mu=mu, var=var, corrfun = cfun,
             ...,
             nsim=nsim, drop=drop)
}

rGRFmatern <- function(W=owin(), mu=0, var=1, scale,
                       nu,
                       ...,
                       nsim=1, drop=TRUE) {
  ## Gaussian random field with Matern covariance
  check.1.real(nu)
  mfun <- function(x) {
    z <- (x/scale) * sqrt(2 * nu)
    ifelse(x == 0,
           1, 
           (z^nu) * besselK(z, nu) * (2^(1-nu))/gamma(nu))
  }
  rGRFcircembed(W=W, mu=mu, var=var, corrfun = mfun,
             ...,
             nsim=nsim, drop=drop)
}


## test functions - not for distribution!

## pvar <- function(Zlist, X=ppp(0.5, 0.5, 0:1, 0:1)) {
##   zvals <- sapply(Zlist, "[", i=X)
##   var(zvals)
## }

## pcor <- function(Zlist,
##                  lag=0.1, X=ppp(0.5, 0.5, 0:1, 0:1),
##                  Y=ppp(0.5, 0.5+lag, 0:1, 0:1)) {
##   zXvals <- sapply(Zlist, "[", i=X)
##   zYvals <- sapply(Zlist, "[", i=Y)
##   cov(zXvals, zYvals)/sqrt(var(zXvals) * var(zYvals))
## }
