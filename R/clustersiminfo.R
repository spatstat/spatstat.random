#'
#'     clustersiminfo.R
#'
#'   Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022-2023
#'   GNU Public Licence >= 2
#'
#' -----------------------------------------------------------------
#' 
#'  Table of additional information for cluster models
#'  for use in Brix-Kendall and Baddeley-Chang simulation algorithms
#'
#'  The simulation algorithm rclusterBKBC uses this information
#'  provided there is an entry in this table for the desired model,
#'  and provided use.special=TRUE. Otherwise it uses a generic algorithm.
#'
#' -------------------------------------------------------------------
#'
#'  D = bounding disc
#'  rD = radius of D (numeric > 0)
#'  r = distance from parent point to centre of disc (numeric vector)
#'  mod = list containing all model parameters:
#'        $par = original (native) parameters of cluster model e.g. (kappa, sigma2) 
#'        $mu = mean number of offspring per parent in stationary process
#'        $margs = list of shape parameters
#' 
#'  Entries:
#'
#' (A) cluster model
#'
#'     roffspring = function(n, mod)
#'          generate x, y coordinates of n offspring of parent at origin
#'
#' (B) Brix-Kendall dominating process
#'
#'     hdom = function(r, mod, rD)
#'          Value inside D of the dominating offspring density, for a parent at distance r from centre of D
#' 
#'     Eplus = function(r, mod, rD)
#'          Integral over D of dominating kernel, for a parent at distance r from centre of D
#'          Eplus(r) = pi * rD^2 * mu * hdom(r)
#' 
#'     rhoplus = function(r, mod, rD)
#'          Intensity of Brix-Kendall dominating parent process, given distance to origin
#'          rhoplus(r) = kappa * (1 - exp(-Eplus(r)))
#' 
#'     Mplus = function(r, mod, rD, ...)
#'          Radial cumulative integral of intensity of Brix-Kendall dominating intensity
#'          Mplus(r) = \int_0^r 2 pi t rhoplus(t) dt
#' 
#'     MplusInf = function(mod, rD)
#'          Total integral Mplus(infty)
#'          Expected total number of parents  in the dominating process
#'
#'     invMplus = function(v, mod, rD, Minfty)
#'          Inverse function of Mplus (IF KNOWN)
#'
#' (C) Baddeley-Chang super-dominating process
#' 
#'     rhoplusplus = function(r, mod, rD)
#'          Intensity of superdominating process
#'          rhoplusplus(r) = kappa * Eplus(r)
#' 
#'     Mplusplus = function(r, mod, rD, ...)
#'          Radial cumulative integral of superdominating intensity
#'          Mplusplus(r) = \int_0^r 2 pi t rhoplusplus(t) dt
#'                       = kappa \int_0^r 2 pi t Eplus(t) dt
#' 
#'     MplusplusInf = function(mod, rD)
#'          Total integral Mplusplus(infty)
#'          Expected total number of parents in the superdominating process
#' 
#'     invMplusplus = function(v, mod, rD, Minfty)
#'          Inverse function of Mplusplus (IF KNOWN)
#'
#'     inflate = function(mod, rD)
#'          Rule for determining optimal inflated radius rE
#'          according to Baddeley and Chang (2023) section 6.6


#' >>>>>>>>>>>>>>>>   THOMAS PROCESS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

.Thomas.roffspring <- function(n, mod) {
  sigma <- sqrt(mod$par[["sigma2"]])
  list(x=rnorm(n, sd=sigma), y=rnorm(n, sd=sigma))
}

.Thomas.hdom <- function(r, mod, rD) {
  inv2sig2 <- 1/(2 * mod$par[["sigma2"]])
  B <- inv2sig2/pi
  z <- numeric(length(r))
  above <- (r > rD)
  z[!above] <- B
  z[above]  <- B * exp(- inv2sig2 * (r[above] - rD)^2)
  return(z)
}
  
.Thomas.Eplus <- function(r, mod, rD) {
  mu <- mod$mu
  inv2sig2 <- 1/(2 * mod$par[["sigma2"]])
  A <- mu * rD^2 * inv2sig2
  z <- numeric(length(r))
  above <- (r > rD)
  z[!above] <- A
  z[above]  <- A * exp(- inv2sig2 * (r[above]-rD)^2)
  return(z)
}

.Thomas.rhoplus <- function(r, mod, rD) {
  kappa <- mod$par[["kappa"]]
  z <- numeric(length(r))
  above <- (r > rD)
  z[!above] <- kappa * (1 - exp(-.Thomas.Eplus(0,        mod, rD)))
  z[above]  <- kappa * (1 - exp(-.Thomas.Eplus(r[above], mod, rD)))
  return(z)
}

.Thomas.MplusIntegrand <- function(r, mod, rD) {
  2 * pi * r * .Thomas.rhoplus(r, mod, rD)
}

.Thomas.Mplus <- function(r, mod, rD, ..., method="q") {
  ## integrand is linear on [0, rD]
  z <- pi * pmin(r, rD)^2 * .Thomas.rhoplus(0, mod, rD)
  high <- (r > rD)
  if(any(high)) 
    z[high] <- z[high] + indefinteg(.Thomas.MplusIntegrand,
                                    r[high],
                                    lower=rD,
                                    mod=mod, rD=rD, method=method)
  return(z)
}

.Thomas.MplusInf <- function(mod, rD) { .Thomas.Mplus(Inf, mod, rD) }

## inverse function requires numerical root-finding

.Thomas.rhoplusplus <- function(r, mod, rD) {
  mod$par[["kappa"]] * .Thomas.Eplus(r, mod, rD)
}

.Thomas.Mplusplus <- function(r, mod, rD, ...) {
  mu <- mod$mu
  kappa <- mod$par[["kappa"]]
  sigma2 <- mod$par[["sigma2"]]
  inv2sig2 <- 1/(2 * sigma2)
  z <- pi * pmin(r, rD)^2
  high <- (r > rD)
  z[high] <- z[high] +
    2 * pi * sigma2 * (1 - exp(-inv2sig2 * (r[high]-rD)^2)) +
    2 * pi * rD * sqrt(2 * pi * sigma2) * (
      pnorm(r[high]-rD, sd=sqrt(sigma2)) - 1/2
    )
  z <- kappa * mu * rD^2 * inv2sig2 * z
  return(z)
}

.Thomas.MplusplusInf <- function(mod, rD) {
  mu <- mod$mu
  kappa <- mod$par[["kappa"]]
  sigma2 <- mod$par[["sigma2"]]
  inv2sig2 <- 1/(2 * sigma2)
  kappa * mu * pi * rD^2 * (
    rD^2 * inv2sig2 +
    1 +
    rD * sqrt(2 * pi * sigma2) * inv2sig2
  )
}

## inverse function requires numerical root-finding

.Thomas.inflate <- function(mod, rD) {
  mu <- mod$mu
  sigma2 <- mod$par[["sigma2"]]
  a <- if(mu == 0) 1 else (1 + (1-exp(-mu))/mu)
  b <- (rD^2)/(2*a*sigma2)
  if(b <= 1) return(rD)
  delta <- 2 * sqrt(sigma2) * sqrt(log(b)/2)
  return(rD + delta)
}
  
#' >>>>>>>>>>>>>>>>   MATERN CLUSTER PROCESS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

.MatClust.roffspring <- function(n, mod) {
  R <- mod$par[["R"]]
  rad <- R * sqrt(runif(n))
  theta <- runif(n, max=2*pi)
  list(x = rad * cos(theta), y = rad * sin(theta))
}

.MatClust.hdom <- function(r, mod, rD) {
  R <- mod$par[["R"]]
  (r - rD < R)/(pi * R^2)
}
  
.MatClust.Eplus <- function(r, mod, rD) {
  R <- mod$par[["R"]]
  mu <- mod$mu
  mu * (r - rD < R) * (rD/R)^2
}

.MatClust.rhoplus <- function(r, mod, rD) {
  mod$par[["kappa"]] * (1 - exp(-.MatClust.Eplus(r, mod, rD)))
}

.MatClust.Mplus <- function(r, mod, rD, ...) {
  R <- mod$par[["R"]]
  kappa <- mod$par[["kappa"]]
  mu <- mod$mu
  pi * pmin(r, R+rD)^2 * kappa * ( 1 - exp(-mu * (rD/R)^2))
}

.MatClust.MplusInf <- function(mod, rD) {
  R <- mod$par[["R"]]
  kappa <- mod$par[["kappa"]]
  mu <- mod$mu
  pi * (R+rD)^2 * kappa * ( 1 - exp(-mu * (rD/R)^2))
}

.MatClust.inverseMplus <- function(v, mod, rD, Minfty) {
  R <- mod$par[["R"]]
  kappa <- mod$par[["kappa"]]
  mu <- mod$mu
  B <- pi * kappa * (1 - exp(- mu * rD^2/R^2))
  sqrt(v/B)
}

.MatClust.rhoplusplus <- function(r, mod, rD) {
  R <- mod$par[["R"]]
  kappa <- mod$par[["kappa"]]
  mu <- mod$mu
  kappa * mu * (r < (R+rD)) * (rD/R)^2
}

.MatClust.Mplusplus <- function(r, mod, rD, ...) {
  R <- mod$par[["R"]]
  kappa <- mod$par[["kappa"]]
  mu <- mod$mu
  kappa * mu * pi * (pmin(r, R+rD) * rD/R)^2
}

.MatClust.MplusplusInf <- function(mod, rD) {
  R <- mod$par[["R"]]
  kappa <- mod$par[["kappa"]]
  mu <- mod$mu
  kappa * mu * pi * ((R+rD) * rD/R)^2
}

.MatClust.inverseMplusplus <- function(v, mod, rD, Minfty) {
  R <- mod$par[["R"]]
  kappa <- mod$par[["kappa"]]
  mu <- mod$mu
  (R/rD) * sqrt(v/(pi * kappa * mu))
}

.MatClust.inflate <- function(mod, rD) {
  R <- mod$par[["R"]]
  rE <- if(R < rD) (rD + R) else rD
  return(rE)
}
  
  
#' >>>>>>>>>>>>>>>>   CAUCHY CLUSTER PROCESS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#' idiosyncratic parameter eta2 = 4 * scale^2
#'  scale <-  sqrt(mod$par[["eta2"]])/2

.Cauchy.roffspring <- function(n, mod) {
  rate <- mod$par[["eta2"]]/8 
  b <- 1/sqrt(rgamma(n, shape=1/2, rate=rate))
  list(x = b * rnorm(n), y = b * rnorm(n))
}

.Cauchy.hdom <- function(r, mod, rD) {
  scale2 <-  mod$par[["eta2"]]/4
  z <- numeric(length(r))
  above <- (r > rD)
  z[!above] <- 1/(2*pi*scale2)
  z[above]  <- (1/(2*pi*scale2)) * (1 + ((r[above]-rD)^2)/scale2)^(-3/2)
  return(z)
}
  
.Cauchy.Eplus <- function(r, mod, rD) {
  scale2 <-  mod$par[["eta2"]]/4
  mu <- mod$mu
  z <- numeric(length(r))
  above <- (r > rD)
  A <- mu * (rD^2/(2*scale2))
  z[!above] <- A
  z[above]  <- A * (1 + ((r[above]-rD)^2)/scale2)^(-3/2)
  return(z)
}

.Cauchy.rhoplus <- function(r, mod, rD) {
  kappa <- mod$par[["kappa"]]
  z <- numeric(length(r))
  above <- (r > rD)
  z[!above] <- kappa * (1 - exp(-.Cauchy.Eplus(0, mod, rD)))
  z[above]  <- kappa * (1 - exp(-.Cauchy.Eplus(r[above], mod, rD)))
  return(z)
}

.Cauchy.MplusIntegrand <- function(r, mod, rD) {
  2 * pi * r * .Cauchy.rhoplus(r, mod, rD)
}

.Cauchy.Mplus <- function(r, mod, rD, ..., method="q") {
  ## integrand is linear on [0, rD]
  z <- pi * pmin(r, rD)^2 * .Cauchy.rhoplus(0, mod, rD)
  high <- (r > rD)
  if(any(high))
    z[high] <- z[high] + indefinteg(.Cauchy.MplusIntegrand,
                                    r[high],
                                    lower=rD,
                                    mod=mod, rD=rD, method=method)
  return(z)
}

.Cauchy.MplusInf <- function(mod, rD) { .Cauchy.Mplus(Inf, mod, rD) }

.Cauchy.rhoplusplus <- function(r, mod, rD) {
  mod$par[["kappa"]] * .Cauchy.Eplus(r, mod, rD)
}

.Cauchy.Mplusplus <- function(r, mod, rD, ...) {
  mu <- mod$mu
  kappa <- mod$par[["kappa"]]
  lambda <- kappa * mu
  scale2 <-  mod$par[["eta2"]]/4
  distD <- pmax(0, r - rD)
  lambda * pi * rD^2 * ifelse(r <= rD,
                              r^2/(2*scale2),
                              rD^2/(2*scale2) + 1 + (rD * distD/scale2 - 1)/sqrt(1+distD^2/scale2))
}

.Cauchy.MplusplusInf <- function(mod, rD) {
  mu <- mod$mu
  kappa <- mod$par[["kappa"]]
  lambda <- kappa * mu
  scale2 <-  mod$par[["eta2"]]/4
  rdons2 <- rD^2/scale2
  lambda * pi * rD^2 * (1 + sqrt(rdons2) + rdons2/2)
}

.Cauchy.inflate <- function(mod, rD) {
  mu <- mod$mu
  scale2 <- mod$par[["eta2"]]/4
  a <- if(mu == 0) 1 else (1 + (1-exp(-mu))/mu)
  b <- (rD^2)/(2*a*scale2)
  if(b <= 1) return(rD)
  delta <- sqrt(scale2) * (b^(2/3) - 1)
  return(rD + delta)
}
  


#' ........................................................................

.spatstat.clustersim.InfoTable <- list(
  Thomas = list(
    iscompact    = FALSE,
    roffspring   = .Thomas.roffspring,
    hdom         = .Thomas.hdom,
    Eplus        = .Thomas.Eplus,
    rhoplus      = .Thomas.rhoplus,
    Mplus        = .Thomas.Mplus,
    MplusInf     = .Thomas.MplusInf,
    invMplus     = NULL,
    rhoplusplus  = .Thomas.rhoplusplus,
    Mplusplus    = .Thomas.Mplusplus,
    MplusplusInf = .Thomas.MplusplusInf,
    invMplusplus = NULL,
    inflate      = .Thomas.inflate
  ),
  MatClust = list(
    iscompact    = TRUE,
    roffspring   = .MatClust.roffspring,
    hdom         = .MatClust.hdom,
    Eplus        = .MatClust.Eplus,
    rhoplus      = .MatClust.rhoplus,
    Mplus        = .MatClust.Mplus,
    MplusInf     = .MatClust.MplusInf,
    invMplus     = .MatClust.inverseMplus,
    rhoplusplus  = .MatClust.rhoplusplus,
    Mplusplus    = .MatClust.Mplusplus,
    MplusplusInf = .MatClust.MplusplusInf,
    invMplusplus = .MatClust.inverseMplusplus,
    inflate      = .MatClust.inflate
  ),
  Cauchy = list(
    iscompact    = FALSE,
    roffspring   = .Cauchy.roffspring,
    hdom         = .Cauchy.hdom,
    Eplus        = .Cauchy.Eplus,
    rhoplus      = .Cauchy.rhoplus,
    Mplus        = .Cauchy.Mplus,
    MplusInf     = .Cauchy.MplusInf,
    invMplus     = NULL,
    rhoplusplus  = .Cauchy.rhoplusplus,
    Mplusplus    = .Cauchy.Mplusplus,
    MplusplusInf = .Cauchy.MplusplusInf,
    invMplusplus = NULL,
    inflate      = .Cauchy.inflate
  )
)

spatstatClusterSimModelMatch <- function(name, verbose=TRUE) {
  if(!is.character(name) || length(name) != 1)
    stop("Argument must be a single character string", call.=FALSE)
  TheTable <- .spatstat.clustersim.InfoTable
  nama2 <- names(TheTable)
  mat <- pmatch(name, nama2)
  if(!is.null(mat))
    return(nama2[mat])
  if(verbose)
    warning(sQuote(name), "is not supported;",
            "available options are", commasep(sQuote(nama2)))
  return(NULL)
}

spatstatClusterSimInfo <- function(name) {
  known <- spatstatClusterSimModelMatch(name, FALSE)
  if(is.null(known)) return(NULL)
  return(.spatstat.clustersim.InfoTable[[known]])
}
  
