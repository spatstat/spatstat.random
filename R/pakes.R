##  pakes.R
##
## The probability distribution that satisfies the distributional equivalence
##              X == U^(1/zeta) (1 + X)
## where X, U are independent r.v.'s and U is uniform [0,1]
##
## Solution due to A.G. Pakes presented in
##     Baddeley, Moller & Pakes (2008)
##     Properties of residuals for spatial point processes,
##     Annals of the Institute of Statistical Mathematics 60 (2008) 627-649
##
## Implementation (c) Adrian Baddeley 2021
## GNU Public Licence >= 2.0


pakesCalc <- function(zeta, dx=0.001, xmax=5) {
  check.1.real(zeta)
  stopifnot(zeta > 0)
  EulerGamma <- -digamma(1) # Euler-Mascheroni constant
  C <- exp(-EulerGamma * zeta)/gamma(1 + zeta)
  ## cat(paste("C=", C, "\n"))
  x <- seq(0, xmax, by=dx)
  xzeta <- x^zeta
  inv1x1zeta <- 1/(1 + x)^(1+zeta)
  nx <- length(x)
  shifted <- seq_len(nx) - floor(1/dx)
  shifted[shifted < 1] <- NA
  Fx <- numeric(nx)
  ceilx <- ceiling(x)
  slice1 <- (ceilx == 1)
  Fx[slice1] <- C * xzeta[slice1]
  nmax <- ceiling(xmax)
  if(nmax > 1) {
    for(k in 2:nmax) {
      integrand <- Fx * inv1x1zeta
      indefInteg <- zeta * dx * cumsum(integrand)
      ## plot(x, indefInteg, type="l")
      shiftedIndefInt <- ifelse(is.na(shifted), 0, indefInteg[shifted])
      ## plot(x, shiftedIndefInt, type="l")
      slicek <- (ceilx == k)
      Fx[slicek] <- xzeta[slicek] * (C - shiftedIndefInt[slicek])
      ## suppress numerical glitch close to F(x) = 1
      if(any(exceed <- (Fx >= 1))) {
        last <- min(which(exceed))
        Fx[last:length(Fx)] <- 1
        break;
      }
    }
  }
  ## suppress numerical glitches
  if(any(dip <- (Fx < cummax(Fx)))) {
    ok <- !dip
    Fx[dip] <- approx(x=x[ok], y=Fx[ok],
                      xout=x[dip],
                      rule=2, yright=1)$y
  }
  ## return
  data.frame(x=x, Fx=Fx)
}

ppakes <- function(q, zeta) {
  q <- as.numeric(q)
  zeta <- as.numeric(zeta)
  if(length(zeta) > 1)
    return(mapply(ppakes, q=q, zeta=zeta))
  a <- pakesCalc(zeta, xmax=max(ceiling(q)))
  Fpakes <- approxfun(a$x, a$Fx, rule=2)
  p <- Fpakes(q)
  return(p)
}

dpakes <- function(x, zeta) {
  x <- as.numeric(x)
  zeta <- as.numeric(zeta)
  if(length(zeta) > 1) 
    return(mapply(dpakes, x=x, zeta=zeta))
  C <- exp(digamma(1) * zeta)/gamma(1 + zeta)
  y <- numeric(length(x))
  if(any(low <- (x <= 1))) 
    y[low] <- C * zeta * x[low]^(zeta-1)
  if(any(high <- !low)) {
    xhigh <- x[high]
    a <- pakesCalc(zeta, xmax=max(xhigh))
    Fxhigh <- approx(x=a$x, y=a$Fx,
                     xout=xhigh,
                     rule=2, yright=1)$y
    Fxhigh1 <- approx(x=a$x, y=a$Fx,
                     xout=xhigh - 1,
                     rule=2, yright=1)$y
    y[high] <- (zeta/xhigh) * (Fxhigh - Fxhigh1)
  }
  return(y)
}

qpakes <- function(p, zeta) {
  p <- as.numeric(p)
  zeta <- as.numeric(zeta)
  if(length(zeta) > 1) 
    return(mapply(qpakes, p=p, zeta=zeta))
  ## find 'xmax' such that F(xmax) = 1 within numerical error
  zetabreaks <- c(0, 0.1, 0.5, 1, 2, 10, 20, Inf)
  xmaxvalues <- c(  3,  6,  8,  10, 20, 30, 1.5 * zeta)
  xmax <- xmaxvalues[findInterval(zeta, zetabreaks, all.inside=TRUE)]
  ## compute CDF
  a <- pakesCalc(zeta, xmax=xmax)
  a <- a[!duplicated(a$Fx), ]
  ## invert
  q <- approx(x=a$Fx, y=a$x,
              xout=p,
              rule=2, yleft=0, yright=1)$y
  return(q)
}

rpakes <- function(n, zeta) {
  qpakes(runif(n), zeta)
}
