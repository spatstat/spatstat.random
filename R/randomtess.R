#
# randomtess.R
#
# Random tessellations
#
# $Revision: 1.8 $  $Date: 2023/10/20 03:38:34 $
#

# Poisson line tessellation

rpoislinetess <- function(lambda, win=owin()) {
  win <- as.owin(win)
  ## determine circumcircle
  xr <- win$xrange
  yr <- win$yrange
  xmid <- mean(xr)
  ymid <- mean(yr)
  width <- diff(xr)
  height <- diff(yr)
  rmax <- sqrt(width^2 + height^2)/2
  boundbox <- owin(xmid + c(-1,1) * rmax, ymid + c(-1,1) * rmax)
  ## generate poisson lines through circumcircle
  n <- rpois(1, lambda * 2 * pi * rmax)
  if(n == 0) {
    ## single tile
    if(is.mask(win)) {
      M <- as.im(factor(1), W=win)
      return(tess(image=M))
    } else {
      return(tess(tiles=list(win)))
    }
  }
  theta <- runif(n, max= 2 * pi)
  p <- runif(n, max=rmax) + xmid * cos(theta) + ymid * sin(theta)
  Y <- infline(p=p, theta=theta)
  # form the induced tessellation in bounding box
  Z <- chop.tess(boundbox, Y)
  # clip to window
  Z <- intersect.tess(Z, win)
  attr(Z, "lines") <- Y
  return(Z)
}

rMosaicSet <- function(X, p=0.5) {
  stopifnot(is.tess(X))
  Y <- tiles(X)
  Y <- Y[runif(length(Y)) < p]
  if(length(Y) == 0)
    return(NULL)
  Z <- NULL
  for(i in seq_along(Y))
    Z <- union.owin(Z, Y[[i]])
  return(Z)
}

rMosaicField <- function(X,
                    rgen=function(n) { sample(0:1, n, replace=TRUE)},
                    ..., 
                    rgenargs=NULL ) {
  stopifnot(is.tess(X))
  Y <- as.im(X, ...)
  ntiles <- length(levels(Y))
  values <- do.call(rgen, append(list(ntiles),rgenargs))
  Z <- eval.im(values[as.integer(Y)])
  return(Z)
}

