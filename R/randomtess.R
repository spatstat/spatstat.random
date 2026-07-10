#
# randomtess.R
#
# Random tessellations
#
# $Revision: 1.10 $  $Date: 2026/07/10 07:02:03 $
#

# Poisson line tessellation

rpoislinetess <- function(lambda, win=owin()) {
  win <- as.owin(win)
  check.1.real(lambda)
  stopifnot(lambda >= 0)
  ## determine circumcircle
  xr <- win$xrange
  yr <- win$yrange
  xmid <- mean(xr)
  ymid <- mean(yr)
  width <- diff(xr)
  height <- diff(yr)
  rmax <- sqrt(width^2 + height^2)/2
  boundbox <- owinInternalRect(xmid + c(-1,1) * rmax, ymid + c(-1,1) * rmax)
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

rpoisDirichletTess <- function(lambda, win=owin()) {
  win <- as.owin(win)
  B <- Frame(win)
  X <- rpoispp(lambda, win=B)
  nX <- npoints(X)
  ## determine maximum possible distance from B at which
  ## a random point outside B could affect the result inside B
  if(nX == 0) {
    dmax <- diameter(B)
  } else {
    A <- tiles(dirichlet(X))
    dmax <- 0
    for(i in 1:nX) 
      dmax <- max(dmax, nncross(vertices(A[[i]]), X[i], what="dist"))
  }
  ## generate extended realisation of Poisson process
  Bplus    <- grow.rectangle(B, dmax)
  Boutside <- setminus.owin(Bplus, B)
  Xoutside <- rpoispp(lambda, win=Boutside)
  Xplus <- superimpose(X, Xoutside, W=Bplus)
  ## now compute Dirichlet tessellation
  Dplus <- dirichlet(Xplus)
  ## finally clip it 
  tilesDW <- tiles(intersect.tess(Dplus, win, keepempty=TRUE))
  retain <- !sapply(tilesDW, is.empty)
  Xplus <- Xplus[retain]
  D <- tess(tiles=tilesDW[retain], window=win)
  ## return with attributes
  Xplus <- Xplus[ripras(Xplus, shape="rectangle")]
  attr(D, "X") <- Xplus
  attr(D, "dmax") <- dmax
  return(D)
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

