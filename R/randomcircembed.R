#'
#' randomcircembed.R
#'
#' Circulant embedding for simulating Gaussian random field
#'
#' Originally derived from code provided in:
#'    Tilman M. Davies and David Bryant
#'    On Circulant Embedding for Gaussian Random Fields in R
#'    Journal of Statistical Software 55 (2013) issue 9
#'    DOI 10.18637/jss.v055.i09
#'
#' Modified by Adrian Baddeley
#'
#'   $Revision: 1.7 $ $Date: 2024/02/26 05:33:18 $
#' 
#'   Copyright (c) 2023 Tilman M. Davies, David Bryant and Adrian Baddeley
#'   GNU Public Licence (>= 2.0)

rGRFcircembed <- local({

  rGRFcircembed <- function(W=owin(), mu=0, var=1, corrfun, ...,
                            nsim=1, drop=TRUE) {
    rCircEmbedEngine(W=W, mu=mu, var=var, corrfun=corrfun,
                     ..., nsim=nsim, drop=drop)
  }

  rCircEmbedEngine <- function(W=owin(), mu=0, var=1, corrfun,
                               ...,
                               dimyx=NULL, eps=NULL, xy=NULL,
                               rule.eps = c("adjust.eps", 
                                            "grow.frame", "shrink.frame"),
                               warn=TRUE, maxrelerr=1e-6,
                               nsim=1, drop=TRUE) {
    ## determine discretisation
    M <- as.mask(W, dimyx=dimyx, eps=eps, xy=xy, rule.eps=rule.eps)
    if(needclip <- !is.rectangle(W))
      outsideM <- !as.vector(as.matrix(M))
    gp <- grid.prep(W, ncol(M), nrow(M))
    ## convert 'mu' to an image on same raster (unless it is a single number)
    if(is.function(mu) || is.im(mu))
      mu <- as.im(mu, W=M, ...)
    ## calculate covariance matrix
    Sigma <- covariance.prep(gp=gp,
                             var=var,
                             corr.func=corrfun,
                             wrap=TRUE,
                             ...)
    ## generate (centred) pixel values for all nsim realisations
    z <- generate.normals(nsim, Sigma, warn, maxrelerr)
    ## reshape
    z <- array(as.numeric(z), dim=c(gp$N.ext, gp$M.ext, nsim))
    z <- z[1:(gp$N), 1:(gp$M), , drop=FALSE]
    ## pack up and add 'mu'
    results <- vector(mode="list", length=nsim)
    R0 <- as.im(0, W=M)
        
    for(isim in 1:nsim) {
      zmat <- z[,,isim]
      if(needclip) zmat[outsideM] <- NA
      R <- R0
      R[drop=FALSE] <- zmat
      results[[isim]] <- R+mu
    }
    results <- simulationresult(results, nsim, drop)
    return(results)
  }

  grid.prep <- function(W, M, N, ext = 2) {
    cell.width <- diff(W$xrange)/M
    cell.height <- diff(W$yrange)/N
    
    mgrid <- seq(W$xrange[1], W$xrange[2], by = cell.width)
    ngrid <- seq(W$yrange[1], W$yrange[2], by = cell.height)
    mcens <- (mgrid + 0.5 * cell.width)[-(M + 1)]
    ncens <- (ngrid + 0.5 * cell.height)[-(N + 1)]
    
    if (ext <= 1) {
      mgrid.ext <- ngrid.ext <- mcens.ext <- ncens.ext <-
        M.ext <- N.ext <- NULL
    } else {
      M.ext <- ext * M
      N.ext <- ext * N
      mgrid.ext <- seq(W$xrange[1],
                       W$xrange[2] + (ext - 1) * diff(W$xrange),
                       by = cell.width)
      ngrid.ext <- seq(W$yrange[1],
                       W$yrange[2] + (ext - 1) * diff(W$yrange),
                       by = cell.height)
      mcens.ext <- (mgrid.ext + 0.5 * cell.width)[-(M.ext + 1)]
      ncens.ext <- (ngrid.ext + 0.5 * cell.height)[-(N.ext + 1)]
    }
  
    return(list(M = M,
                N = N,
                mgrid = mgrid,
                ngrid = ngrid,
                mcens = mcens,
                ncens = ncens, 
                cell.width = cell.width,
                cell.height = cell.height,
                M.ext = M.ext,
                N.ext = N.ext, 
                mgrid.ext = mgrid.ext,
                ngrid.ext = ngrid.ext,
                mcens.ext = mcens.ext,
                ncens.ext = ncens.ext))
  }

  ## get.EIG <- function(Sigma){
  ##   eigs <- eigen(Sigma,symmetric=TRUE)
  ##   lambda <- eigs$values
  ##   lambda[lambda < .Machine$double.eps] <- 0
  ##  result <- eigs$vectors %*% diag(sqrt(lambda))
  ##  return(result)
  ## }

  covariance.prep <- function(gp, var, corr.func, wrap=FALSE, ...){
    if(!wrap){
      cent <- expand.grid(gp$mcens,gp$ncens)
      mmat <- matrix(rep(cent[,1],gp$M*gp$N),gp$M*gp$N,gp$M*gp$N)
      nmat <- matrix(rep(cent[,2],gp$M*gp$N),gp$M*gp$N,gp$M*gp$N)
      D <- sqrt((mmat-t(mmat))^2+(nmat-t(nmat))^2)
      covmat <- var*corr.func(D,...)
      return(covmat)
    } else {
      Rx <- gp$M.ext*gp$cell.width
      Ry <- gp$N.ext*gp$cell.height
      m.abs.diff.row1 <- abs(gp$mcens.ext[1]-gp$mcens.ext)
      m.diff.row1 <- pmin(m.abs.diff.row1,Rx-m.abs.diff.row1)
      n.abs.diff.row1 <- abs(gp$ncens.ext[1]-gp$ncens.ext)
      n.diff.row1 <- pmin(n.abs.diff.row1,Ry-n.abs.diff.row1)
      cent.ext.row1 <- expand.grid(m.diff.row1,n.diff.row1)
      D.ext.row1 <- matrix(sqrt(cent.ext.row1[,1]^2+cent.ext.row1[,2]^2),
                           gp$M.ext,gp$N.ext)
      C.tilde <- var*corr.func(D.ext.row1,...)
      return(C.tilde)
    }
  }

  generate.normals <- function(nsim, Sigma, warn=TRUE, maxrelerr=1e-6) {
    ## calculate fft of covariance matrix
    refft <- Re(fft(Sigma, inverse=TRUE))
    if(min(refft) < 0) {
      ## in theory this is impossible because Sigma is positive definite, 
      ## but small negative values do occur
      if(warn) {
        ra <- range(refft)
        if(-ra[1]/ra[2] > maxrelerr) {
          bad <- (refft < 0)
          warning(paste(sum(bad), "out of", length(bad), "terms",
                        paren(percentage(mean(bad))),
                        "in FFT calculation of matrix square root",
                        "were negative, and were set to zero.",
                        "Range:", prange(signif(ra, 3))),
                  call.=FALSE)
        }
      }
      refft <- pmax(0, refft)
    }
    ## fft of square root of matrix Sigma
    sqrtFFTsigma <- sqrt(refft)
    ## set up simulation 
    nr <- nrow(Sigma)
    nc <- ncol(Sigma)
    ncell <- prod(dim(Sigma))
    sqrtncell <- sqrt(ncell)
    ## run
    realisations <- matrix(, ncell, nsim)
    noise <- array(rnorm(ncell * nsim), dim=c(nr, nc, nsim))
    for(i in 1:nsim) {
      field <- sqrtFFTsigma * fft(noise[,,i])/sqrtncell
      realisations[,i] <- Re(fft(field,inverse=TRUE)/sqrtncell)
    }
    return(realisations)
  }

  rGRFcircembed

})

