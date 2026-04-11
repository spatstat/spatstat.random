#'
#'  rdiffuse.R
#'
#'  Random diffusion by random walk on raster
#'
#'  $Revision: 1.4 $ $Date: 2026/04/11 02:40:38 $
#'

rdiffuse <- function(X, sigma, ...) {
  UseMethod("rdiffuse")
}

rdiffuse.ppp <- function(X, sigma, ..., connect=8,
                         method=c("C", "interpreted"),
                         unround=TRUE) {
  stopifnot(is.ppp(X))
  method <- match.arg(method)
  ## >>>>. construct transition matrix etc <<<<<<<<<<<<<<<<<<<<<<
  stuff <- resolveHeatTransitions(x=X, sigma=sigma, connect=connect, ...)
  if(!is.null(stuff$sigmaX))
    stop("lagged arrivals are not yet implemented", call.=FALSE)
  ## discretised point pattern - determines the grid 
  Ximage <- stuff$Y
  ## position of each X[i] as serial number in grid
  Xpos <- stuff$Xpos
  ## transition matrix
  P <- stuff$P
  ## k-step transition matrix
  Pk <- stuff$Pk
  ## number of blocks of size k
  Nblock <- stuff$Nblock
  ## remaining number of steps
  Nrump  <- stuff$Nrump
  ## map from serial numbers 'u' to grid positions (NULL if window is rectangle)
  backmap  <- stuff$backmap
  ## convert transition matrices to efficient row-major format
  P <- as(P, "RsparseMatrix")
  Pk <- as(Pk, "RsparseMatrix")
  ## >>>>>>>>>>>>>   run parallel Markov chains <<<<<<<<<<<<<<<<<<<<<
  Z0 <- Xpos
  Zb <- runSparseMarkovChain(Pk, Z0, Nblock,
                             result="last", check=FALSE, method=method)
  Zf <- runSparseMarkovChain(P,  Zb, Nrump, 
                             result="last", check=FALSE, method=method)
  ## >>>>>>>>>>>>>  map pixel serial numbers to spatial positions <<<<<<<
  if(!is.null(backmap))
    Zf <- backmap[Zf]
  xy <- rasterxy.im(Ximage)
  xyZ <- xy[Zf,]
  ## point pattern result
  Z <- as.ppp(xyZ, W=Window(X))
  ## de-discretise
  if(unround) {
    Z <- rUnround(Z, xstep=Ximage$xstep, ystep=Ximage$ystep)
  }
  return(Z)
}

resolveHeatTransitions <-  function(x, sigma=NULL, ..., connect=8,
                                    symmetric=FALSE, sigmaX=NULL, k=1,
                                    weights=NULL, verbose=TRUE) {
  stopifnot(is.ppp(x))
  nX <- npoints(x)
  check.1.integer(k)
  stopifnot(k >= 1)
  check.1.integer(connect)
  if(!(connect %in% c(4,8)))
    stop("connectivity must be 4 or 8", call.=FALSE)
  if(length(weights)) check.nvector(weights, nX) else weights <- NULL

  ## determine initial state for diffusion
  if(delayed <- !is.null(sigmaX)) {
    #' smoothing bandwidths attributed to each data point
    check.nvector(sigmaX, nX)
    stopifnot(all(is.finite(sigmaX)))
    stopifnot(all(sigmaX >= 0))
    if(is.null(sigma)) sigma <- max(sigmaX) else check.1.real(sigma)
    #' sort in decreasing order of bandwidth
    osx <- order(sigmaX, decreasing=TRUE)
    sigmaX <- sigmaX[osx]
    x <- x[osx]
    #' discretise window
    W <- do.call.matched(as.mask,
                         resolve.defaults(list(...),
                                          list(w=Window(x))))
    #' initial state is zero
    Y <- as.im(W, value=0)
    #' discretised coordinates
    Xpos <- nearest.valid.pixel(x$x, x$y, Y)
  } else {
    #' pixellate pattern
    Y <- pixellate(x, ..., weights=weights, preserve=TRUE, savemap=TRUE)
    Xpos <- attr(Y, "map")
    osx <- NULL
  } 

  #' validate sigma
  if(is.im(sigma)) {
    # ensure Y and sigma are on the same grid
    A <- harmonise(Y=Y, sigma=sigma)
    Y     <- A$Y
    sigma <- A$sigma
  } else if(is.function(sigma)) {
    sigma <- as.im(sigma, as.owin(Y))
  } else {
    sigma <- as.numeric(sigma)
    check.1.real(sigma)
  }

  #' normalise as density 
  pixelarea <- with(Y, xstep * ystep)
  Y <- Y / pixelarea
  v <- as.matrix(Y)
  #' initial state
  u <- as.vector(v)

  #' map (row, col) to serial number 
  serial <- matrix(seq_len(length(v)), nrow(v), ncol(v))
  Xpos <- serial[as.matrix(as.data.frame(Xpos))]
  
  #' symmetric random walk?
  if(symmetric) {
    asprat <- with(Y, ystep/xstep)
    if(abs(asprat-1) > 0.01)
      warning(paste("Symmetric random walk on a non-square grid",
                    paren(paste("aspect ratio", asprat))),
              call.=FALSE)
  }
  #' determine appropriate jump probabilities & time step
  pmax <- 1/(connect+1) # maximum permitted jump probability
  xstep <- Y$xstep
  ystep <- Y$ystep
  minstep <- min(xstep, ystep)
  if(symmetric) {
    #' all permissible transitions have the same probability 'pjump'.
    #' Determine Nstep, and dt=sigma^2/Nstep, such that
    #' Nstep >= 16 and M * pjump * minstep^2 = dt
    M <- if(connect == 4) 2 else 6
    Nstep <- max(16, ceiling(max(sigma)^2/(M * pmax * minstep^2)))    
    dt <- sn <- (sigma^2)/Nstep
    px <- py <- pxy <- sn/(M * minstep^2)
  } else {
    #' px is the probability of jumping 1 step to the right
    #' py is the probability of jumping 1 step up
    #' if connect=4, horizontal and vertical jumps are exclusive.
    #' if connect=8, horizontal and vertical increments are independent
    #' Determine Nstep, and dt = sigma^2/Nstep, such that
    #' Nstep >= 16 and 2 * pmax * minstep^2 = dt
    Nstep <- max(16, ceiling(max(sigma)^2/(2 * pmax * minstep^2)))
    dt <- sn <- (sigma^2)/Nstep
    px <- sn/(2 * xstep^2)
    py <- sn/(2 * ystep^2)
    if(max(px) > pmax) stop("Internal error: px exceeds pmax")
    if(max(py) > pmax) stop("Internal error: py exceeds pmax")
    if(connect == 8) pxy <- px * py
  }
  #' arrival times
  if(!is.null(sigmaX)) 
    iarrive <- pmax(1, pmin(Nstep, Nstep - round((sigmaX^2)/sn)))
  #' construct adjacency matrices
  dimv <- dim(v)
  my <- gridadjacencymatrix(dimv, across=FALSE, down=TRUE, diagonal=FALSE)
  mx <- gridadjacencymatrix(dimv, across=TRUE,  down=FALSE, diagonal=FALSE)
  if(connect == 8)
    mxy <- gridadjacencymatrix(dimv, across=FALSE,  down=FALSE, diagonal=TRUE)
  #' restrict to window
  if(anyNA(u)) {
    ok <- !is.na(u)
    u <- u[ok]
    #' adjust serial numbers
    Xpos <- cumsum(ok)[Xpos]
    backmap <- which(ok)
    mx <- mx[ok,ok,drop=FALSE]
    my <- my[ok,ok,drop=FALSE]
    if(connect == 8) 
      mxy <- mxy[ok,ok,drop=FALSE]
    if(is.im(sigma)) {
      px <- px[ok]
      py <- py[ok]
      if(connect == 8) 
        pxy <- pxy[ok]
    }
  } else {
    ok <- TRUE
    backmap <- NULL
    if(is.im(sigma)) {
      px <- px[]
      py <- py[]
      if(connect == 8) pxy <- pxy[]
    }
  }
  #' construct iteration matrix
  if(connect == 4) {
    P <- px * mx + py * my
  } else {
    P <- px * (1 - 2 * py) * mx + py * (1 - 2 * px) * my + pxy * mxy
  }

  #' construct one-step transition probability matrix
  if(any(P < 0)) 
    stop("Internal error: negative jump probabilities", call.=FALSE)
  totaljump <- rowSums(P)
  if(max(totaljump) > 1)
    stop("Internal error: jump probability exceeds 1", call.=FALSE)
  diag(P) <- 1 - totaljump
  
  #' k-step transition probabilities
  k <- as.integer(k)
  if(k == 1) {
    Pk = P
  } else {
    Pk <- P
    for(j in 2:k) Pk <- Pk %*% P
  } 
  Nstep <- as.integer(Nstep)
  Nblock <- Nstep/k
  Nrump  <- Nstep - Nblock * k

  return(list(Y       = Y,          # discretised point pattern (-> grid)
              u       = u,          # initial state vector
              Xpos    = Xpos,       # position of each X[i] as serial number
              backmap = backmap,    # map from serial numbers to grid positions
              sigma   = sigma,      # numeric or image
              sigmaX  = sigmaX,     # numeric vector or NULL
              osx     = osx,        # integer vector or NULL
              weights = weights,    # numeric vector or NULL
              P       = P,          # transition matrix (sparse)
              Pk      = Pk,         # k-step transition matrix (sparse)
              k       = k,          # 
              Nstep   = Nstep,      # total number of steps
              Nblock  = Nblock,     # number of blocks of size k steps
              Nrump   = Nrump,      # Nstep - k * Nblock
              dx      = xstep,      # pixel width
              dy      = ystep,      # pixel height
              dt      = dt))        # time increment
}
