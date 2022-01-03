#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.random
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.random)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
#'  tests/randoms.R
#'   Further tests of random generation code
#'  $Revision: 1.14 $ $Date: 2021/09/09 10:02:00 $


local({
  if(FULLTEST) {
    #' cases not covered in examples
    A <- runifdisc(6, nsim=2)
    A <- runifpoispp(5, nsim=2)
    A <- runifpoispp(0, nsim=2)
    A <- rSSI(0.05, 6, nsim=2)
    A <- rSSI(0.05, 10, win=square(c(-0.5, 1.5)), x.init=A[[1]], nsim=2)  
    A <- rstrat(nx=4, nsim=2)
    A <- rcell(square(1), nx=5, nsim=2)
  }
  if(ALWAYS) { # involves C code etc
    A <- rthin(cells, P=0.5, nsim=2)
    A <- rthin(cells, runif(42))
    A <- rthin(cells[FALSE], P=0.5, nsim=2)
  }
  f <- function(x,y) { 10*x }
  Z <- as.im(f, square(1))
  if(ALWAYS) {
    A <- rpoint(n=6, f=f, fmax=10, nsim=2)
    A <- rpoint(n=6, f=Z, fmax=10, nsim=2)
    A <- rpoint(n=0, f=f, fmax=10, nsim=2)
    A <- rpoint(n=0, f=Z, fmax=10, nsim=2)

    op <- spatstat.options(fastpois=FALSE)
    A <- runifpoispp(5, nsim=2)
    A <- rpoispp(Z)
    spatstat.options(op)
  }
  if(FULLTEST) {
    b3 <- box3(c(0,1))
    b4 <- boxx(c(0,1), c(0,1), c(0,1), c(0,1))
    b5 <- c(0, 2, 0, 2)
    X <- rMaternInhibition(2, kappa=20, r=0.1, win=b3)
    Y <- rMaternInhibition(2, kappa=20, r=0.1, win=b4)
    Y <- rMaternInhibition(2, kappa=20, r=0.1, win=b5, nsim=2)

    X <- rSSI(0.05, 6)
    Y <- rSSI(0.05, 6, x.init=X) # no extra points

    Z <- rlabel(finpines)
  }

  f1 <- function(x,y){(x^2 + y^3)/10}
  f2 <- function(x,y){(x^3 + y^2)/10}
  ZZ <- solist(A=as.im(f1, letterR),
               B=as.im(f2, letterR))
  g <- function(x,y,m){(10+as.integer(m)) * (x^2 + y^3)}
  if(FULLTEST) {
    XX <- rmpoispp(ZZ, nsim=3)
    YY <- rmpoint(10, f=ZZ, nsim=3)
    VV <- rpoint.multi(10, f=g,
                       marks=factor(sample(letters[1:3], 10, replace=TRUE)),
                       nsim=3)
  }
  if(ALWAYS) { # depends on C code
    L <- edges(letterR)
    E <- runifpoisppOnLines(5, L)
    G <- rpoisppOnLines(ZZ, L)
    G2 <- rpoisppOnLines(list(A=f1, B=f2), L, lmax=max(sapply(ZZ, max)))
  }

  if(FULLTEST) {
    #' cluster models + bells + whistles
    X <- rThomas(10, 0.2, 5, saveLambda=TRUE)
    if(is.null(attr(X, "Lambda")))
      stop("rThomas did not save Lambda image")
    Y <- rThomas(0, 0.2, 5, saveLambda=TRUE)
    if(is.null(attr(Y, "Lambda")))
      stop("rThomas did not save Lambda image when kappa=0")
    X <- rMatClust(10, 0.05, 4, saveLambda=TRUE)
    X <- rCauchy(30, 0.01, 5, saveLambda=TRUE)
    X <- rVarGamma(30, 2, 0.02, 5, saveLambda=TRUE)
    Z <- as.im(function(x,y){ 5 * exp(2 * x - 1) }, owin())
    Y <- rThomas(10, 0.2, Z, saveLambda=TRUE)
    Y <- rMatClust(10, 0.05, Z, saveLambda=TRUE)
    Y <- rCauchy(30, 0.01, Z, saveLambda=TRUE)
    Y <- rVarGamma(30, 2, 0.02, Z, saveLambda=TRUE)
  }

  if(FULLTEST) {
    #' perfect simulation code infrastructure
    expandwinPerfect(letterR, 2, 3)

    #' trivial cases of random generators for ppx
    B4 <- boxx(0:1, 0:1, 0:1, 0:1)
    Z0 <- runifpointx(0, domain=B4, nsim=2)
    Z1 <- runifpointx(1, domain=B4, nsim=2)
  }

})
reset.spatstat.options()


