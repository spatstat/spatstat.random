##        clusterinfo.R
## 
##   Lookup table of explicitly-known K functions and pcf
##   and algorithms for computing sensible starting parameters
##
##   $Revision: 1.37 $ $Date: 2022/02/11 13:23:47 $


.Spatstat.ClusterModelInfoTable <- 
  list(
       Thomas=list(
         ## Thomas process: old par = (kappa, sigma2) (internally used everywhere)
         ## Thomas process: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Thomas process", # In modelname field of mincon fv obj.
         descname = "Thomas process", # In desc field of mincon fvb obj.
         modelabbrev = "Thomas process", # In fitted obj.
         printmodelname = function(...) "Thomas process", # Used by print.kppm
         parnames = c("kappa", "sigma2"),
         clustargsnames = NULL,
         checkpar = function(par, old = TRUE, ...){
           ## 'par' is in either format
           if(is.null(par))
             par <- c(kappa=1,scale=1)
           if(any(par<=0))
             stop("par values must be positive.", call.=FALSE)
           nam <- check.named.vector(par, c("kappa","sigma2"),
                                     onError="null")
           if(is.null(nam)) {
             check.named.vector(par, c("kappa","scale"))
             names(par)[2L] <- "sigma2"
             par[2L] <- par[2L]^2
           }
           if(!old){
             names(par)[2L] <- "scale"
             par[2L] <- sqrt(par[2L])
           }
           return(par)
         },
         checkclustargs = function(margs, old = TRUE, ...) list(),
         resolvedots = function(...){
           return(list(...))
         },
         # density function for the distance to offspring
         ddist = function(r, scale, ...) {
           ## 'scale' is generic format
           2 * pi * r * dnorm(r, 0, scale)/sqrt(2*pi*scale^2)
         },
         ## Practical range of clusters
         range = function(..., par=NULL, thresh=NULL){
           ## 'par' is in generic format
           scale <- retrieve.param("scale", "sigma", ..., par=par)
           if(!is.null(thresh)){
             ## The squared length of isotropic Gaussian (sigma)
             ## is exponential with mean 2 sigma^2
             rmax <- scale * sqrt(2 * qexp(thresh, lower.tail=FALSE))
           } else {
             rmax <- 4*scale
           }
           return(rmax)
         },
         kernel = function(par, rvals, ...) {
           ## 'par' is in idiosyncratic ('old') format
           scale <- sqrt(par[2L])
           dnorm(rvals, 0, scale)/sqrt(2*pi*scale^2)
         },
         isPCP=TRUE,
         ## K-function
         K = function(par,rvals, ...){
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           pi*rvals^2+(1-exp(-rvals^2/(4*par[2L])))/par[1L]
         },
         ## pair correlation function
         pcf= function(par,rvals, ...){
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           1 + exp(-rvals^2/(4 * par[2L]))/(4 * pi * par[1L] * par[2L])
         },
         ## gradient of pcf (contributed by Chiara Fend)
         Dpcf= function(par,rvals, ...){
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0)){
             dsigma2 <- rep.int(Inf, length(rvals))
             dkappa <- rep.int(Inf, length(rvals))
           } else {
             dsigma2 <- exp(-rvals^2/(4 * par[2L])) * (rvals/(4^2 * pi * par[1L] * par[2L]^3) - 1/(4 * pi * par[1L] * par[2L]^2))
             dkappa <- -exp(-rvals^2/(4 * par[2L]))/(4 * pi * par[1L]^2 * par[2L])
           }
           out <- rbind(dkappa, dsigma2)
           rownames(out) <- c("kappa","sigma2")
           return(out)
         },
         ## sensible starting parameters
         selfstart = function(X) {
           ## return 'par' in idiosyncratic ('old') format
           kappa <- intensity(X)
           sigma2 <- 4 * mean(nndist(X))^2
           c(kappa=kappa, sigma2=sigma2)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           sigma <- sqrt(par[["sigma2"]])
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, sigma=sigma, mu=mu)
         }
       ),
       ## ...............................................
       MatClust=list(
         ## Matern Cluster process: old par = (kappa, R) (internally used everywhere)
         ## Matern Cluster process: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Matern cluster process", # In modelname field of mincon fv obj.
         descname = "Matern cluster process", # In desc field of mincon fv obj.
         modelabbrev = "Matern cluster process", # In fitted obj.
         printmodelname = function(...) "Matern cluster process", # Used by print.kppm
         parnames = c("kappa", "R"),
         clustargsnames = NULL,
         checkpar = function(par, old = TRUE, ...){
           ## 'par' is in either format
           if(is.null(par))
                 par <- c(kappa=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.", call.=FALSE)
             nam <- check.named.vector(par, c("kappa","R"), onError="null")
             if(is.null(nam)) {
               check.named.vector(par, c("kappa","scale"))
               names(par)[2L] <- "R"
             }
             if(!old){
                 names(par)[2L] <- "scale"
             }
             return(par)
         },
         # density function for the distance to offspring
         ddist = function(r, scale, ...) {
           ## 'scale' is generic format
           ifelse(r>scale, 0, 2 * r / scale^2)
         },
         ## Practical range of clusters
         range = function(..., par=NULL, thresh=NULL){
           ## 'par' is in generic format
           scale <- retrieve.param("scale", "R", ..., par=par)
           if(!is.null(thresh))
             warn.once("thresh.Matern",
                       "Argument", sQuote("thresh"), "is ignored for Matern Cluster model")
           return(scale)
         },
         checkclustargs = function(margs, old = TRUE, ...) list(),
         resolvedots = function(...){
           return(list(...))
         },
         kernel = function(par, rvals, ...) {
           ## 'par' is in idiosyncratic ('old') format
           scale <- par[2L]
           ifelse(rvals>scale, 0, 1/(pi*scale^2))
         },
         isPCP=TRUE,
         K = function(par,rvals, ..., funaux){
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           kappa <- par[1L]
           R <- par[2L]
           Hfun <- funaux$Hfun
           y <- pi * rvals^2 + (1/kappa) * Hfun(rvals/(2 * R))
           return(y)
         },
         pcf= function(par,rvals, ..., funaux){
           ## 'par' is in idiosyncratic ('old') format
             if(any(par <= 0))
               return(rep.int(Inf, length(rvals)))
             kappa <- par[1L]
             R <- par[2L]
             g <- funaux$g
             y <- 1 + (1/(pi * kappa * R^2)) * g(rvals/(2 * R))
             return(y)
         },
         Dpcf= function(par,rvals, ..., funaux){
           ## 'par' is in idiosyncratic ('old') format
           kappa <- par[1L]
           R <- par[2L]
           g <- funaux$g
           gprime <- funaux$gprime
           if(any(par <= 0)){
             dkappa <- rep.int(Inf, length(rvals))
             dR <- rep.int(Inf, length(rvals))
           } else {
             dkappa <- -g(rvals/(2 * R)) / (pi * kappa^2 * R^2)
             dR <- -2*g(rvals/(2 * R))/(pi * kappa * R^3) - (1/(pi * kappa * R^2)) * gprime(rvals/(2 * R))*rvals/(2*R^2)
           }
           out <- rbind(dkappa, dR)
           rownames(out) <- c("kappa","R")
           return(out)
         },         
         funaux=list(
           Hfun=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 1
             z <- zz[ok]
             h[ok] <- 2 + (1/pi) * (
                                    (8 * z^2 - 4) * acos(z)
                                    - 2 * asin(z)
                                    + 4 * z * sqrt((1 - z^2)^3)
                                    - 6 * z * sqrt(1 - z^2)
                                    )
             return(h)
           },
           DOH=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 0
             z <- zz[ok]
             h[ok] <- (16/pi) * (z * acos(z) - (z^2) * sqrt(1 - z^2))
             return(h)
           },
           ## g(z) = DOH(z)/z has a limit at z=0.
           g=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 0
             z <- zz[ok]
             h[ok] <- (2/pi) * (acos(z) - z * sqrt(1 - z^2))
             return(h)
           },
           gprime=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 0
             z <- zz[ok]
             h[ok] <- -(2/pi) * 2 * sqrt(1 - z^2)
             return(h)
           }),
         ## sensible starting paramters
         selfstart = function(X) {
           ## return 'par' in idiosyncratic ('old') format
           kappa <- intensity(X)
           R <- 2 * mean(nndist(X)) 
           c(kappa=kappa, R=R)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           R     <- par[["R"]]
           mu    <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA           
           c(kappa=kappa, R=R, mu=mu)
         }
         ),
       ## ...............................................
       Cauchy=list(
         ## Neyman-Scott with Cauchy clusters: old par = (kappa, eta2) (internally used everywhere)
         ## Neyman-Scott with Cauchy clusters: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Neyman-Scott process with Cauchy kernel", # In modelname field of mincon fv obj.
         descname = "Neyman-Scott process with Cauchy kernel", # In desc field of mincon fv obj.
         modelabbrev = "Cauchy process", # In fitted obj.
         printmodelname = function(...) "Cauchy process", # Used by print.kppm
         parnames = c("kappa", "eta2"),
         clustargsnames = NULL,
         checkpar = function(par, old = TRUE, ...){
           ## 'par' is in either format
           if(is.null(par))
                 par <- c(kappa=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.", call.=FALSE)
             nam <- check.named.vector(par, c("kappa","eta2"), onError="null")
             if(is.null(nam)) {
                 check.named.vector(par, c("kappa","scale"))
                 names(par)[2L] <- "eta2"
                 par[2L] <- (2*par[2L])^2
             }
             if(!old){
                 names(par)[2L] <- "scale"
                 par[2L] <- sqrt(par[2L])/2
             }
             return(par)
         },
         checkclustargs = function(margs, old = TRUE, ...) list(),
         resolvedots = function(...){
           return(list(...))
         },
         # density function for the distance to offspring
         ddist = function(r, scale, ...) {
           ## 'scale' is generic format
           r/(scale^2) *  (1 + (r / scale)^2)^(-3/2)
         },
         ## Practical range of clusters
         range = function(..., par=NULL, thresh=0.01){
           ## 'par' is in generic format
           thresh <- as.numeric(thresh %orifnull% 0.01)
           scale <- retrieve.param("scale", character(0), ..., par=par)
           ## integral of ddist(r) dr is 1 - (1+(r/scale)^2)^(-1/2)
           ## solve for integral = 1-thresh:
           rmax <- scale * sqrt(1/thresh^2 - 1)
           return(rmax)
         },
         kernel = function(par, rvals, ...) {
           ## 'par' is in idiosyncratic ('old') format
           scale <- sqrt(par[2L])/2
           1/(2*pi*scale^2)*((1 + (rvals/scale)^2)^(-3/2))
         },
         isPCP=TRUE,
         K = function(par,rvals, ...){
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           pi*rvals^2 + (1 - 1/sqrt(1 + rvals^2/par[2L]))/par[1L]
         },
         pcf= function(par,rvals, ...){
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           1 + ((1 + rvals^2/par[2L])^(-1.5))/(2 * pi * par[2L] * par[1L])
         },
         Dpcf= function(par,rvals, ...){
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0)){
             dkappa <- rep.int(Inf, length(rvals))
             deta2 <- rep.int(Inf, length(rvals))
           } else {
             dkappa <- -(1 + rvals^2/par[2L])^(-1.5)/(2 * pi * par[2L] * par[1L]^2)
             deta2 <- 1.5 * rvals^2 * (1 + rvals^2/par[2L])^(-2.5)/(2 * par[2L]^3 * par[1L] * pi) - (1 + rvals^2/par[2L])^(-1.5)/(2*pi*par[1L]*par[2L]^2)
           }
           out <- rbind(dkappa, deta2)
           rownames(out) <- c("kappa","eta2")
           return(out)
         },
         selfstart = function(X) {
           ## return 'par' in idiosyncratic ('old') format
           kappa <- intensity(X)
           eta2 <- 4 * mean(nndist(X))^2
           c(kappa = kappa, eta2 = eta2)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           omega <- sqrt(par[["eta2"]])/2
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, omega=omega, mu=mu)
         }
         ),
       ## ...............................................
       VarGamma=list(
         ## Neyman-Scott with VarianceGamma/Bessel clusters: old par = (kappa, eta) (internally used everywhere)
         ## Neyman-Scott with VarianceGamma/Bessel clusters: new par = (kappa, scale) (officially recommended for input/output)
         modelname = "Neyman-Scott process with Variance Gamma kernel", # In modelname field of mincon fv obj.
         descname = "Neyman-Scott process with Variance Gamma kernel", # In desc field of mincon fv obj.
         modelabbrev = "Variance Gamma process", # In fitted obj.
         printmodelname = function(obj){ # Used by print.kppm
             paste0("Variance Gamma process (nu=",
                    signif(obj$clustargs[["nu"]], 2), ")")
         },
         parnames = c("kappa", "eta"),
         clustargsnames = "nu",
         checkpar = function(par, old = TRUE, ...){
           ## 'par' is in either format
           if(is.null(par))
                 par <- c(kappa=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.", call.=FALSE)
             nam <- check.named.vector(par, c("kappa","eta"), onError="null")
             if(is.null(nam)) {
               check.named.vector(par, c("kappa","scale"))
               names(par)[2L] <- "eta"
             }
             if(!old) names(par)[2L] <- "scale"
             return(par)
         },
         checkclustargs = function(margs, old = TRUE, ...){
             if(!old)
                 margs <- list(nu=margs$nu.ker)
             return(margs)
         },
         resolvedots = function(...){
           ## resolve dots for kppm and friends allowing for old/new par syntax
           dots <- list(...)
           out <- list()
           nu <- dots$nu
           if(is.null(nu)){
             nu <- try(resolve.vargamma.shape(nu.ker=dots$nu.ker,
                                              nu.pcf=dots$nu.pcf)$nu.ker,
                       silent = TRUE)
             if(inherits(nu, "try-error"))
               nu <- -1/4
           } else {
             check.1.real(nu)
             stopifnot(nu > -1/2)
           }
           out$margs <- list(nu.ker=nu, nu.pcf=2*nu+1)
           out$covmodel <- list(type="Kernel", model="VarGamma",
                                margs=out$margs)
           return(out)
         },
         # density function for the distance to offspring
         ddist = function(r, scale, nu, ...) {
           ## 'scale' is generic format
           numer <- ((r/scale)^(nu+1)) * besselK(r/scale, nu)
           numer[r==0] <- 0
           denom <- (2^nu) * scale * gamma(nu + 1)
           numer/denom
         },
         ## Practical range of clusters
         range = function(..., par=NULL, thresh=0.001){
           ## 'par' is in generic format
           thresh <- as.numeric(thresh %orifnull% 0.001)
           scale <- retrieve.param("scale", character(0), ..., par=par)
           ## Find value of nu:
           extra <- .Spatstat.ClusterModelInfoTable$VarGamma$resolvedots(...)
           nu <- .Spatstat.ClusterModelInfoTable$VarGamma$checkclustargs(extra$margs, old=FALSE)$nu
           if(is.null(nu))
             stop(paste("Argument ", sQuote("nu"), " must be given."),
                  call.=FALSE)
           ddist <- .Spatstat.ClusterModelInfoTable$VarGamma$ddist
           f1 <- function(rmx) {
             integrate(ddist, 0, rmx, scale=scale, nu=nu)$value - (1 - thresh)
           }
           f <- Vectorize(f1)
           rmax <- uniroot(f, lower = scale, upper = 1000 * scale)$root
           return(rmax)
         },
         ## kernel function in polar coordinates (no angular argument).
         kernel = function(par, rvals, ..., margs) {
           ## 'par' is in idiosyncratic ('old') format
             scale <- as.numeric(par[2L])
             nu <- margs$nu
             if(is.null(nu))
               stop(paste("Argument ", sQuote("nu"), " is missing."),
                    call.=FALSE)
             numer <- ((rvals/scale)^nu) * besselK(rvals/scale, nu)
             numer[rvals==0] <- ifelse(nu>0, 2^(nu-1)*gamma(nu), Inf)
             denom <- pi * (2^(nu+1)) * scale^2 * gamma(nu + 1)
             numer/denom
         },
         isPCP=TRUE,
         K = local({
           ## K function requires integration of pair correlation
           xgx <- function(x, par, nu.pcf) {
             ## x * pcf(x) without check on par values
             numer <- (x/par[2L])^nu.pcf * besselK(x/par[2L], nu.pcf)
             denom <- 2^(nu.pcf+1) * pi * par[2L]^2 * par[1L] * gamma(nu.pcf + 1)
             return(x * (1 + numer/denom))
           }
           vargammaK <- function(par,rvals, ..., margs){
             ## 'par' is in idiosyncratic ('old') format
             ## margs = list(.. nu.pcf.. ) 
             if(any(par <= 0))
               return(rep.int(Inf, length(rvals)))
             nu.pcf <- margs$nu.pcf
             out <- numeric(length(rvals))
             ok <- (rvals > 0)
             rvalsok <- rvals[ok]
             outok <- numeric(sum(ok))
             for (i in 1:length(rvalsok))
               outok[i] <- 2 * pi * integrate(xgx,
                                              lower=0, upper=rvalsok[i],
                                              par=par, nu.pcf=nu.pcf)$value
             out[ok] <- outok
             return(out)
           }
           ## Initiated integration in sub-subintervals, but it is unfinished!
           ## vargammaK <- function(par,rvals, ..., margs){
           ##   ## margs = list(.. nu.pcf.. ) 
           ##   if(any(par <= 0))
           ##     return(rep.int(Inf, length(rvals)))
           ##   nu.pcf <- margs$nu.pcf
           ##   out <- numeric(length(rvals))
           ##   out[1L] <- if(rvals[1L] == 0) 0 else 
           ##   integrate(xgx, lower=0, upper=rvals[1L],
           ##             par = par, nu.pcf=nu.pcf)$value
           ##   for (i in 2:length(rvals)) {
           ##     delta <- integrate(xgx,
           ##                        lower=rvals[i-1L], upper=rvals[i],
           ##                        par=par, nu.pcf=nu.pcf)
           ##     out[i]=out[i-1L]+delta$value
           ##   }
           ##   return(out)
           ## }
           vargammaK
           }), ## end of 'local'
         pcf= function(par,rvals, ..., margs){
           ## 'par' is in idiosyncratic ('old') format
           ## margs = list(..nu.pcf..)
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           nu.pcf <- margs$nu.pcf
           sig2 <- 1 / (4 * pi * (par[2L]^2) * nu.pcf * par[1L])
           denom <- 2^(nu.pcf - 1) * gamma(nu.pcf)
           rr <- rvals / par[2L]
           ## Matern correlation function
           fr <- ifelseXB(rr > 0,
                        (rr^nu.pcf) * besselK(rr, nu.pcf) / denom,
                        1)
           return(1 + sig2 * fr)
         },
         Dpcf = NULL,
         parhandler = function(..., nu.ker = -1/4) {
           check.1.real(nu.ker)
           stopifnot(nu.ker > -1/2)
           nu.pcf <- 2 * nu.ker + 1
           return(list(type="Kernel",
                       model="VarGamma",
                       margs=list(nu.ker=nu.ker,
                                  nu.pcf=nu.pcf)))
         },
         ## sensible starting values
         selfstart = function(X) {
           ## return 'par' in idiosyncratic ('old') format
           kappa <- intensity(X)
           eta <- 2 * mean(nndist(X))
           c(kappa=kappa, eta=eta)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           omega <- par[["eta"]]
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, omega=omega, mu=mu)
         }
         ),
       ## ...............................................
       LGCP=list(
         ## Log Gaussian Cox process: old par = (sigma2, alpha) (internally used everywhere)
         ## Log Gaussian Cox process: new par = (var, scale) (officially recommended for input/output)
         modelname = "Log-Gaussian Cox process", # In modelname field of mincon fv obj.
         descname = "LGCP", # In desc field of mincon fv obj.
         modelabbrev = "log-Gaussian Cox process", # In fitted obj.
         printmodelname = function(...) "log-Gaussian Cox process", # Used by print.kppm
         parnames = c("sigma2", "alpha"),
         checkpar = function(par, old = TRUE, ...){
           ## 'par' is in either format
           if(is.null(par))
                 par <- c(var=1,scale=1)
             if(any(par<=0))
                 stop("par values must be positive.", call.=FALSE)
             nam <- check.named.vector(par, c("sigma2","alpha"), onError="null")
             if(is.null(nam)) {
                 check.named.vector(par, c("var","scale"))
                 names(par) <- c("sigma2", "alpha")
             }
             if(!old) names(par) <- c("var", "scale")
             return(par)
         },
         checkclustargs = function(margs, old = TRUE, ...) return(margs),
         resolvedots = function(...){
           ## resolve dots for kppm and friends allowing for old/new par syntax
           dots <- list(...)
           nam <- names(dots)
           out <- list()
           cmod <- dots$covmodel
           model <- cmod$model %orifnull% dots$model %orifnull% "exponential"
           margs <- NULL
           shortcut <- existsSpatstatVariable("RFshortcut") && isTRUE(getSpatstatVariable("RFshortcut"))
           if((model %in% c("exponential", "fastGauss", "fastStable", "fastGencauchy")) ||
              (shortcut && (model %in% c("gauss", "stable", "cauchy")))) {
             ## avoid RandomFields package
             ## extract shape parameters and validate them
             switch(model,
                    stable = ,
                    fastStable = {
                      stuff <- cmod %orifnull% dots
                      ok <- "alpha" %in% names(stuff)
                      if(!ok) stop("Parameter 'alpha' is required")
                      margs <- stuff["alpha"]
                      with(margs, {
                        check.1.real(alpha)
                        stopifnot(0 < alpha && alpha <= 2)
                      })
                    },
                    gencauchy = ,
                    fastGencauchy = {
                      stuff <- cmod %orifnull% dots
                      ok <- c("alpha", "beta") %in% names(stuff)
                      if(!ok[1]) stop("Parameter 'alpha' is required")
                      if(!ok[2]) stop("Parameter 'beta' is required")
                      margs <- stuff[c("alpha", "beta")]
                      with(margs, {
                        check.1.real(alpha)
                        check.1.real(beta)
                        stopifnot(0 < alpha && alpha <= 2)
                        stopifnot(beta > 0)
                      })
                    })
           } else {
             ## get the 'model generator' 
             modgen <- getRandomFieldsModelGen(model)
             attr(model, "modgen") <- modgen
             if(is.null(cmod)){
               margsnam <- names(formals(modgen))
               margsnam <- margsnam[!(margsnam %in% c("var", "scale"))]
               margs <- dots[nam %in% margsnam]
             } else{
               margs <- cmod[names(cmod)!="model"]
             }
           }
           if(length(margs)==0) {
             margs <- NULL
           } else {
             ## detect anisotropic model
             if("Aniso" %in% names(margs))
               stop("Anisotropic covariance models cannot be used",
                    call.=FALSE)
           }
           out$margs <- margs
           out$model <- model
           out$covmodel <- list(type="Covariance", model=model, margs=margs)
           return(out)
         },
         isPCP=FALSE,
         ## calls relevant covariance function from RandomFields package
         K = function(par, rvals, ..., model, margs) {
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           shortcut <- existsSpatstatVariable("RFshortcut") && isTRUE(getSpatstatVariable("RFshortcut"))
           if((model %in% c("exponential", "fastGauss", "fastStable", "fastGencauchy")) || 
              (shortcut && (model %in% c("gauss", "stable", "cauchy")))) {
             ## For efficiency and to avoid need for RandomFields package
             switch(model,
                    exponential = {
                      integrand <- function(r) { 2*pi*r*exp(par[1L]*exp(-r/par[2L])) }
                    },
                    gauss = ,
                    fastGauss = {
                      integrand <- function(r) { 2*pi*r*exp(par[1L]*exp(-(r/par[2L])^2)) }
                    },
                    stable = ,
                    fastStable = {
                      alpha <- margs[["alpha"]]
                      integrand <- function(r) { 2*pi*r*exp(par[1L]*exp(-(r/par[2L])^alpha)) }
                    },
                    gencauchy = ,
                    fastGencauchy = {
                      alpha <- margs[["alpha"]]
                      beta  <- margs[["beta"]]
                      integrand <- function(r) { 2*pi*r*exp(par[1L] * (1 + (r/par[2L])^alpha)^(-beta/alpha)) }
                    })
           } else {
             ## Use RandomFields package
             kraeverRandomFields()
             modgen <- attr(model, "modgen")
             ## create the model with the specified parameters
             if(length(margs) == 0) {
               mod <- modgen(var=par[1L], scale=par[2L])
             } else {
               mod <- do.call(modgen,
                              append(list(var=par[1L], scale=par[2L]),
                                     margs))
             }
             ## Encode integrand 
             ## RandomFields does not like evaluating covariance at r=0
             integrand <- function(r) {
               z <- numeric(length(r))
               if(any(ok <- (r != 0))) {
                 rok <- r[ok]
                 z[ok] <- 2*pi*rok*exp(RandomFields::RFcov(model=mod, x=rok))
               }
               return(z)
             }
           }
           ## compute indefinite integral
           imethod <- if(spatstat.options("fastK.lgcp")) "trapezoid" else "quadrature"
           th <- indefinteg(integrand, rvals, lower=0, method=imethod)
           return(th)
         },
         pcf= function(par, rvals, ..., model, margs) {
           ## 'par' is in idiosyncratic ('old') format
           if(any(par <= 0))
             return(rep.int(Inf, length(rvals)))
           shortcut <- existsSpatstatVariable("RFshortcut") && isTRUE(getSpatstatVariable("RFshortcut"))
           if((model %in% c("exponential", "fastGauss", "fastStable", "fastGencauchy")) ||
              (shortcut && (model %in% c("gauss", "stable", "cauchy")))) {
             ## For efficiency and to avoid need for RandomFields package
             switch(model,
                    exponential = {
                      gtheo <- exp(par[1L]*exp(-rvals/par[2L]))
                    },
                    gauss = ,
                    fastGauss = {
                      gtheo <- exp(par[1L]*exp(-(rvals/par[2L])^2))
                    },
                    stable = ,
                    fastStable = {
                      alpha <- margs[["alpha"]]
                      gtheo <- exp(par[1L]*exp(-(rvals/par[2L])^alpha))
                    },
                    gencauchy = ,
                    fastGencauchy = {
                      alpha <- margs[["alpha"]]
                      beta  <- margs[["beta"]]
                      gtheo <- exp(par[1L] * (1 + (rvals/par[2L])^alpha)^(-beta/alpha))
                    })
           } else {
             ## use RandomFields
             kraeverRandomFields()
             modgen <- attr(model, "modgen")
             ## create the model with the specified parameters
             if(length(margs) == 0) {
               mod <- modgen(var=par[1L], scale=par[2L])
             } else {
               mod <- do.call(modgen,
                              append(list(var=par[1L], scale=par[2L]),
                                     margs))
             }
             ## calculate pcf
             gtheo <- exp(RandomFields::RFcov(model=mod, x=rvals))
           }
           return(gtheo)
         },
         Dpcf= function(par,rvals, ..., model){
           ## 'par' is in idiosyncratic ('old') format
           if(!identical(model, "exponential")) {
             stop("Gradient of the pcf not available for this model.")
           } 
           dsigma2 <- exp(-rvals/par[2L]) * exp(par[1L]*exp(-rvals/par[2L]))
           dalpha <- rvals * par[1L] * exp(-rvals/par[2L]) * exp(par[1L]*exp(-rvals/par[2L]))/par[2L]^2
           out <- rbind(dsigma2, dalpha)
           rownames(out) <- c("sigma2","alpha")
           return(out)
         },
         parhandler=function(model = "exponential", ...) {
           if(!is.character(model))
             stop("Covariance function model should be specified by name",
                  call.=FALSE)
           margs <- c(...)
           if(!identical(model, "exponential")) {
             ## get the 'model generator' 
             modgen <- getRandomFieldsModelGen(model)
             attr(model, "modgen") <- modgen
           }
           return(list(type="Covariance", model=model, margs=margs))
         },
         ## sensible starting values
         selfstart = function(X) {
           ## return 'par' in idiosyncratic ('old') format
           alpha <- 2 * mean(nndist(X))
           c(sigma2=1, alpha=alpha)
         },
         ## meaningful model parameters
         interpret = function(par, lambda) {
           sigma2 <- par[["sigma2"]]
           alpha  <- par[["alpha"]]
           mu <- if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0)
             log(lambda) - sigma2/2 else NA
           c(sigma2=sigma2, alpha=alpha, mu=mu)
         }
       )
  )

spatstatClusterModelInfo <- function(name, onlyPCP = FALSE) {
  if(inherits(name, "detpointprocfamily")) {
    if(requireNamespace("spatstat.core")) {
      return(spatstat.core::spatstatDPPModelInfo(name))
    } else {
      message("The package 'spatstat.core' is required")
      return(NULL)
    }
  }
  if(!is.character(name) || length(name) != 1)
    stop("Argument must be a single character string", call.=FALSE)
  TheTable <- .Spatstat.ClusterModelInfoTable
  nama2 <- names(TheTable)
  if(onlyPCP){
    ok <- sapply(TheTable, getElement, name="isPCP")
    nama2 <- nama2[ok]
  } 
  if(!(name %in% nama2))
    stop(paste(sQuote(name), "is not recognised;",
               "valid names are", commasep(sQuote(nama2))),
         call.=FALSE)
  out <- TheTable[[name]]
  return(out)
}

## ................. helper functions (user-callable) ....................

## The following function simplifies code maintenance
## (due to changes in subscripting behaviour in recent versions of R)

retrieve.param <- function(desired, aliases, ..., par=NULL) {
  ## Retrieve the generic parameter named <desired> (or one of its <aliases>)
  ## from (...) or from 'par'
  dots <- list(...)
  par  <- as.list(par) # may be empty
  dnames <- names(dots)
  pnames <- names(par)
  for(key in c(desired, aliases)) {
    if(key %in% dnames) return(dots[[key]])
    if(key %in% pnames) return(par[[key]])
  }
  ## failed
  nali <- length(aliases)
  if(nali == 0) {
    explain <- NULL
  } else {
    explain <- paste("also tried", ngettext(nali, "alias", "aliases"), commasep(sQuote(aliases)))
  }
  mess <- paste("Unable to retrieve argument", sQuote(desired), paren(explain))
  stop(mess, call.=FALSE)
}

resolve.vargamma.shape <- function(...,
                                   nu.ker=NULL, nu.pcf=NULL, default = FALSE) {
  if(is.null(nu.ker) && is.null(nu.pcf)){
    if(!default)
        stop("Must specify either nu.ker or nu.pcf", call.=FALSE)
    nu.ker <- -1/4
  }
  if(!is.null(nu.ker) && !is.null(nu.pcf))
    stop("Only one of nu.ker and nu.pcf should be specified",
         call.=FALSE)
  if(!is.null(nu.ker)) {
    check.1.real(nu.ker)
    stopifnot(nu.ker > -1/2)
    nu.pcf <- 2 * nu.ker + 1
  } else {
    check.1.real(nu.pcf)
    stopifnot(nu.pcf > 0)
    nu.ker <- (nu.pcf - 1)/2
  }
  return(list(nu.ker=nu.ker, nu.pcf=nu.pcf))
}

