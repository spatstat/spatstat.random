#'
#'   indefinteg.R
#'
#'   Indefinite integral
#'
#'   $Revision: 1.8 $ $Date: 2022/02/12 02:56:33 $

indefinteg <- function (f, x, ...,
                        method=c("trapezoid", "quadrature"),
                        lower=min(x), nfine=8192) {
  method <- match.arg(method)
  if(length(x) == 0) return(numeric(0))
  adjust <- !missing(lower)
  if(method == "trapezoid" && (any(is.infinite(x)) ||
                               (adjust && is.infinite(lower)) ||
                               (diff(ra <- range(x)) < sqrt(.Machine$double.eps)))) {
    method <- "quadrature"
  }
  switch(method,
         trapezoid = {
           ## indefinite integral using trapezoidal rule
           ## Determine range for numerical calculation
           if(adjust) {
             check.1.real(lower)
             raplus <- ra + c(-1,1) * diff(ra)/2
             included <- inside.range(lower, raplus)
             if(included) ra <- range(ra, lower)
           }
           ## Make a fine sequence of x values
           xfine <- seq(ra[1L], ra[2L], length.out=nfine)
           delta <- diff(ra)/(nfine - 1)
           ## Evaluate integrand on finer sequence
           yfine <- f(xfine, ...) 
           ## Apply trapezoidal rule
           zfine <- c(0, cumsum(delta * (yfine[-1L] + yfine[-nfine]))/2)
           ## Evaluate at 'x'
           Intf <- approxfun(xfine, zfine, rule=2)
           z <- Intf(x)
           ## Adjust for different lower limit
           if(adjust) {
             ## calculate indefinite integral from 'lower' to min(xfine)
             x0 <- ra[1L]
             deltaI <- if(included) {
                         Intf(x0) - Intf(lower)
                       } else {
                         integrate(f, lower=lower, upper=x0, ...)$value
                       }
             ## adjust
             z <- z + deltaI
           }
         },
         quadrature = {
           ## indefinite integral using 'integrate' at each value
           n <- length(x)
           z <- numeric(n)
           for(i in 1:n)
             z[i] <- integrate(f, lower=lower, upper=x[i], ...)$value
         })
  return(z)
}
