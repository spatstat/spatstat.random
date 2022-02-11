#'
#'   indefinteg.R
#'
#'   Indefinite integral
#'
#'   $Revision: 1.5 $ $Date: 2022/02/11 13:57:28 $

indefinteg <- function (f, x, ..., method=c("trapezoid", "quadrature"), start=min(x), nfine=8192) {
  method <- match.arg(method)
  switch(method,
         trapezoid = {
           ## indefinite integral using trapezoidal rule
           ## First make a finer sequence of x values
           ra <- range(x, start)
           xfine <- seq(ra[1L], ra[2L], length.out=nfine)
           delta <- diff(ra)/nfine
           ## Evaluate integrand on finer sequence
           yfine <- f(xfine, ...) 
           ## trapezoidal rule
           zfine <- c(0, cumsum(delta * (yfine[-1L] + yfine[-nfine]))/2)
           ## make approxfun
           Intf <- approxfun(xfine, zfine, rule=2)
           z <- Intf(x)
         },
         quadrature = {
           ## indefinite integral using 'integrate' at each value
           n <- length(x)
           z <- numeric(n)
           for(i in 1:n)
             z[i] <- integrate(f, lower=start, upper=x[i], ...)$value
         })
  return(z)
}
