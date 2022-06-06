#'
#'   randompp3.R
#'
#'   $Revision: 1.2 $ $Date: 2022/05/23 02:33:06 $
#'

runifpoint3 <- function(n, domain=box3(), nsim=1, drop=TRUE) {
  domain <- as.box3(domain)
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  result <- vector(mode="list", length=nsim)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for(i in seq_len(nsim)) {
    x <- with(dd, runif(n, min=xrange[1], max=xrange[2]))
    y <- with(dd, runif(n, min=yrange[1], max=yrange[2]))
    z <- with(dd, runif(n, min=zrange[1], max=zrange[2]))
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  if(nsim > 0) names(result) <- paste("Simulation", seq_len(nsim))
  return(result)
}

rpoispp3 <- function(lambda, domain=box3(), nsim=1, drop=TRUE) {
  domain <- as.box3(domain)
  v <- volume(domain)
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  if(!(is.numeric(lambda) && length(lambda) == 1))
    stop("lambda must be a single numeric value")
  np <- rpois(nsim, lambda * v)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  result <- vector(mode="list", length=nsim)
  for(i in seq_len(nsim)) {
    ni <- np[i]
    x <- with(dd, runif(ni, min=xrange[1], max=xrange[2]))
    y <- with(dd, runif(ni, min=yrange[1], max=yrange[2]))
    z <- with(dd, runif(ni, min=zrange[1], max=zrange[2]))
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  if(nsim > 0) names(result) <- paste("Simulation", seq_len(nsim))
  return(result)
}

