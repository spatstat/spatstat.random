#'
#'   randomppx.R
#'
#'   $Revision: 1.2 $ $Date: 2022/04/06 07:16:18 $
#'

runifpointx <- function(n, domain, nsim=1, drop=TRUE) {
  check.1.integer(n)
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  stopifnot(inherits(domain, "boxx"))
  ra <- domain$ranges
  d <- length(ra)
  result <- vector(mode="list", length=nsim)
  for(i in seq_len(nsim)) {
    if(n == 0) {
      coo <- matrix(numeric(0), nrow=0, ncol=d)
    } else {
      coo <- mapply(runif,
                    n=rep(n, d),
                    min=ra[1,],
                    max=ra[2,])
      if(!is.matrix(coo)) coo <- matrix(coo, ncol=d)
    }
    colnames(coo) <- colnames(ra)
    df <- as.data.frame(coo)
    result[[i]] <- ppx(df, domain, coord.type=rep("s", d))
  }
  if(nsim == 1 && drop)
    return(result[[1]])
  result <- as.anylist(result)
  if(nsim > 0) names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

rpoisppx <- function(lambda, domain, nsim=1, drop=TRUE) {
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  stopifnot(inherits(domain, "boxx"))
  stopifnot(is.numeric(lambda) && length(lambda) == 1 && lambda >= 0)
  n <- rpois(nsim, lambda * volume.boxx(domain))
  result <- vector(mode="list", length=nsim)
  for(i in seq_len(nsim))
    result[[i]] <- runifpointx(n[i], domain)
  if(nsim == 1 && drop)
    return(result[[1]])
  result <- as.anylist(result)
  if(nsim > 0) names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

