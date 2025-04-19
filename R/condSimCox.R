#'
#'    condSimCox.R
#'
#'    $Revision: 1.2 $ $Date: 2025/04/19 05:30:14 $
#'
#'    Conditional simulation for Cox models

condSimCox <- function(object, nsim=1,
                       ..., win=NULL, window=win,
                       n.cond=NULL, w.cond=NULL,
                       giveup=1000, maxchunk=100,
                       saveLambda=FALSE, LambdaOnly=FALSE,
                       verbose=TRUE, drop=FALSE) {
  stopifnot(inherits(object, c("kppm", "clusterprocess", "zclustermodel")))

  w.sim <- as.owin(window)
  fullwindow <- is.null(w.cond)
  if(fullwindow) {
    w.cond <- w.sim
    w.free <- NULL
  } else {
    stopifnot(is.owin(w.cond))
    w.free <- setminus.owin(w.sim, w.cond)
  }
  
  nremaining <- nsim
  ntried <- 0
  accept <- FALSE
  nchunk <- 1
  phistory <- mhistory <- numeric(0)
  results <- list()
  while(nremaining > 0) {
    ## increase chunk length
    nchunk <- min(maxchunk, giveup - ntried, 2 * nchunk)
    ## bite off next chunk of simulations
    lamlist <- simulate(object, nsim=nchunk,
                        LambdaOnly=TRUE,
                        ..., drop=FALSE, verbose=FALSE)
    ## validate/recover
    if(!all(sapply(lamlist, is.im))) {
      if(all(sapply(lamlist, is.ppp))) {
        lamlist <- lapply(unname(lamlist), attr, which="Lambda", exact=TRUE)
        if(!all(sapply(lamlist, is.im)))
          stop("Internal error: intensity images were not saved", call.=FALSE)
      } else stop("Internal error: corrupted result from simulate method")
    }
    ## compute acceptance probabilities
    lamlist <- lapply(lamlist, "[", i=w.sim, drop=FALSE, tight=TRUE)
    if(fullwindow) {
      mu <- sapply(lamlist, integral)
    } else {
      mu <- sapply(lamlist, integral, domain=w.cond)
    }
    p <- exp(n.cond * log(mu/n.cond) + n.cond - mu)
    phistory <- c(phistory, p)
    mhistory <- c(mhistory, mu)
    ## accept/reject
    accept <- (runif(length(p)) < p)
    if(any(accept)) {
      jaccept <- which(accept)
      if(length(jaccept) > nremaining)
        jaccept <- jaccept[seq_len(nremaining)]
      naccepted <- length(jaccept)
      if(verbose)
        splat("Accepted the",
              commasep(ordinal(ntried + jaccept)),
              ngettext(naccepted, "proposal", "proposals"))
      nremaining <- nremaining - naccepted
      for(j in jaccept) {
        lamj <- lamlist[[j]]
        if(min(lamj) < 0)
          lamj <- eval.im(pmax(lamj, 0))
        if(LambdaOnly) {
          ## return only the intensity
          Y <- lamj
        } else {
          ## generate point pattern
          if(fullwindow) {
            Y <- rpoint(n.cond, lamj, win=w.sim, forcewin=TRUE)
          } else {
            lamj.cond <- lamj[w.cond, drop=FALSE, tight=TRUE]
            lamj.free <- lamj[w.free, drop=FALSE, tight=TRUE]
            Ycond <- rpoint(n.cond, lamj.cond, win=w.cond)
            Yfree <- rpoispp(lamj.free)
            Y <- superimpose(Ycond, Yfree, W=w.sim)
          }
          if(saveLambda) attr(Y, "Lambda") <- lamj
        }
        results <- append(results, list(Y))
      }
    }
    ntried <- ntried + nchunk
    if(ntried >= giveup && nremaining > 0) {
      message(paste("Gave up after", ntried,
                    "proposals with", nsim - nremaining, "accepted"))
      message(paste("Mean acceptance probability =",
                    signif(mean(phistory), 3)))
      break
    }
  }
  nresults <- length(results)
  results <- simulationresult(results, nresults, drop)
  attr(results, "history") <- data.frame(mu=mhistory, p=phistory)
  if(verbose && nresults == nsim)
    splat("Mean acceptance probability", signif(mean(phistory), 3))
  return(results)
}
