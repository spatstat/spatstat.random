#'
#'    condSimCox.R
#'
#'    $Revision: 1.8 $ $Date: 2025/05/20 08:21:36 $
#'
#'    Conditional simulation for Cox models

CondSimCox <- function(object, nsim=1,
                       ..., win=NULL, window=win,
                       n.cond=NULL, w.cond=NULL,
                       giveup=1000, maxchunk=100,
                       saveLambda=FALSE, LambdaOnly=FALSE,
                       verbose=TRUE, drop=FALSE) {
  stopifnot(inherits(object, c("kppm", "clusterprocess", "zclustermodel")))

  #' argument used only for cluster processes
  saveparents <- isTRUE(list(...)$saveparents)
  
  if(is.null(window)) {
    if(!inherits(object, "kppm"))
      stop("Argument 'window' or 'win' is required", call.=FALSE)
    ## default to window of original data
    window <- Window(object)
  }
  w.sim <- as.owin(window)
  fullwindow <- is.null(w.cond)
  if(fullwindow) {
    w.cond <- w.sim
    w.free <- NULL
  } else {
    stopifnot(is.owin(w.cond))
    w.free <- setminus.owin(w.sim, w.cond)
  }

#  if(verbose) {
#    ## calculate order-of-magnitude approximate mean acceptance probability
#    EN <- integral(predict(object, locations=w.sim), domain=w.cond)
#    phat <- exp(n.cond * log(EN/n.cond) + n.cond - EN)
#    message("Rough estimate of mean acceptance probability:", signif(phat, 2))
#  }
  
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
    lamlist <- simulate(object, window=w.sim, nsim=nchunk,
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
    ## extract cluster parents if required
    if(saveparents) {
      parentlist <- lapply(lamlist, attr, which="parents")
      if(any(sapply(parentlist, is.null))) {
        warning("Internal error: cluster parents were not returned", call.=FALSE)
        saveparents <- FALSE
      } else {
        clusker <- clusterkernel(object)
      }
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
          if(saveparents) {
            genitori <- parentlist[[j]]
            attr(Y, "parents") <- genitori
            attr(Y, "parentid") <- assignParents(Y, genitori, clusker)
            attr(Y, "cost") <- npoints(genitori) + npoints(Y)
          }
        }
        results <- append(results, list(Y))
      }
    }
    ntried <- ntried + nchunk
    naccepted <- nsim - nremaining
    if(ntried >= giveup && nremaining > 0) {
      message(paste("Gave up after", ntried,
                    "proposals with", naccepted, "accepted"))
      message(paste("Average acceptance probability =",
                    signif(mean(phistory), 3)))
      break
    }
    
    if(verbose && ((ntried >= 30 && naccepted == 0) || ((pbar <- mean(phistory)) < 1e-3)))
      message(paste("Accepted", naccepted,
                    "of the first", ntried, "proposals;",
                    "average acceptance probability =",
                    signif(mean(phistory), 3)))
  }
  nresults <- length(results)
  results <- simulationresult(results, nresults, drop)
  attr(results, "history") <- data.frame(mu=mhistory, p=phistory)
  if(verbose && nresults == nsim)
    message(paste("Average acceptance probability", signif(mean(phistory), 3)))
  return(results)
}

assignParents <- function(offspring, parents, kernfun) {
  #' assign parents to offspring in Neyman-Scott Cox model
  dx <- outer(offspring$x, parents$x, "-")
  dy <- outer(offspring$y, parents$y, "-")
  #' offspring pdf value for each (offspring, parent) pair
  kerval <- matrix(0, nrow=nrow(dx), ncol=ncol(dx))
  kerval[] <- kernfun(as.numeric(dx), as.numeric(dy))
  #' probabilities of offspring are proportional to pdf values
  kersum <- rowSums(kerval)
  probs <- kerval/kersum
  if(all(ok <- is.finite(probs) & probs >= 0 & probs <= 1)) {
    #' select at random 
    parentmap <- apply(probs, 1, function(p) { sample.int(length(p), 1, prob=p) })
  } else {
    #' handle 0/0 etc, using l'Hopital's rule
    badrow <- matrowany(!ok)
    d2 <- dx^2 + dy^2
    parentmap <- integer(nrow(probs))
    parentmap[badrow] <- apply(d2[badrow, , drop=FALSE], 1, which.min)
    if(!all(badrow))
      parentmap[!badrow] <- apply(probs[!badrow, , drop=FALSE], 1,
                                  function(p) { sample.int(length(p), 1, prob=p) })
  }
  return(parentmap)
}
