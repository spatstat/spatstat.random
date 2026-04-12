#' Truncated Poisson random variables
#'
#'  $Revision: 1.8 $ $Date: 2026/04/12 06:23:47 $
#' 
#'   Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022 
#'   GNU Public Licence >= 2

#' --------- Poisson X given X > 0 ------------------------

dpoisnonzero <- function(x, lambda, log=FALSE) {
  ## P(X = x)
  f <- dpois(x=x, lambda=lambda, log=log)
  ## P(X > 0) 
  PG0 <- ppois(q=0, lambda=lambda, lower.tail=FALSE, log.p=log)
  ## P(X = x | X > 0)
  y <- if(log) f-PG0 else f/PG0
  y[x == 0] <- if(log) -Inf else 0
  return(y)
}

ppoisnonzero <- function(q, lambda, lower.tail=TRUE, log.p=FALSE) {
  if(lower.tail) {
    ## P(X <= q)
    p  <- ppois(q=q, lambda=lambda, lower.tail=TRUE,  log.p=FALSE)
    ## P(X = 0)
    p0 <- dpois(x=0, lambda=lambda, log=FALSE)
    ## P(X <= q | X > 0)
    y <- (p - p0)/(1 - p0)
    if(log.p) y <- log(y)
  } else {
    ## P(X > q)
    p   <- ppois(q=q, lambda=lambda, lower.tail=FALSE, log.p=log.p)
    ## P(X > 0)
    PG0 <- ppois(q=0, lambda=lambda, lower.tail=FALSE, log.p=log.p)
    ## P(X > q | X > 0)
    y <- if(log.p) p - PG0 else p/PG0
  }
  return(y)
}

qpoisnonzero <- function(p, lambda, lower.tail=TRUE, log.p=FALSE) {
  if(lower.tail) {
    ## find x such that P(X <= x | X > 0) = p 
    ## P(X = 0)
    p0 <- dpois(x=0, lambda=lambda, log=FALSE)
    ## P(X <= x) = P(X = 0) + P(X > 0) P(X <= x | X > 0)
    pun <- p0 + (1-p0) * (if(log.p) exp(p) else p)
    ## quantiles of Poisson
    y <- qpois(p=pun, lambda, lower.tail=TRUE, log.p=FALSE)
  } else {
    ## find x such that P(X > x | X > 0) = p 
    ## P(X > 0)
    PG0 <- ppois(q=0, lambda=lambda, lower.tail=FALSE, log.p=log.p)
    ## P(X > x) = P(X > 0) P(X > x | X > 0)
    pun <- if(log.p) PG0 + p else PG0 * p
    ## quantiles of Poisson
    y <- qpois(p=pun, lambda, lower.tail=FALSE, log.p=log.p)
  }
  y[] <- pmax(y[], 1)
  return(y)
}

rpoisnonzero <- function(n, lambda, method=c("harding", "transform"), implem=c("R", "C")) {
  ## Poisson random variable, conditioned to be nonzero
  method <- match.arg(method)
  implem <- match.arg(implem)
  switch(implem,
         R = {
           switch(method,
                  harding = {
                    ## From a post by Ted Harding (2005)
                    lam1 <- lambda + log(runif(n, min=exp(-lambda), max=1))
                    lam1 <- pmax(0, lam1) ## avoid numerical glitches (lam1 is theoretically > 0)
                    y <- rpois(n, lam1) + 1L
                  },
                  transform = {
                    ## From a post by Peter Dalgaard (2005) in response to Harding
                    ## Surprisingly, this is 3 times slower!
                    y <- as.integer(qpois(runif(n, min=exp(-lambda), max=1), lambda))
                  })
         },
         C = {
           storage.mode(n) <- "integer"
           storage.mode(lambda) <- "double"
           switch(method,
                  harding = {
                    y <- .Call(SR_RrnzpoisHarding,
                               n, lambda,
                               PACKAGE="spatstat.random")
                  },
                  transform = {
                    y <- .Call(SR_RrnzpoisDalgaard,
                               n, lambda,
                               PACKAGE="spatstat.random")
                  })
         })
  return(y)
}

#' --------- Poisson X given X >= minimum  ------------------------

dpoistrunc <- function(x, lambda, minimum=1, log=FALSE) {
  minimum <- as.integer(minimum)
  ## P(X = x)
  f <- dpois(x=x, lambda=lambda, log=log)
  ## P(X >= m) 
  PG <- ppois(q = minimum - 1, lambda=lambda, lower.tail=FALSE, log.p=log)
  ## P(X = x | X >= m)
  y <- if(log) f-PG else f/PG
  y[x < minimum] <- if(log) -Inf else 0
  return(y)
}

ppoistrunc <- function(q, lambda, minimum=1, lower.tail=TRUE, log.p=FALSE) {
  if(lower.tail) {
    ## P(X <= q)
    p  <- ppois(q=q,         lambda=lambda, lower.tail=TRUE, log.p=FALSE)
    ## P(X < m)
    PL <- ppois(q=minimum-1, lambda=lambda, lower.tail=TRUE, log.p=FALSE)
    ## P(X <= q | X >= m)
    y <- (p - PL)/(1 - PL)
    y[q < minimum] <- 0
    if(log.p) y <- log(y)
  } else {
    ## P(X > q)
    p   <- ppois(q=q, lambda=lambda, lower.tail=FALSE, log.p=log.p)
    ## P(X >= m)
    PG <- ppois(q = minimum - 1, lambda=lambda, lower.tail=FALSE, log.p=log)
    ## P(X > q | X >= m)
    y <- if(log.p) p - PG else p/PG
  }
  return(y)
}

qpoistrunc <- function(p, lambda, minimum=1, lower.tail=TRUE, log.p=FALSE) {
  if(lower.tail) {
    ## find x such that P(X <= x | X >= m) = p 
    ## P(X < m)
    PL <- ppois(q = minimum - 1, lambda=lambda, lower.tail=TRUE, log.p=log.p)
    ## P(X <= x) = P(X < m) + P(X >= m) P(X <= x | X >= m)
    pun <- PL + (1-PL) * (if(log.p) exp(p) else p)
    ## quantiles of Poisson
    y <- qpois(p=pun, lambda, lower.tail=TRUE, log.p=FALSE)
  } else {
    ## find x such that P(X > x | X >= m) = p 
    ## P(X >= m)
    PG <- ppois(q = minimum - 1, lambda=lambda, lower.tail=FALSE, log.p=log.p)
    ## P(X > x) = P(X >= m) P(X > x | X >= m)
    pun <- if(log.p) PG + p else PG * p
    ## quantiles of Poisson
    y <- qpois(p=pun, lambda, lower.tail=FALSE, log.p=log.p)
  }
  y[] <- pmax(y[], minimum)
  return(y)
}

rpoistrunc <- function(n, lambda, minimum=1, method=c("harding", "transform"), implem=c("R", "C")) {
  ## Poisson random variable, conditioned to be at least 'minimum'
  stopifnot(all(is.finite(minimum)))
  minimum <- pmax(0L, as.integer(minimum))
  method <- match.arg(method)
  implem <- match.arg(implem)
  switch(implem,
         R = {
           switch(method,
                  transform = {
                    y <- qpois(runif(n, min=ppois(minimum-1L, lambda), max=1), lambda)
                  },
                  harding = {
                    if(length(minimum) == 1) {
                      for(k in seq_len(minimum)) 
                        lambda <- pmax(0, lambda + log(1 - runif(n) * (1 - exp(-lambda))))
                    } else if(length(minimum) == n) {
                      if(length(lambda) == 1) lambda <- rep(lambda, n)
                      remaining <- minimum
                      while(any(todo <- (remaining > 0))) {
                        lambda[todo] <- pmax(0, lambda[todo] + log(1 - runif(sum(todo)) * (1 - exp(-lambda[todo]))))
                        remaining[todo] <- remaining[todo] - 1L
                      }
                    } else stop("Argument 'minimum' should be a vector of length 1 or n", call.=FALSE)
                    y <- rpois(n, lambda) + minimum
                  })
         },
         C = {
           storage.mode(n) <- "integer"
           storage.mode(lambda) <- "double"
           storage.mode(minimum) <- "integer"
           switch(method,
                  harding = {
                    y <- .Call(SR_RrtruncpoisHarding,
                               n, lambda, minimum, 
                               PACKAGE="spatstat.random")
                  },
                  transform = {
                    y <- .Call(SR_RrtruncpoisDalgaard,
                               n, lambda, minimum, 
                               PACKAGE="spatstat.random")
                  })
         })
  return(y)
}

recipEnzpois <- function(mu, exact=TRUE) {
  ## first reciprocal moment of nzpois
  if(exact && isNamespaceLoaded("gsl")) {
    gamma <- -digamma(1)
    ans <- (gsl::expint_Ei(mu) - log(mu) - gamma) * exp(-mu) /(1 - exp(-mu))
    return(ans)
  } else {
    n <- length(mu)
    ans <- numeric(n)
    xx <- 1:max(ceiling(mu + 6 * sqrt(mu)), 100)
    for(i in 1:n) 
      ans[i] <- sum(dpois(xx, mu[i])/xx)/(1 - exp(-mu[i]))
  }
  return(ans)
}



