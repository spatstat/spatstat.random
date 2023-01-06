/*

  rtruncpois.c

  $Revision: 1.2 $ $Date: 2023/01/06 11:26:25 $

  Generate random variate with zero-truncated Poisson distribution

  Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022
  Licence: GNU Public Licence >= 2

 */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

/* ============== Functions for use in C code, to generate a single realisation =================== */

/*
   If used alone, these functions should be preceded by GetRNGState() and succeeded by PutRNGState()
*/

int rnzpoisHarding(
     double lambda  /* mean parameter */
) {
  /* from a post by Ted Harding (2005) */
  int x;
  lambda = lambda + log(runif(exp(-lambda), (double) 1.0));
  if(lambda < 0.0) return((int) 1);
  x = 1 + rpois(lambda);
  return(x);
}

int rnzpoisDalgaard(
  double lambda  /* mean parameter */
) {
  /* from a post by Peter Dalgaard (2005) in response to Harding */
  int x;
  x = qpois(runif(exp(-lambda), (double) 1.0), lambda, (int) 1, (int) 0);
  return(x);
}

int rtruncpoisHarding(
  double lambda,  /* mean parameter */
  int k           /* truncation value (minimum value of x) */
) {
  /* Adrian Baddeley, after Harding (2005) */
  int i, x;
  for(i = 0; i < k; i++) {
    lambda = lambda + log(runif(exp(-lambda), (double) 1.0));
    if(lambda < 0.0) {
      return(k);
    }
  }
  x = k + rpois(lambda);
  return(x);
}

int rtruncpoisDalgaard(
  double lambda,  /* mean parameter */
  int k           /* truncation value (minimum value of x) */
) {
  /* Adrian Baddeley, after Dalgaard (2005) */
  double pk;
  int x;
  pk = ppois(k-1, lambda, (int) 1, (int) 0);
  x = qpois(runif(pk, (double) 1.0), lambda, (int) 1, (int) 0);
  return(x);
}

/* ============  Interface to R, with efficiencies ================ */

SEXP RrnzpoisHarding(SEXP N, SEXP LAMBDA) 
{
  int n, i, nlambda;
  double *lambdavector;
  double lambda, lambdadash, expmlam;
  SEXP y;
  int *yp;

  PROTECT(N = AS_INTEGER(N));
  PROTECT(LAMBDA = AS_NUMERIC(LAMBDA));

  GetRNGstate();

  n             = *(INTEGER_POINTER(N));
  lambdavector  = NUMERIC_POINTER(LAMBDA);
  nlambda       = LENGTH(LAMBDA);

  PROTECT(y = NEW_INTEGER(n));
  yp = INTEGER_POINTER(y);
  
  if(nlambda == 1) {
    /* common value of lambda */
    lambda = lambdavector[0];
    expmlam = exp(-lambda);
    for(i = 0; i < n; i++) {
      lambdadash = lambda + log(runif(expmlam, (double) 1.0));
      yp[i] = 1 + rpois(lambdadash);
    }
  } else {
    /* vector of lambda values */
    for(i = 0; i < n; i++) {
      lambda = lambdavector[i];
      lambdadash = lambda + log(runif(exp(-lambda), (double) 1.0));      
      yp[i] = 1 + rpois(lambdadash);
    }
  }
  
  PutRNGstate();
  UNPROTECT(3);
  return(y);
}

SEXP RrnzpoisDalgaard(SEXP N, SEXP LAMBDA) 
{
  int n, i, nlambda;
  double *lambdavector;
  double lambda, lambdadash, expmlam;
  SEXP y;
  int *yp;

  PROTECT(N = AS_INTEGER(N));
  PROTECT(LAMBDA = AS_NUMERIC(LAMBDA));

  GetRNGstate();

  n             = *(INTEGER_POINTER(N));
  lambdavector  = NUMERIC_POINTER(LAMBDA);
  nlambda       = LENGTH(LAMBDA);

  PROTECT(y = NEW_INTEGER(n));
  yp = INTEGER_POINTER(y);
  
  if(nlambda == 1) {
    /* common value of lambda */
    lambda = lambdavector[0];
    expmlam = exp(-lambda);
    for(i = 0; i < n; i++) {
      yp[i] = qpois(runif(expmlam, (double) 1.0), lambda, (int) 1, (int) 0);
    }
  } else {
    /* vector of lambda values */
    for(i = 0; i < n; i++) {
      lambda = lambdavector[i];
      yp[i] = qpois(runif(exp(-lambda), (double) 1.0), lambda, (int) 1, (int) 0);
    }
  }
  
  PutRNGstate();
  UNPROTECT(3);
  return(y);
}

SEXP RrtruncpoisHarding(SEXP N, SEXP LAMBDA, SEXP TRUNC) 
{
  int n, i, k, nlambda, ntrunc;
  double *lambdavector;
  int *truncvector;
  int trunc;
  double lambda;
  SEXP y;
  int *yp;

  PROTECT(N      = AS_INTEGER(N));
  PROTECT(LAMBDA = AS_NUMERIC(LAMBDA));
  PROTECT(TRUNC  = AS_INTEGER(TRUNC));

  GetRNGstate();

  n             = *(INTEGER_POINTER(N));
  lambdavector  = NUMERIC_POINTER(LAMBDA);
  truncvector   = INTEGER_POINTER(TRUNC);
  nlambda       = LENGTH(LAMBDA);
  ntrunc        = LENGTH(TRUNC);

  PROTECT(y = NEW_INTEGER(n));
  yp = INTEGER_POINTER(y);

  lambda = lambdavector[0];
  trunc  = truncvector[0];
  
  if(nlambda == 1 && ntrunc == 1) {
    lambda = lambdavector[0];
    trunc  = truncvector[0];
    for(i = 0; i < n; i++) {
      for(k = 0; k < trunc; k++) {
	lambda = lambda + log(runif(exp(-lambda), (double) 1.0));
	if(lambda <= 0.0) {
	  yp[i] = trunc;
	  break;
	}
      }
      if(lambda > 0.0) yp[i] = trunc + rpois(lambda);
    }
  } else if(nlambda == 1 && ntrunc == n) {
    lambda = lambdavector[0];
    for(i = 0; i < n; i++) {
      trunc = truncvector[i];
      for(k = 0; k < trunc; k++) {
	lambda = lambda + log(runif(exp(-lambda), (double) 1.0));
	if(lambda <= 0.0) {
	  yp[i] = trunc;
	  break;
	}
      }
      if(lambda > 0.0) yp[i] = trunc + rpois(lambda);
    }
  } else if(nlambda == n && ntrunc == 1) {
    trunc  = truncvector[0];
    for(i = 0; i < n; i++) {
      lambda = lambdavector[i];      
      for(k = 0; k < trunc; k++) {
	lambda = lambda + log(runif(exp(-lambda), (double) 1.0));
	if(lambda <= 0.0) {
	  yp[i] = trunc;
	  break;
	}
      }
      if(lambda > 0.0) yp[i] = trunc + rpois(lambda);
    }
  } else if(nlambda == n && ntrunc == n) {
    for(i = 0; i < n; i++) {
      lambda = lambdavector[i];
      trunc = truncvector[i];
      for(k = 0; k < trunc; k++) {
	lambda = lambda + log(runif(exp(-lambda), (double) 1.0));
	if(lambda <= 0.0) {
	  yp[i] = trunc;
	  break;
	}
      }
      if(lambda > 0.0) yp[i] = trunc + rpois(lambda);
    }
  }
  
  PutRNGstate();
  UNPROTECT(4);
  return(y);
}

SEXP RrtruncpoisDalgaard(SEXP N, SEXP LAMBDA, SEXP TRUNC) 
{
  int n, i, k, nlambda, ntrunc;
  double *lambdavector;
  int *truncvector;
  int trunc;
  double lambda, ptrunc;
  SEXP y;
  int *yp;

  PROTECT(N      = AS_INTEGER(N));
  PROTECT(LAMBDA = AS_NUMERIC(LAMBDA));
  PROTECT(TRUNC  = AS_INTEGER(TRUNC));

  GetRNGstate();

  n             = *(INTEGER_POINTER(N));
  lambdavector  = NUMERIC_POINTER(LAMBDA);
  truncvector   = INTEGER_POINTER(TRUNC);
  nlambda       = LENGTH(LAMBDA);
  ntrunc        = LENGTH(TRUNC);

  PROTECT(y = NEW_INTEGER(n));
  yp = INTEGER_POINTER(y);

  if(nlambda == 1 && ntrunc == 1) {
    lambda = lambdavector[0];
    trunc  = truncvector[0];
    for(i = 0; i < n; i++) {
      ptrunc = ppois(trunc-1, lambda, (int) 1, (int) 0);
      yp[i] = qpois(runif(ptrunc, (double) 1.0), lambda, (int) 1, (int) 0);
    }
  } else if(nlambda == 1 && ntrunc == n) {
    lambda = lambdavector[0];
    for(i = 0; i < n; i++) {
      trunc = truncvector[i];
      ptrunc = ppois(trunc-1, lambda, (int) 1, (int) 0);
      yp[i] = qpois(runif(ptrunc, (double) 1.0), lambda, (int) 1, (int) 0);
    }
  } else if(nlambda == n && ntrunc == 1) {
    trunc  = truncvector[0];
    for(i = 0; i < n; i++) {
      lambda = lambdavector[i];      
      ptrunc = ppois(trunc-1, lambda, (int) 1, (int) 0);
      yp[i] = qpois(runif(ptrunc, (double) 1.0), lambda, (int) 1, (int) 0);
    }
  } else if(nlambda == n && ntrunc == n) {
    for(i = 0; i < n; i++) {
      lambda = lambdavector[i];
      trunc = truncvector[i];
      ptrunc = ppois(trunc-1, lambda, (int) 1, (int) 0);
      yp[i] = qpois(runif(ptrunc, (double) 1.0), lambda, (int) 1, (int) 0);
    }
  }
  PutRNGstate();
  UNPROTECT(4);
  return(y);
}

 
