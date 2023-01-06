/*

  rmatclus.h

  $Revision: 1.1 $ $Date: 2023/01/06 10:48:37 $

  Generate realisation of stationary Matern cluster process in a disc D

  Baddeley-Chang hybrid algorithm

  This file is included multiple times in rmatclus.c
  Macros:
     FNAME        name of C function
     BUGGER       print debug messages
     SAVEPARENTS  save coordinates of parents, and map from offspring to parents

  Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022
  Licence: GNU Public Licence >= 2

 */

#include <math.h>


SEXP FNAME(SEXP KAPPA,     
              SEXP MU,        
	      SEXP CLUSTERSCALE,         
	      SEXP DISCRADIUS,
              SEXP INFLATE
	      )
{
  /* generic inputs */
  double kappa, mu, scale, rD, inflate;
  
  /* generic outputs */
  double *xo, *yo;     /* offspring locations */
  SEXP Sout, Sxo, Syo;
  double *xoffspring, *yoffspring;
#ifdef SAVEPARENTS
  double *xp, *yp;     /* parent locations */
  int *ip;            /* map from offspring to parents */
  int *parentid;
  SEXP Sxp, Syp, Sip;
  double *xparent, *yparent;
#endif  
  
  /* quantities/variables used in generic algorithm */
  double lambda, kappadag, rE, rD2, rE2, Minf, MrE, diffM, p0, p0plus;
  double rpi, xpi, ypi, mi, roj, xoj, yoj, theta, muplus, dx, dy;
  int NoMax, newmax, no, i, j, n, m;
#ifdef SAVEPARENTS
  int np, added, ipcurrent;
#endif  
  
  /* model parameters (for readability) */
  double R;

  /* model-specific quantities */
  double R2, RrD, A, B;

  PROTECT(KAPPA = AS_NUMERIC(KAPPA));
  PROTECT(MU = AS_NUMERIC(MU));
  PROTECT(CLUSTERSCALE = AS_NUMERIC(CLUSTERSCALE));
  PROTECT(DISCRADIUS = AS_NUMERIC(DISCRADIUS));
  PROTECT(INFLATE = AS_NUMERIC(INFLATE));
  /* That's 5 protected */
#define NINPUTS 5
  
  GetRNGstate();
  
  /* get values */
  kappa   = *(NUMERIC_POINTER(KAPPA));
  mu      = *(NUMERIC_POINTER(MU));
  scale   = *(NUMERIC_POINTER(CLUSTERSCALE));
  rD      = *(NUMERIC_POINTER(DISCRADIUS));
  inflate = *(NUMERIC_POINTER(INFLATE));

#ifdef BUGGER
  Rprintf("INPUT: kappa = %lf, mu = %lf, scale = %lf\n", kappa, mu, scale);
  Rprintf("rD = %lf, inflate = %lf\n", rD, inflate);
#endif    

  /* model-specific translation of inputs */
  R = scale;
  
  /* specific to kernels with compact support */
  RrD      = R + rD; /* maximum distance from origin to parent, if parent has offspring in D */
  rE       = inflate * rD;
#ifdef BUGGER
  Rprintf("R + rD = %lf,\t rE = %lf\n", RrD, rE);
#endif  
  if(rE > RrD) {
    /* no need to generate parents in disc larger than RrD */
    rE = RrD;
#ifdef BUGGER
    Rprintf("Trimmed rE to %lf\n", RrD);
#endif    
  } 
  
  /* calculate some constants */
  lambda   = kappa * mu;             /* intensity of cluster process */
  kappadag = kappa * (1 - exp(-mu)); /* intensity of parents which have offspring anywhere */
  p0       = exp(-mu);               /* P(X == 0) where X ~ Pois(mu) */
  rD2      = rD * rD;
  rE2      = rE * rE;

  /* model-specific constants */
  R2       = R * R;
  muplus   = mu * rD2/R2; /* integral of dominating kernel over D (for parents in b(0, RrD)) */
  p0plus   = exp(-muplus); 
  A        = kappa * (1 - exp(- muplus)); /* intensity of dominating parents (for parents in b(0, RrD)) */
  B        = M_PI * A;
  MrE      = B * rE2;
  Minf     = B * RrD * RrD;

#ifdef BUGGER
  Rprintf("p0 = %lf, p0plus = %lf\n", p0, p0plus);
#endif    
  
  /* Guess amount of storage required */
  NoMax = (int) ceil(2.0 * M_PI * lambda * rD2);
  if(NoMax < 2048) NoMax = 2048;
  xo = (double *) R_alloc(NoMax, sizeof(double));
  yo = (double *) R_alloc(NoMax, sizeof(double));
  no = 0;
#ifdef SAVEPARENTS
  ip = (int *) R_alloc(NoMax, sizeof(int));
  xp = (double *) R_alloc(NoMax, sizeof(double));
  yp = (double *) R_alloc(NoMax, sizeof(double));
  np = 0;
#endif  

  /* -----------  parents inside E ------------------- */
  n = rpois(M_PI * rE2 * kappadag);
#ifdef BUGGER
  Rprintf("Generating %d parents inside E\n", n);
#endif
  if(n > 0) {
    for(i = 0; i < n; i++) {
      R_CheckUserInterrupt();
      /* generate parent position uniform in E */
      rpi = sqrt(runif((double) 0.0, rE2));
      theta = runif((double) 0.0, M_2PI);
      xpi = rpi * cos(theta);
      ypi = rpi * sin(theta);
#ifdef SAVEPARENTS      
      added = 0;
#endif      
      /* number of offspring of parent i: zero truncated Poisson (mu) */
      m = (int) qpois(runif(p0, (double) 1.0), mu, (int) 1, (int) 0);
#ifdef BUGGER
      Rprintf("Generating %d offspring of parent %d\n", m, i);
#endif
      /* generate offspring positions */
      for(j = 0; j < m; j++) {
	/* model specific: displacement radius */
	roj = sqrt(runif((double) 0.0, R2));
	theta = runif((double) 0.0, M_2PI);
	xoj = xpi + roj * cos(theta);
	yoj = ypi + roj * sin(theta);
	if(xoj * xoj + yoj * yoj < rD2) {
	  /* offspring point will be retained */
#ifdef SAVEPARENTS	  
	  if(added == 0) {
#ifdef BUGGER
	    Rprintf("Adding proposed parent %d to result, as parent %d\n", i, np);
#endif
	    /* add parent point  */
	    xp[np] = xpi;
	    yp[np] = ypi;
	    ipcurrent = np;
	    np++;
	    added = 1;
	  }
#endif	  
	  /* add offspring point */
#ifdef BUGGER
      Rprintf("\tAdding offspring %d to result\n", j);
#endif
	  xo[no] = xoj;
	  yo[no] = yoj;
#ifdef SAVEPARENTS	  
	  ip[no] = ipcurrent;
#endif	  
	  no++;
	  /* check data overflow */
	  if(no > NoMax) {
#ifdef BUGGER
	    Rprintf("OVERFLOW\n");
#endif
	    newmax = 2 * NoMax;
	    xo = (double *) S_realloc((char *) xo,
				      newmax, NoMax, sizeof(double));
	    yo = (double *) S_realloc((char *) yo,
				      newmax, NoMax, sizeof(double));
#ifdef SAVEPARENTS	    
	    xp = (double *) S_realloc((char *) xp,
				      newmax, NoMax, sizeof(double));
	    yp = (double *) S_realloc((char *) yp,
				      newmax, NoMax, sizeof(double));
	    ip = (int *) S_realloc((char *) ip,
				   newmax, NoMax, sizeof(int));
#endif	    
	    NoMax = newmax;
	  }
	}
      }
    }
  }

#ifdef BUGGER
#ifdef SAVEPARENTS  
  Rprintf("\n\nRunning total %d parents, %d offspring\n\n", np, no);
#else  
  Rprintf("\n\nRunning total %d offspring\n\n", no);
#endif
#endif

  /* -----------  parents outside E ------------------- */

  /* number of dominating parents */
  if(RrD <= rE || Minf <= MrE) {
    n = 0;
#ifdef BUGGER
    Rprintf("No dominating parents outside E because R+rD <= rE\n");
#endif
  } else {
    n = rpois(Minf - MrE);
#ifdef BUGGER
  Rprintf("Expect %lf dominating parents outside E\n", Minf - MrE);
  Rprintf("Generated %d dominating parents outside E\n", n);
#endif
  }

  if(n > 0) {
    for(i = 0; i < n; i++) {
      R_CheckUserInterrupt();
      /* generate parent position using dominating intensity */
      mi = runif(MrE, Minf);
      /* solve M(r) = mi for radius r */
      rpi = sqrt(mi/B);
      /* make coordinates */
      theta = runif((double) 0.0, M_2PI);
      xpi = rpi * cos(theta);
      ypi = rpi * sin(theta);
#ifdef SAVEPARENTS      
      added = 0;
#endif      
      /* number of dominating offspring */
      if(rpi > RrD) {
	/* This should not be reached */
	m = 0;
      } else {
	/* zero truncated Poisson (muplus) */
	m = (int) qpois(runif(p0plus, (double) 1.0), muplus, (int) 1, (int) 0);
      }
#ifdef BUGGER
      Rprintf("Generated %d offspring of dominating parent %d\n", m, i);
#endif
      if(m > 0) {
	for(j = 0; j < m; j++) {
	  /* generate dominating offspring uniformly in D */
	  roj = sqrt(runif((double) 0.0, rD2));
	  theta = runif((double) 0.0, M_2PI);
	  xoj = roj * cos(theta);
	  yoj = roj * sin(theta);
	  /* thin according to true kernel */
	  dx = xoj - xpi;
	  dy = yoj - ypi;
	  /* model specific */
	  /* true kernel: k(u|x) = 1/(pi R2) if |u-x| < R, 0 otherwise
	  /* dominating kernel: ktil(u|x) = 1/(pi R2) if |x| < R+rD, 0 otherwise */
	  if(dx * dx + dy * dy < R2) {
	    /* offspring will be retained */
#ifdef SAVEPARENTS	    
	    if(added == 0) {
	      /* add parent point  */
#ifdef BUGGER
	      Rprintf("Adding proposed parent %d to the output list as parent %d\n", i, np);
#endif
	      xp[np] = xpi;
	      yp[np] = ypi;
	      ipcurrent = np;
	      np++;
	      added = 1;
	    }
#endif
	    /* add offspring point */
	    xo[no] = xoj;
	    yo[no] = yoj;
#ifdef SAVEPARENTS	    
	    ip[no] = ipcurrent;
#endif	    
	    no++;
#ifdef BUGGER
	    Rprintf("\tAdded offspring %d to the output list\n", j);
#endif
	    /* check data overflow */
	    if(no > NoMax) {
#ifdef BUGGER
	    Rprintf("OVERFLOW\n");
#endif
	      newmax = 2 * NoMax;
	      xo = (double *) S_realloc((char *) xo,
					newmax, NoMax, sizeof(double));
	      yo = (double *) S_realloc((char *) yo,
					newmax, NoMax, sizeof(double));
#ifdef SAVEPARENTS	      
	      xp = (double *) S_realloc((char *) xp,
					newmax, NoMax, sizeof(double));
	      yp = (double *) S_realloc((char *) yp,
					newmax, NoMax, sizeof(double));
	      ip = (int *) S_realloc((char *) ip,
				     newmax, NoMax, sizeof(int));
#endif	      
	      NoMax = newmax;
	    }
	  }
	}
      }
    }
  }

#ifdef BUGGER
#ifdef SAVEPARENTS  
  Rprintf("Total %d parents, %d offspring\n", np, no);
#else
  Rprintf("Total %d offspring\n", no);
#endif
#endif
  
  /* copy to result */

  /* create output list */
#ifdef SAVEPARENTS
#define NOUT 5
#else
#define NOUT 2
#endif  
  PROTECT(Sout = NEW_LIST(NOUT));
  
  /* create entries in output list */
  PROTECT(Sxo = NEW_NUMERIC(no));
  PROTECT(Syo = NEW_NUMERIC(no));
#ifdef SAVEPARENTS  
  PROTECT(Sxp = NEW_NUMERIC(np));
  PROTECT(Syp = NEW_NUMERIC(np));
  PROTECT(Sip = NEW_INTEGER(no));
#endif
#define NPROTECTED (NINPUTS + 1 + NOUT)

  /* create pointers to list components */
  xoffspring = NUMERIC_POINTER(Sxo);
  yoffspring = NUMERIC_POINTER(Syo);
#ifdef SAVEPARENTS  
  xparent = NUMERIC_POINTER(Sxp);
  yparent = NUMERIC_POINTER(Syp);
  parentid = INTEGER_POINTER(Sip);
#endif

  /* copy */
#ifdef SAVEPARENTS  
  for(i = 0; i < np; i++) {
    xparent[i] = xp[i];
    yparent[i] = yp[i];
  }
#endif  
  for(j = 0; j < no; j++) {
    xoffspring[j] = xo[j];
    yoffspring[j] = yo[j];
#ifdef SAVEPARENTS    
    parentid[j] = ip[j] + 1;
#endif    
  }

  SET_VECTOR_ELT(Sout, 0, Sxo);
  SET_VECTOR_ELT(Sout, 1, Syo);
#ifdef SAVEPARENTS  
  SET_VECTOR_ELT(Sout, 2, Sxp);
  SET_VECTOR_ELT(Sout, 3, Syp);
  SET_VECTOR_ELT(Sout, 4, Sip);
#endif  

  PutRNGstate();
  
  UNPROTECT(NPROTECTED);
  return(Sout);
}

#undef NINPUTS
#undef NOUT
#undef NPROTECTED
