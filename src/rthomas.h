/*

  rthomas.h

  $Revision: 1.3 $ $Date: 2023/01/25 00:59:28 $

  Generate realisation of stationary Thomas cluster process in a disc D

  Baddeley-Chang hybrid algorithm

  This file is included multiple times in rthomas.c
  Macros used:
     FNAME        name of C function
     BUGGER       activate debugging code
     SAVEPARENTS  save coordinates of parents, and map from offspring to parents

  Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022
  Licence: GNU Public Licence >= 2

 */

#define RUNIF01 runif((double) 0.0, (double) 1.0)
#define RUNIFPOS runif(DBL_EPSILON, (double) 1.0)

#define PNORM(X, MEAN, SD) pnorm(X, MEAN, SD, (int) 1, (int) 0)

#define RTRUNCPOIS(MU) (int) qpois(runif(exp(-(MU)), (double) 1.0), MU, (int) 1, (int) 0)

#define MPLUSPLUS(R) \
  lambda * M_PI * rD * rD * inv2sig2 * ( \
  rD2 + \
  2.0 * sigma2 * (1 - exp(-inv2sig2 * ((R)-rD) * ((R)-rD))) + \
  2.0 * rD * sigr2pi * ( PNORM((R)-rD, (double) 0.0, sigma) - 0.5 ))


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
  SEXP Sxp, Syp, Sip;
  double *xparent, *yparent;
  int *parentid;
#endif  

  /* quantities/variables used in generic algorithm */
  double rE, rD2, rE2, areaD;
  double rhoplus, rhoplusplus, muplus;
  double lambda, kappadag, edag, Minf, MrE, diffM, p0;
  double rpi, xpi, ypi, mi, roj, xoj, yoj, theta, dx, dy;
  int NoMax, Npmax, newmax, no, i, j, k, n, m;
  double rhi, rlo, rtry, mhi, mtry, tmp;
  double dxD, ktrue, kdom;
#ifdef SAVEPARENTS  
  int np, added, ipcurrent;
#endif  

  /* model parameters (for readability) */
  double sigma, sigma2;

  /* model-specific quantities */
  double inv2sig2, B, sigr2pi;
  
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

  /* inflation */
  rE = inflate * rD;
  
  /* model-specific translation of inputs */
  sigma = scale;

  /* calculate some constants */
  lambda   = kappa * mu;             /* intensity of cluster process */
  kappadag = kappa * (1 - exp(-mu)); /* intensity of parents which have offspring anywhere */
  p0       = exp(-mu);               /* P(X == 0) where X ~ Pois(mu) */
  rD2      = rD * rD;
  rE2      = rE * rE;
  areaD    = M_PI * rD2;

  /* model-specific constants */
  sigma2 = sigma * sigma;
  sigr2pi = sigma * sqrt(2.0 * M_PI);
  inv2sig2 = 1.0/(2.0 * sigma2);
  B = inv2sig2/M_PI;
  /* A = mu * rD2 * inv2sig2; */
  /* rDE = rE - rD; */
  
  /* superdominating intensity */
  Minf = lambda * M_PI * rD2 * inv2sig2 * (rD2 + 2.0 * sigma2 + rD * sigr2pi);
  MrE  = MPLUSPLUS(rE);
  
#ifdef BUGGER
  Rprintf("Minf = %lf, MrE = %lf\n", Minf, MrE);
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
  edag = M_PI * rE2 * kappadag;
  tmp = rpois(edag);
  n = (tmp > 2147483647.0) ? 2147483647 : ((int) tmp);
#ifdef BUGGER
  Rprintf("Expect %lf parents inside E\n", edag);
  Rprintf("Generated %d parents inside E\n", n);
#endif
  if(n > 0) {
    for(i = 0; i < n; i++) {
      R_CheckUserInterrupt();
      /* generate parent position uniform in E */
      rpi = rE * sqrt(RUNIF01);
      theta = M_2PI * RUNIF01;
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
	/* model specific: displacement radius (Box-Muller) */
	roj = sigma * sqrt( - 2.0 * log(RUNIFPOS));
	theta = M_2PI * RUNIF01;
	xoj = xpi + roj * cos(theta);
	yoj = ypi + roj * sin(theta);
	if(xoj * xoj + yoj * yoj < rD2) {
	  /* offspring point will be retained */
#ifdef SAVEPARENTS	  
	  if(added == 0) {
#ifdef BUGGER
	    Rprintf("\tAdding proposed parent %d to result, as parent %d\n", i, np);
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
      Rprintf("\t\tAdding offspring %d to result\n", j);
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
  Rprintf("\n\nRunning total %d parents, %d offspring\n\n", np, no);
#endif

  /* -----------  parents outside E ------------------- */

  diffM = Minf - MrE;
  if(diffM < 0.0) diffM = 0.0;
  
#ifdef BUGGER
  Rprintf("Expect %lf super-dominating parents outside E\n", diffM);
#endif

  /* Generate super-parents in descending order of distance */

  /* Generate values of 'mi' in descending order 
     using unit rate Poisson process */

  mi = Minf;

  /* Ensure we don't get trapped */
  Npmax = (int) ceil(diffM + 10.0 * sqrt(diffM));
#ifdef BUGGER
    Rprintf("Npmax = %d\n", Npmax);
#endif
  
  for(i = 0; i < Npmax; i++) {
    R_CheckUserInterrupt();

    mi = mi - rexp((double) 1.0);
    if(mi <= MrE) break;

#ifdef BUGGER
    Rprintf("Generated mi = %lf\n", mi);
#endif

    /* determine upper bound on solution of M(r) = mi */
    if(i == 0) {
      rhi = 2 * rE;
      for(k = 0; k < 256; k++) {
	mhi = MPLUSPLUS(rhi);
	if(mhi > mi) break;
	rhi = 2.0 * rhi;
      }
    } else {
      /* use previous value of rhi */
      mhi = MPLUSPLUS(rhi);
    }
    
    /* solve M(r) = mi */
    if(mhi <= mi) {
      /* numerical problem - failed to find upper bound */
      rpi = rhi;
#ifdef BUGGER
	Rprintf("\tFailed to find upper bound on radius\n");
#endif
    } else {
#ifdef BUGGER
      Rprintf("\tSeeking solution to Mplusplus(r) = mi on [%lf, %lf]\n",
	      rE, rhi);
#endif
      rlo = rE;
      /* mlo = MrE; */
      for(k = 0; k < 512; k++) {
	rtry = (rlo + rhi)/2.0;
	mtry = MPLUSPLUS(rtry);
	if(fabs(mtry - mi) < 0.000001) break;
	if(mtry > mi) {
	  rhi = rtry;
	} else {
	  rlo = rtry;
	}
      }
      rpi = rtry;
    }
#ifdef BUGGER      
    Rprintf("\tUsing rpi = %lf\n", rpi);
#endif

    /* compute intensities at this parent */
    dxD = rpi - rD;
    if(dxD < 0.0) dxD = 0.0;
    /* model specific */
    /* dominating kernel (for this parent, for offspring in D) */
    kdom = B * exp(-inv2sig2 * dxD * dxD);
    /* expected number of offspring using dominating kernel */
    muplus = mu * areaD * kdom;
    /* intensity of dominating parents */
    rhoplus = kappa * (1 - exp(-muplus));
    /* intensity of super-dominating parents */
    rhoplusplus = kappa * muplus;
#ifdef BUGGER
    Rprintf("\tmuplus = %lf; rhoplus = %lf; rhoplusplus = %lf\n",
	    muplus, rhoplus, rhoplusplus);
#endif
    /* THIN PARENTS TO ACHIEVE DOMINATING INTENSITY */
    if(rhoplusplus * RUNIF01 < rhoplus) {
      /* accepted */
#ifdef BUGGER
      Rprintf("\tSuper-parent %d is accepted as a dominating parent\n", i);
#endif	
      /* make coordinates */
      theta = M_2PI * RUNIF01;
      xpi = rpi * cos(theta);
      ypi = rpi * sin(theta);
      /* offspring */
#ifdef SAVEPARENTS      
      added = 0;
#endif      
      /* number of dominating offspring */
      /* zero truncated Poisson (muplus) */
      m = RTRUNCPOIS(muplus);
#ifdef BUGGER
      Rprintf("\tGenerated %d offspring of dominating parent\n", m);
#endif
      if(m > 0) {
	for(j = 0; j < m; j++) {
	  /* generate dominating offspring uniformly in D */
	  roj = rD * sqrt(RUNIF01);
	  theta = M_2PI * RUNIF01;
	  xoj = roj * cos(theta);
	  yoj = roj * sin(theta);
	  /* thin according to true kernel */
	  dx = xoj - xpi;
	  dy = yoj - ypi;
	  /* model specific */	  
	  /* true kernel: k(u|x) = 2Dgaussian(u-x) */
	  ktrue = B * exp(-inv2sig2 * (dx * dx + dy * dy));
#ifdef BUGGER
	  Rprintf("\t\tRetain offspring %d with probability %lf\n",
		  j, ktrue/kdom);
#endif
	  if(RUNIF01 * kdom < ktrue) {
	    /* offspring will be retained */
#ifdef SAVEPARENTS	    
	    if(added == 0) {
	      /* add parent point  */
#ifdef BUGGER
	      Rprintf("\t\tAdding proposed parent %d to the output list as parent %d\n", i, np);
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
	    Rprintf("\t\t\tAdded offspring %d to the output list\n", j);
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
	  } /* thinning on kernel */
	} /* loop over offspring j */
      }
    } /* thinning superparent */
  } /* loop over superparents */

#ifdef BUGGER
#ifdef SAVEPARENTS  
  Rprintf("Final total %d parents, %d offspring\n", np, no);
#else
  Rprintf("Final total %d offspring\n", no);
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

  /* create vector entries in output list */
  PROTECT(Sxo = NEW_NUMERIC(no));
  PROTECT(Syo = NEW_NUMERIC(no));
  
#ifdef SAVEPARENTS  
  PROTECT(Sxp = NEW_NUMERIC(np));
  PROTECT(Syp = NEW_NUMERIC(np));
  PROTECT(Sip = NEW_INTEGER(no));
#endif
#define NPROTECTED (NINPUTS + 1 + NOUT)

  /* create pointers to output vectors */
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

#undef RUNIF01
#undef RUNIFPOS
#undef PNORM
#undef RTRUNCPOIS
#undef MPLUSPLUS
