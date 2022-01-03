#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include "methas.h"

#include "mhsnoopdef.h"

/*
  mhsnoop.c

  $Revision: 1.11 $  $Date: 2021/12/24 04:28:15 $

  support for visual debugger in RMH

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

/* 
   To switch on debugging code, 
   insert the line: #define MH_DEBUG YES

   To switch off debugging code, 
   insert the line: #define MH_DEBUG NO
*/
#ifndef MH_DEBUG
#define MH_DEBUG NO
#endif


void initmhsnoop(Snoop *s, SEXP env) {
  s->active = isEnvironment(env);
  s->nextstop = 0;         /* stop at iteration 0 */
  s->nexttype = NO_TYPE;   /* deactivated */
  if(s->active) {
    s->env = env;
    s->expr = findVar(install("callbackexpr"), env);
  } else {
    s->env = s->expr = R_NilValue;
  }
}

void mhsnoop(Snoop *s, 
	     int irep, 
	     Algor *algo,
	     State *state,
	     Propo *prop,
	     double numer,
	     double denom,
	     int *itype) 
{
  SEXP e;
  int npts, j, proptype, accepted, fateMH, fateUser;
  /* passed from C to R before debugger */
  SEXP Sirep, Sx, Sy, Sm, Sproptype, Sproplocn, Spropmark, Spropindx;
  SEXP Snumer, Sdenom, Sitype;
  double *Px, *Py, *Pproplocn;
  int *Pm;
  /* passed from R to C after debugger */
  SEXP Sinxt, Stnxt, SitypeUser;

#if MH_DEBUG
  Rprintf("mhsnoop called at iteration %d\n", irep);
#endif

  if(!(s->active)) return;

#if MH_DEBUG
  Rprintf("mhsnoop is active\n");
#endif

  /* 
     execute when the simulation reaches the next stopping time
 */

  if(s->nextstop >= 0) {
    /*  specified iteration number 'nextstop' or later */
    if(irep < s->nextstop) return;
  } else if(s->nexttype >= 0) {
    /*  specified proposal type 'nexttype' */
    if(prop->itype != s->nexttype) return;
  } else {
    /* no stopping rule - skip all */
    return;
  }

#if MH_DEBUG
  Rprintf("debug triggered\n");
#endif

  /* environment for communication with R */
  e = s->env;
  /* 
     copy data to R
  */
  /* copy iteration number */
  PROTECT(Sirep = NEW_INTEGER(1));
  *(INTEGER_POINTER(Sirep)) = irep;
  setVar(install("irep"), Sirep, e);
  UNPROTECT(1);
  /* copy (x,y) coordinates */
  npts = state->npts;
  PROTECT(Sx = NEW_NUMERIC(npts));
  PROTECT(Sy = NEW_NUMERIC(npts));
  Px = NUMERIC_POINTER(Sx);
  Py = NUMERIC_POINTER(Sy);
  for(j = 0; j < npts; j++) {
    Px[j] = state->x[j];
    Py[j] = state->y[j];
  }
  setVar(install("xcoords"), Sx, e);
  setVar(install("ycoords"), Sy, e);
  UNPROTECT(2);
  /* copy marks */
  if(state->ismarked) {
    PROTECT(Sm = NEW_INTEGER(npts));
    Pm = INTEGER_POINTER(Sm);
    for(j = 0; j < npts; j++) {
      Pm[j] = state->marks[j];
    }
    setVar(install("mcodes"), Sm, e);
    UNPROTECT(1);
  }
  /* proposal type */
  PROTECT(Sproptype = NEW_INTEGER(1));
  *(INTEGER_POINTER(Sproptype)) = proptype = prop->itype;
  setVar(install("proptype"), Sproptype, e);
  UNPROTECT(1);
  /* proposal coordinates */
  PROTECT(Sproplocn = NEW_NUMERIC(2));
  Pproplocn = NUMERIC_POINTER(Sproplocn);
  Pproplocn[0] = prop->u;
  Pproplocn[1] = prop->v;
  setVar(install("proplocn"), Sproplocn, e);
  UNPROTECT(1);
  /* proposal mark value */
  if(state->ismarked) {
    PROTECT(Spropmark = NEW_INTEGER(1));
    *(INTEGER_POINTER(Spropmark)) = prop->mrk;
    setVar(install("propmark"), Spropmark, e);
    UNPROTECT(1);
  }
  /* proposal point index */
  PROTECT(Spropindx = NEW_INTEGER(1));
  *(INTEGER_POINTER(Spropindx)) = prop->ix;
  setVar(install("propindx"), Spropindx, e);
  UNPROTECT(1);
  /* Metropolis-Hastings numerator and denominator */
  PROTECT(Snumer = NEW_NUMERIC(1));
  PROTECT(Sdenom = NEW_NUMERIC(1));
  *(NUMERIC_POINTER(Snumer)) = numer;
  *(NUMERIC_POINTER(Sdenom)) = denom;
  setVar(install("numerator"), Snumer, e);
  setVar(install("denominator"), Sdenom, e);
  UNPROTECT(2);
  /* tentative outcome of proposal (0 = reject, other=accept) */
  PROTECT(Sitype = NEW_INTEGER(1));
  *(INTEGER_POINTER(Sitype)) = fateMH = *itype;
  setVar(install("itype"), Sitype, e);
  UNPROTECT(1);

  /* ..... call visual debugger */

#if MH_DEBUG
  Rprintf("executing [callback]\n");
#endif

  eval(s->expr, s->env);
  
#if MH_DEBUG
  Rprintf("exited [callback]\n");
#endif

  /* update outcome of proposal */
  SitypeUser = findVar(install("itype"), e);
  fateUser = *(INTEGER_POINTER(SitypeUser));
  if(fateUser != fateMH)
    *itype = fateUser;

#if MH_DEBUG
  Rprintf("Returned itype = %d\n", fateUser);
  if(fateUser == fateMH) {
    if(fateMH == REJECT) {
      Rprintf("Confirmed: Proposal rejected\n");
    } else {
      Rprintf("Confirmed: Proposal accepted\n");
    }
  } else {
    if(fateUser == REJECT) {
      Rprintf("User changed fate of proposal to REJECTED\n");
    } else {
      Rprintf("User changed fate of proposal to ACCEPTED\n");
    }
  }
  Rprintf("Assigned itype = %d\n", *itype);
#endif
    
  /* update stopping time */
  Sinxt = findVar(install("inxt"), e);
  s->nextstop = *(INTEGER_POINTER(Sinxt));
  Stnxt = findVar(install("tnxt"), e);
  s->nexttype = *(INTEGER_POINTER(Stnxt));

#if MH_DEBUG
  if(s->nextstop >= 0)
    Rprintf("Next stop: iteration %d\n", s->nextstop);
  if(s->nexttype >= 0) {
    if(s->nexttype == BIRTH) 
      Rprintf("Next stop: first birth proposal\n");
    if(s->nexttype == DEATH) 
      Rprintf("Next stop: first death proposal\n");
    if(s->nexttype == SHIFT) 
      Rprintf("Next stop: first shift proposal\n");
  }
  Rprintf("Exiting mhsnoop\n");
#endif

  return;
}
	     
