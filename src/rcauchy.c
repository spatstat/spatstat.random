/*

  rcauchy.c

  $Revision: 1.1 $ $Date: 2023/01/06 10:48:26 $

  Generate realisation of stationary Cluster cluster process in a disc D

  Baddeley-Chang hybrid algorithm

  Parameter: scale = sqrt(eta2)/2
             eta2 = 4 * scale^2

  Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022
  Licence: GNU Public Licence >= 2

 */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

/* debug activated if this is #defined */
#undef BUGGER

/* macros used */
#undef FNAME
#undef SAVEPARENTS

/* return offspring, parents and offspring-parent map */
#define FNAME rcauchyAll
#define SAVEPARENTS
#include "rcauchy.h"

#undef FNAME
#undef SAVEPARENTS

/* return offspring only */
#define FNAME rcauchyOff
#undef SAVEPARENTS
#include "rcauchy.h"


