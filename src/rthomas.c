/*

  rthomas.c

  $Revision: 1.1 $ $Date: 2023/01/06 10:48:01 $

  Generate realisation of stationary Thomas cluster process in a disc D

  Baddeley-Chang hybrid algorithm

  Copyright (C) Adrian Baddeley and Ya-Mei Chang 2022
  Licence: GNU Public Licence >= 2

 */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

/* disable debugging code */
#undef BUGGER

/* macros used */
#undef FNAME
#undef SAVEPARENTS

/* return offspring, parents and offspring-parent map */
#define FNAME rthomasAll
#define SAVEPARENTS
#include "rthomas.h"

#undef FNAME
#undef SAVEPARENTS

/* return offspring only */
#define FNAME rthomasOff
#undef SAVEPARENTS
#include "rthomas.h"


