/*  File src/init.c in package ergm.ego, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2015-2021 Statnet Commons
 */
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};

void R_init_ergm_ego(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
