#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "plugdensity.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _typ)/sizeof(name ## _typ[0]), name ##_typ}


static R_NativePrimitiveArgType plugin_dens_typ[6] = {
    REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(plugin_dens),
    {NULL, NULL, 0}
};

void R_init_plugdensity(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
