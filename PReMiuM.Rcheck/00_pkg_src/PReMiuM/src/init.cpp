#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>  // optional

#include "PReMiuM.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
   CALLDEF(profRegr, 1),
   CALLDEF(calcDisSimMat, 7),
   CALLDEF(pZpX, 10),
   CALLDEF(pYGivenZW, 14),
   CALLDEF(GradpYGivenZW, 9),
   {NULL, NULL, 0}
};

void R_init_PReMiuM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}






