#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void genD(void *, void *);

static const R_CMethodDef CEntries[] = {
  {"genD", (DL_FUNC) &genD, 2},
  {NULL, NULL, 0}
};

void R_init_gamlss(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
