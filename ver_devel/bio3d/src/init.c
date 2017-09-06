#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP bio3d_read_cif(SEXP, SEXP, SEXP);
extern SEXP bio3d_read_crd(SEXP);
extern SEXP bio3d_read_pdb(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bio3d_read_prmtop(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bio3d_read_cif",    (DL_FUNC) &bio3d_read_cif,    3},
    {"bio3d_read_crd",    (DL_FUNC) &bio3d_read_crd,    1},
    {"bio3d_read_pdb",    (DL_FUNC) &bio3d_read_pdb,    5},
    {"bio3d_read_prmtop", (DL_FUNC) &bio3d_read_prmtop, 1},
    {NULL, NULL, 0}
};

void R_init_bio3d(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
