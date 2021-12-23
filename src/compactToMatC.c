#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>
SEXP compactToMatC(SEXP compactMat, SEXP len, SEXP n, SEXP diagVal, SEXP lowerTri)
{
  int In, lTri, i, j, index;
  double dVal;
  double *cMat, *cans;

  // cast R variables to C variables
  In = INTEGER(n)[0];
  lTri = INTEGER(lowerTri)[0];
  dVal = REAL(diagVal)[0];
  cMat = REAL(compactMat);
  SEXP ans = PROTECT(allocMatrix(REALSXP, In, In));
  cans = REAL(ans);

  // set upper or lower triangle of output matrix
  index = 0;
  if (lTri)
  {
    for (i = 0; i < In; i++)
    {
      for (j = i + 1; j < In; j++)
      {
        cans[i * In + j] = cMat[index];
        index++;
      }
    }
  }
  else
  {
    for (i = 0; i < In; i++)
    {
      for (j = i + 1; j < In; j++)
      {
        cans[j * In + i] = cMat[index];
        index++;
      }
    }
  }

  // set diagonal values of output matrix
  for (i = 0; i < In; i++)
  {
    cans[i * In + i] = dVal;
  }

  UNPROTECT(1);
  return ans;
}
