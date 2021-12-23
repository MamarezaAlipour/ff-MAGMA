#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>

SEXP ExponentialUpperC(SEXP distMat, SEXP n, SEXP alpha, SEXP phi)
{
  int In, i, j;
  double dAlpha, dPhi;
  double *dMat, *cans;
  
  //cast R variables to C variables
  In = INTEGER(n)[0];
  dAlpha = REAL(alpha)[0];
  dPhi = REAL(phi)[0];
  dMat = REAL(distMat);
  SEXP ans = PROTECT(allocMatrix(REALSXP, In, In));
  cans = REAL(ans);
  
  //set upper triangle of output matrix
  for(i = 0; i < In; i++) {
    for(j=0; j<= i; j++) {
      if(i == j)
        cans[i*In+j] = dPhi;
      else
        cans[i*In+j] = dPhi*exp(-1*dMat[i*In+j]*dAlpha);
    }
  }
  
  UNPROTECT(1);
  return ans;
}
