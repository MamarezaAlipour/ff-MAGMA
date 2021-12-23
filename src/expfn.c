
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>

SEXP expfn(SEXP n, SEXP d2, SEXP par) {
  int In, i;
  double Dpar, par2;
  double *Pd2, *cans;
  
  //caste R variables to C variables, allocate answer vector
  In = INTEGER(n)[0];
  Dpar = REAL(par)[0];
  par2 = Dpar/2;
  Pd2 = REAL(d2);
  
  for(i = 0; i < In; i++) {
    Pd2[i] = exp(-1*pow(Pd2[i], par2));
  }
  
  return(R_NilValue);
}
