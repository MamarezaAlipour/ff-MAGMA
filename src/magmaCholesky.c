#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Print.h>
#include <cuda.h>
#include <cublas.h>
#include <magma.h>
#include <magma_lapack.h>

SEXP magmaCholeskyFinal(SEXP A, SEXP n, SEXP NB, SEXP id, SEXP zeroTri, SEXP lowerTri)
{

	//initialize MAGMA, GPU
	magma_init();
//	magma_print_devices();
	int In, ID;
	In = INTEGER_VALUE(n);
	ID = INTEGER_VALUE(id);
	magma_int_t N, status, info, max_size;
	N = In;
	status = 0;
	magma_setdevice(ID);
	
	//perform Cholesky decomposition with MAGMA
	double *PA = NUMERIC_POINTER(A);
	int lTri;
	lTri = INTEGER_VALUE(lowerTri);
	if(lTri)
		magma_dpotrf(MagmaLower, N, PA, N, &info);
	else
		magma_dpotrf(MagmaUpper, N, PA, N, &info);
	if(info != 0)
	{
		REprintf("magma_dpotrf returned error %d: %s.\n", (int) info, magma_strerror(info));
	}
	
	//shutdown MAGMA
	magma_finalize();
	cublasShutdown();
	
	//set upper or lower triangle to zero if necessary
	int IZeroTri, i, j;
        IZeroTri = INTEGER_VALUE(zeroTri);
	if(IZeroTri & lTri) {
		for(i = 1; i < In; i++)
        	{
       			for(j=0; j< i; j++)
                	{
                       		PA[i*In+j] = 0.0;
                	}
        	}
	}
	else if(IZeroTri)
		for(i = 0; i < In; i++)
                {
                        for(j=i+1; j < In; j++)
                        {
                                PA[i*In+j] = 0.0;
                        }
                }
	return(R_NilValue);
}	
