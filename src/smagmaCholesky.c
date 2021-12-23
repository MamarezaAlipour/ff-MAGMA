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

SEXP smagmaCholeskyFinal(SEXP A, SEXP n, SEXP NB, SEXP id, SEXP zeroTri, SEXP lowerTri)
{
	
	//initialize MAGMA, GPUs
	magma_init();
//	magma_print_devices();
	int In, ID;
	In = INTEGER_VALUE(n);
	ID = INTEGER_VALUE(id);
	magma_int_t N, status, info, max_size;
	N = In;
	status = 0;
	magma_setdevice(ID);
	
	//cast R objects to C data types
	double *PA = NUMERIC_POINTER(A);
	int lTri, i, j;
	lTri = INTEGER_VALUE(lowerTri);
	
	//copy input matrix into single precision matrix (only copy relevant triangle)
	//also run Cholesky decomposition on single precision matrix
	float *sPA = PROTECT(malloc(In*In * sizeof(float)));
	if(lTri) {
		for(i = 0; i < In; i++)
                {   
                        for(j = i ; j < In; j++)
                        {   
                                sPA[i*In + j] = (float) PA[i*In + j]; 
                        }   
                }
		magma_spotrf(MagmaLower, N, sPA, N, &info);
	} else {
		for(i = 0; i < In; i++)
                {   
                        for(j = 0; j <= i; j++)
                        {   
                                sPA[i*In + j] = (float) PA[i*In + j]; 
                        }   
                }
		magma_spotrf(MagmaUpper, N, sPA, N, &info);
	}
	if(info != 0)
	{
		REprintf("magma_spotrf returned error %d: %s.\n", (int) info, magma_strerror(info));
	}
	
	//shutdown MAGMA
	magma_finalize();
	cublasShutdown();
	
	//caste single precision matrix back to double and set upper or lower triangle to zero if necessary
	int IZeroTri = INTEGER_VALUE(zeroTri);
	if(!IZeroTri & lTri) {
                for(i = 0; i< In; i++) {
                        for(j=i; j < In; j++) {
                                PA[i*In + j] = (double) sPA[i*In + j]; 
                        }   
                }   
        } else if(!IZeroTri & !lTri) {
                for(i = 0; i< In; i++) {
                        for(j=0; j <= i; j++) {
                                PA[i*In + j] = (double) sPA[i*In + j]; 
                        }   
                }   
        } else if(IZeroTri & lTri) {
                for(i = 0; i< In; i++) {
                        for(j=0; j < In; j++) {
                                if(i > j)
                                         PA[i*In + j] = 0;
                                else
                                         PA[i*In + j] = (double) sPA[i*In + j]; 
                        }   
                }   
        } else if(IZeroTri & !lTri) {
                for(i = 0; i< In; i++) {
                        for(j=0; j < In; j++) {
                                if(i < j)
                                         PA[i*In + j] = 0;
                                else
                                         PA[i*In + j] = (double) sPA[i*In + j];
                        }
                }
        }
	
	//free memory
	UNPROTECT(1);
	free(sPA);
	
	return(R_NilValue);
}
	
	
		
	
