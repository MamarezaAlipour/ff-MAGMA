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
#include <errno.h>

SEXP smagmaCholeskyFinal_m(SEXP A, SEXP n, SEXP NB, SEXP zeroTri, SEXP ngpu, SEXP lowerTri)
{
	
	//initialize MAGMA, GPUs
	magma_init();
	int ndevices;
	ndevices = INTEGER_VALUE(ngpu);
        int idevice;
        for(idevice=0; idevice < ndevices; idevice++)
        {
                magma_setdevice(idevice);
                if(CUBLAS_STATUS_SUCCESS != cublasInit())
                {
                        REprintf("Error: gpu %d: cublasInit failed\n", idevice);
                        magma_finalize();
                        error("cublasInit failed");
                }
        }
//	magma_print_devices();
	
	//cast R objects to C data types
	int In, ILowerTri;
	In = INTEGER_VALUE(n);
	ILowerTri = INTEGER_VALUE(lowerTri);
	double *PA = NUMERIC_POINTER(A);
	
	//allocate single precision matrix
	float *sPA = PROTECT(malloc(In*In * sizeof(float)));
	magma_int_t N, status, info, nGPUs;
	N = In;
	status = 0;
	nGPUs = ndevices;
	
	//copy input matrix into single precision matrix (only copy relevant triangle)
	//also run Cholesky decomposition on single precision matrix
	int lTri, i, j;
	lTri = INTEGER_VALUE(lowerTri);
	if(lTri) {
		for(i = 0; i < In; i++)
                {
                        for(j = i; j < In; j++)
                        {
                                sPA[i*In + j] = (float) PA[i*In + j];
                        }
                }
		magma_spotrf_m(nGPUs, MagmaLower, N, sPA, N, &info);
	} else {
		for(i = 0; i < In; i++)
                {   
                        for(j = 0; j <= i; j++)
                        {   
                                sPA[i*In + j] = (float) PA[i*In + j]; 
                        }   
                }
		magma_spotrf_m(nGPUs, MagmaUpper, N, sPA, N, &info);
	}
	if(info != 0)
	{
		REprintf("magma_spotrf returned error %d: %s.\n", (int) info, magma_strerror(info));
	}
	
	//shutdown MAGMA, GPUs
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
	
	
		
	
