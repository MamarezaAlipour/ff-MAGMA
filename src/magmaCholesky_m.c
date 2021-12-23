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

SEXP magmaCholeskyFinal_m(SEXP A, SEXP n, SEXP NB, SEXP zeroTri, SEXP ngpu, SEXP lowerTri)
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
	int In;
	In = INTEGER_VALUE(n);
	magma_int_t N, status, info, nGPUs;
	N = In;
	status = 0;
	nGPUs = ndevices;
	
	//perform Cholesky decomposition using MAGMA
	double *PA = NUMERIC_POINTER(A);
	int lTri, i, j;
	lTri = INTEGER_VALUE(lowerTri);
	if(lTri)
		magma_dpotrf_m(nGPUs, MagmaLower, N, PA, N, &info);
	else
		magma_dpotrf_m(nGPUs, MagmaUpper, N, PA, N, &info);
	if(info != 0)
	{
		REprintf("magma_dpotrf returned error %d: %s.\n", (int) info, magma_strerror(info));
	}
	
	//shutdown MAGMA
	magma_finalize();
	cublasShutdown();
	
	//set lower or upper triangle of output matrix to zero if necessary
	int IZeroTri;
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
