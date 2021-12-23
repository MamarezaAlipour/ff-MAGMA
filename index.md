# ff-MAGMA
GPU functionality to the chol function using MAGMA.
A MAGMA-Accelerated finite field Extension to the `fields` Package.

# Overview
## Description
ff-MAGMA is an extension of the fields package, which is freely available on CRAN 
  and is used for spatial statistics and analyzing spatial data. The mKrig and 
  mKrig.MLE functions from fields are reorganized to be able to use newly provided
  covariance functions that avoid recomputing the distance matrix and only compute 
  the upper triangle for symmetric covariance matrices when possible. This package
  also supports the use of the freely available and open source MAGMA library for 
  using GPUs in addition to multi-core CPUs. MAGMA must be installed on your computer 
  and it must be linked along with its prerequesite libraries in your ~/.R/Makevars 
  file.
  MAGMA requires CUDA and LAPACK, which are also free.
  
## Hardware Requirements
Although it may be possible to install ff-MAGMA with older versions of CUDA and 
MAGMA, these instructions assume CUDA version 6.5 or above and MAGMA 1.6.1 or above.  
It is also assumed that the user has an NVIDIA graphics card with CUDA compute 
capability 2.0 or above.  The gcc or icc compilers are recommended for building 
MAGMA.  The clang compiler (the default compiler on Macs) sometimes has problems 
when building MAGMA.  This is because most toolchains for R using clang, which is 
the default compiler used on Macs, do not support openMP.  openMP is a tool used 
for creating parallel code and we do not currently know of a way of building MAGMA 
without openMP.  It is therefore highly recommended if not required to use an 
alternative compiler to clang.

## Supported Operating Systems
I have only been able to install MAGMA on unix-based systems so far.  It is 
possible, however to install MAGMA on a Windows system.  In that case, it may 
also be possible to install this package on that system, but these instructions 
do not detail how to do so.

# Installation Instructions
Detailed installation instructions are given for specific Linux and Mac OS X 
systems in the following technical reports that we plan to submit for publication 
soon:

Parisa Khaleghi, Naser Rahmati, Reza Hayati and Sara Hedayat (2021), "Incorporating MAGMA into the `fields' spatial statistics 
  package" Technical report, National Center for Amirkabir University Research.

Parisa Khaleghi, Naser Rahmati and Reza Hayati (2021), "Using single-precision, 
  MAGMA-accelerated Cholesky decompositions in the `fields' spatial statistics 
  package" Technical report, National Center for Amirkabir University Research.

The following installation instructions will be more general, but it is 
recommended to use the above technical reports as reference.  Until they are 
published, recent versions of the above technical reports will be available on 
github.  The github repository can be cloned with the command:

``` git clone https://github.com/parikhaleghi/ff-MAGMA.git ```

## a) gcc/icc Compiler Suite
Note that most toolchains for R using clang, which is the default compiler used 
on Macs, do not support openMP.  openMP is a tool used for creating parallel 
code and we do not currently know of a way of building MAGMA without openMP.  It 
is therefore highly recommended if not required to use an alternative compiler to 
clang.

The gcc compiler suite is available at https://gcc.gnu.org/.  On a Mac it is very 
easy to install using Macports, which is available at https://www.macports.org/.  
Our installation instructions assume gcc 4.8 and above with gfortran, although 
other versions will likely work as well.

The icc Intel compilers are not free, but can also be used to install MAGMA.  
They are available at https://software.intel.com/en-us/intel-compilers.

## b) NVIDIA CUDA Toolkit
Instructions for installing the NVIDIA CUDA Toolkit on Linux can be found at 
http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-linux/#axzz3gq5mX46P
and similar instructions for installation on Mac OS X are at
https://www.clear.rice.edu/comp422/resources/cuda/html/cuda-getting-started-
guide-for-mac-os-x/index.html

## c) LAPACK Shared Library
Many LAPACK and BLAS libraries can be used when installing MAGMA.  These 
instructions assume either OpenBLAS or Intel's MKL. OpenBLAS is freely 
available at http://www.openblas.net/ (or via Macports on Mac OS X).  MKL is not 
free, but can be purchased at https://software.intel.com/en-us/intel-mkl.  
The installation process will likely be similar when installing MAGMA with 
other libraries (e.g. ATLAS, GotoBLAS), and MAGMA provides their own installation 
supplementary instructions.

## MAGMA Software
Download MAGMA 1.6.1 (or possibly another version but we cannot guarantee the 
success of the installation) from http://icl.cs.utk.edu/magma/software/index.html 
and unpack the .tar file in your home directory.  MAGMA provides some 
installation instructions in the unpacked directory's README file if you are 
using different libraries than those assumed in these instructions.

MAGMA assumes that certain variables have been set in the .bashrc file in your 
home directory (assuming you are using a bash shell).  The first technical 
report given above gives instructions if you are using tcsh instead of bash.  
The following lines should be included in your .bashrc file (with "XXXX" and 
other paths substituted for the appropriate path for your specific system and 
installation):
### Linux with icc 15.0.3, CUDA 6.5, and MKL 11.2.3:
```
export MKLROOT=XXXX/intel/psxe-2015_update3/composer_xe_2015.3.187/mkl/
export CUDADIR=XXXX/cuda/6.5
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:XXXX/intel/psxe-2015_update3/\
composer_xe_2015.3.187/mkl/lib/intel64:XXXX/cuda/6.5/lib64
export PATH=${PATH}:XXXX/cuda/6.5/bin
```
### Mac OS X with gcc 4.8, CUDA 7.0, and OpenBLAS:
```
export PATH=/Developer/NVIDIA/CUDA-7.0/bin:$PATH
export DYLD_LIBRARY_PATH=/usr/local/cuda/lib/:/Developer/NVIDIA/\
CUDA-7.0/lib:$DYLD_LIBRARY_PATH
export CUDADIR=/Developer/NVIDIA/CUDA-7.0
export OPENBLASDIR=/opt/local

source ~/.profile
alias gfortran-4.8='gfortran'
```
Note that in the Mac OS X .bashrc file shown above, the last two lines 
are only required if Macports is used to download gcc or OpenBLAS, and 
the second-to-last line is only required if Macports is used to download 
OpenBLAS.

The following are examples of make.inc files for building MAGMA that are 
used in the technical reports given above for specific Linux and Mac OS X 
systems.  The make.inc files should be put in MAGMA's outermost directory 
created from the untar-ed .tar MAGMA file.  After the make.inc file is 
put there, MAGMA can be installed on your system.  MAGMA has several 
example make files for different BLAS and LAPACK libraries, different 
compilers, and different operating systems.  These example make files 
as well as MAGMA's README file contain instructions for modifications 
that need to be made depending on the specific system on which you build 
MAGMA.

Here are the make.inc files we used to build MAGMA on the specific 
systems used in the above technical reports:
### Linux with icc 15.0.3, CUDA 6.5, and MKL 11.2.3:
```
GPU_TARGET = Kepler

CC        = icc
CXX       = icpc
NVCC      = nvcc
FORT      = ifort

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

# Use -fPIC to make shared (.so) and static (.a) library;
# can be commented out if making only static library.
FPIC      = -fPIC

#CFLAGS    = -O3 $(FPIC) -DADD_ -Wall     -fopenmp -DMAGMA_SETAFFINITY -DMAGMA_WITH_MKL
CFLAGS    = -O3 $(FPIC) -DADD_ -Wall     -fopenmp -DMAGMA_SETAFFINITY -DMAGMA_WITH_MKL -Xlinker -shared
FFLAGS    = -O3 $(FPIC) -DADD_ -warn all -warn nounused
F90FLAGS  = -O3 $(FPIC) -DADD_ -warn all -warn nounused
NVCCFLAGS = -O3         -DADD_           -Xcompiler "-fno-strict-aliasing $(FPIC)"
LDFLAGS   =     $(FPIC)                  -fopenmp

# old MKL
#LIB       = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_lapack -lmkl_core -lguide -lpthread -lcublas -lcudart -lstdc++ -lm

# see MKL Link Advisor at http://software.intel.com/sites/products/mkl/
# icc with MKL 10.3, Intel threads
LIB       = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lcublas -lcudart -lstdc++ -lm

# define library directories preferably in your environment, or here.
# for MKL run, e.g.: source /opt/intel/composerxe/mkl/bin/mklvars.sh intel64
#MKLROOT ?= /opt/intel/composerxe/mkl
#CUDADIR ?= /usr/local/cuda
-include make.check-mkl
-include make.check-cuda

LIBDIR    = -L$(CUDADIR)/lib64 \
            -L$(MKLROOT)/lib/intel64

INC       = -I$(CUDADIR)/include \
            -I$(MKLROOT)/include
```

### Mac OS X with gcc 4.8, CUDA 7.0, and OpenBLAS:
```
GPU_TARGET = Kepler

CC        = gcc
CXX       = g++
NVCC      = nvcc
FORT      = gfortran

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

# use -m32 to compile with 32-bit long & pointers.
# use -m64 to compile with 64-bit long & pointers (lp64). int is still 32-bit.
# add -DNDEBUG to disable asserts and certain error checks.
#
# MacOS veclib has a bug where some single precision functions return
# a double precision result, for instance slange.
# This is observed with -m64, but oddly not with -m32.
# The easiest fix is to replace those routines with correct ones from LAPACK.
# See BLAS_FIX below.
# Alternatively, don't link with the veclib/accelerate framework;
# use a different BLAS and LAPACK library.

# Use -fPIC to make shared (.so) and static (.a) library;
# can be commented out if making only static library.
FPIC      = -fPIC

CFLAGS    = -m64 -O3 $(FPIC) -DADD_ -Wall -fopenmp
FFLAGS    = -m64 -O3 $(FPIC) -DADD_ -Wall -Wno-unused-dummy-argument
F90FLAGS  = -m64 -O3 $(FPIC) -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input
NVCCFLAGS = -m64 -O3         -DADD_       -Xcompiler "-fno-strict-aliasing $(FPIC)"
LDFLAGS   = -m64     $(FPIC) -fopenmp

# MacOS likes the library's path to be set
INSTALL_NAME = -install_name @rpath/

LIB       = -framework Accelerate -lopenblas -lcublas -lcudart -lstdc++ -lm

# define library directories preferably in your environment, or here.
#OPENBLASDIR ?= /usr/local/openblas
#CUDADIR ?= /usr/local/cuda
-include make.check-openblas
-include make.check-cuda

LIBDIR    = -L$(CUDADIR)/lib \
            -L$(OPENBLASDIR)/lib

INC       = -I$(CUDADIR)/include \
            -I$(OPENBLASDIR)/include


# ========================================
# replace single & single-complex BLAS functions with reference versions.
# (i.e., functions that return float; subroutines do not need a fix.)
LIB      := -L$(MAGMA_DIR)/lib -lblas_fix $(LIB)

BLAS_FIX  = $(MAGMA_DIR)/lib/libblas_fix.a

.PHONY: blas_fix

blas_fix:
  @echo "======================================== BLAS fix for MacOS"
	( cd $(MAGMA_DIR)/blas_fix && $(MAKE) )
	@echo

lib: blas_fix
```

After make.inc has been created in MAGMA's directory, run "make" to build 
MAGMA (or "make clean" and then "make" if this is not your first time 
building MAGMA from this directory).

If MAGMA was built without errors, run the following on the command line 
in MAGMA's testing directory to make sure MAGMA was installed successfully:
```
./testing_dpotrf -c --ngpu 1
./testing_dpotrf_m -c --ngpu 1
./testing_spotrf -c --ngpu 1
./testing_spotrf_m -c --ngpu 1
```
Note that the -h option lets the user see specific options that can be 
used with the above testing scripts to get other options the testing 
scripts can be used with.  If they all work correctly, then MAGMA has been 
successfully installed on your system.

## ff-MAGMA
After you have made sure MAGMA works correctly on your system, you can install 
ff-MAGMA by modifying (and possibly creating) a ~/.R/Makevars file using 
your command line.  To do this, open your command line and run:
```
cd ~
ls -a
```
If a ".R" directory is shown then run:
```
cd .R
ls
```
If the "Makevars" file is not shown, run ``` touch Makevars ``` to create it.  The Makevars file contains variables that R uses to build 
packages.  Here are the variables important for this installation along with 
descriptions of what they are for:

``CC``: C compiler

``CXX``: C++ compiler

``F77``: Fortran77 compiler

``FC``: Fortran compiler

``PKG_CFLAGS``: Flags (or options) for the C compiler used when building packages.  
    Also which C headers to include when compiling C code.
    
``PKG_LIBS``: Libraries to include when building packages.  In our case, we need 
    MAGMA, CUDA, and any BLAS and LAPACK libraries you used.

Now make sure the following lines are in your Makevars file if 
you are using icc (the Intel C compiler) and built MAGMA with the MKL BLAS 
library on a 64 bit Linux system with CUDA 6.5:

```
CC=icc
CXX=icc
F77=ifort
FC=ifort
PKG_CFLAGS= -I/XXXX/magma-1.6.1/include -I/XXXX/mkl/include \
-I/XXXX/cuda/6.5/include -DHAVE_CUBLAS -DADD_ -openmp -O3
PKG_LIBS=/XXXX/magma-1.6.1/lib/libmagma.a \ -L/XXXX/mkl/lib/intel64 \
-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_vml_avx -lmkl_def \
-L/XXXX/cuda/6.5/lib64 -lcublas -lcudart -openmp
```

Be sure to substitute the specific paths on you system for "XXXX" along with 
alternative version numbers and other paths where necessary.  If you are 
using a Mac then some of the above paths will likely be slightly different 
depending on where each program was installed.  It also might be necessary 
to give the paths to icc and ifort.  You can find the paths using the 
commands:
```
which icc
which ifort
```
on the command line.

If you are using gcc on a Mac and built MAGMA with OpenBLAS on a 64 bit 
system with CUDA 7.0, the following Makevars would be appropriate:

```
CC=gcc
CXX=gcc
F77=gfortran
FC=gfortran
PKG_CFLAGS= -I/XXXX/magma-1.6.1/include -I/opt/local/include \
-I/usr/local/cuda/include -I/Developer/NVIDIA/CUDA-7.0/include -DHAVE_CUBLAS \
-DADD_ -fopenmp -std=c99
PKG_LIBS=/XXXX/magma-1.6.1/lib/libmagma.a -L/opt/local/lib \
-lopenblas -L/usr/local/cuda/lib/ -L/Developer/NVIDIA/CUDA-7.0/lib -lcudart \
-lcublas -fopenmp
```

Once again, the true paths that you could use depend on where CUDA and your 
BLAS libraries are installed on your specific system.  Make sure to substitute 
the appropriate paths for "XXXX" and the appropriate program versions where 
necessary.  It may also be necessary to use the exact paths to gcc and 
gfortran, which can be found using the commands:
```
which gcc
which gfortran
```
Overall, the PKG_CFLAGS and PKG_LIBS variables should be set to include all 
libraries used to build MAGMA.

After the Makevars file has been created and includes MAGMA and all the 
libraries used to build MAGMA, you will be able to install ff-MAGMA on 
your system.  To install from source, download the package tarball and on 
the command line run:
```
R CMD BUILD ff-MAGMA_1.0.tar.gz
R CMD INSTALL ff-MAGMA
```
substituting the appropriate version number for "1.0".

Assuming you received no errors when running the above command, you now 
have ff-MAGMA installed.

If you use RStudio on a Mac, there may still be some problems when loading 
the package.  Try to load the package from RStudio. You might receive an 
error similar to:
```
Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/Users/jpaige/Library/R/3.1/library/ff-MAGMA/libs/ff-MAGMA.so':
  dlopen(/Users/jpaige/Library/R/3.1/library/ff-MAGMA/libs/ff-MAGMA.so, 6): Library not loaded: @rpath/libcudart.7.0.dylib
  Referenced from: /Users/jpaige/Library/R/3.1/library/ff-MAGMA/libs/ff-MAGMA.so
  Reason: image not found
Error: package or namespace load failed for ‘ff-MAGMA’
```
That may be due to RStudio, since in our experience when that happens the 
package is still able to be loaded when running R in Terminal.  If you 
would like to run ff-MAGMA in RStudio there is an easy fix for this 
error.  Run ``.libPaths()`` in R to see the location(s) of the packages you have installed to R.  
Find the ff-MAGMA package directory under the given library paths 
and go to the "libs" subdirectory.  Use the ``otool -L`` command in Terminal to view the paths to the libraries used by 
ff-MAGMA. The following command will fix the above error:
```
install_name_tool -change "@rpath/libcudart.7.0.dylib" "/Developer/NVIDIA/CUDA-7.0/lib/libcudart.7.0.dylib" ff-MAGMA.so
```
substituting paths and library names as appropriate.  You may need to 
run equivalent commands for other libraries with "@rpath" in the load 
paths given by the "otool" command given above.

After you have completed the above commands and are able to load the 
ff-MAGMA library into R, congratulations! You now have access to 
a substantially accelerated version of the fields package.
