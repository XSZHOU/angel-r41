# FILE EDITED BY SS
# blanc.cscs.ch version

BASEDIR      = /users/steiger/src/release/gnu

FFTWPATH     = /apps/fftw/3.1.3/XL
CBLASPATH    = $(BASEDIR)/lib

# gfortran: -L/usr/lib -llapi -lgfortranbegin -lgfortran -lm -shared-libgcc

CXX          = mpCC -compiler g++
INCDIRS     += -I. -I$(BASEDIR) -I$(BASEDIR)/cblas/src -I$(FFTWPATH)/include 

# switch between XL and GNU versions!
#LDFLAGS     += /apps/lapack/3.1.1/XL/lib/lapack_ppc64.a $(CBLASPATH)/libcblas.a /apps/lapack/3.1.1/XL/lib/blas_ppc64.a $(FFTWPATH)/lib/libfftw3.a
LDFLAGS     += /apps/lapack/3.1.1/GNU/lapack_ppc64.a $(CBLASPATH)/libcblas.a /apps/lapack/3.1.1/GNU/blas_ppc64.a $(FFTWPATH)/lib/libfftw3.a

# <ss> NEW: compiler libraries should be OK anyway
#LDFLAGS    += -L/opt/ibmcmp/xlf/11.1/lib64 -lxlf90_r -lxl -lxlfmath -L/opt/ibmcmp/xlsmp/1.7/lib64 -lxlomp_ser

DYLIB_EXT    = a
CXXDYLIB     = 
# -c flag is inserted individually
# define NOUNDERSCORE if XL version of lapack/blas is used
#CXXFLAGS    += -g -O3 -q64 -Wall -Wno-deprecated -Wno-deprecated-declarations -DNOUNDERSCORE
CXXFLAGS    += -g -O3 -q64 -Wall -Wno-deprecated -Wno-deprecated-declarations

#LINKER      = $(CXX)
LINKER       = ar r
#CXXLIBFLAGS  = $(CXXLIBFLAGS) -o
CXXLIBFLAGS  = 

