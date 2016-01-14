# FILE EDITED BY SS
# palu.cscs.ch version
# be sure that you use gcc:
# module swap PrgEnv-pgi PrgEnv-gnu 

BASEDIR = /users/steiger/src/release/gnu

ACMLPATH  = /opt/acml/4.0.1a/gfortran64
FFTWPATH  = /apps/fftw/fftw-3.1.2
CBLASPATH = $(BASEDIR)/cblas/lib

GCCVERSION = 4.2.4
GCCPATH = /opt/gcc/$(GCCVERSION)/cnos
GCCINC1 = $(GCCPATH)/lib/gcc/x86_64-suse-linux/$(GCC_VERSION)/include/
GCCINC2 = $(GCCPATH)/include/g++/tr1
GCCINC3 = /usr/lib64/gcc-lib/x86_64-suse-linux/3.3.3/include

CXX         = CC
INCDIRS    += -I. -I$(BASEDIR) -I$(ACMLPATH)/include -I$(BASEDIR)/cblas/src -I$(FFTWPATH)/include 
#LDFLAGS    += -L$(ACMLPATH)/lib -lacml -lacml_mv -lgfortran -L$(CBLASPATH) -lcblas -L$(FFTWPATH) -lfftw3 
LDFLAGS    += -L$(ACMLPATH)/lib -lacml_mv -lgfortran -L$(CBLASPATH) -lcblas -L$(FFTWPATH) -lfftw3 

#DYLIB_EXT   = so
DYLIB_EXT   = a
#CXXDYLIB    = -shared
CXXDYLIB    = 
# -c flag is inserted individually
CXXFLAGS    += -Wall -g -O3 -Wall -Wno-deprecated -Wno-deprecated-declarations
LINKER      = ar r
CXXLIBFLAGS += 

