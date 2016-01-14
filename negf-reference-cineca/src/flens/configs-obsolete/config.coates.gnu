# FILE EDITED BY SS
# coates.rcac.purdue.edu version

BASEDIR      = /home/ba01/u125/steiger/src

FFTWPATH     = $(FFTW_HOME)
CBLASPATH    = $(LAPACK_INCLUDE) # is actually wrong....

# gfortran: -L/usr/lib -llapi -lgfortranbegin -lgfortran -lm -shared-libgcc

# already set CXX 
INCDIRS     += -I. -I$(BASEDIR) $(CBLASPATH) -I$(FFTWPATH)/include 

LDFLAGS     += $(FFTW_LOADLIBES) -Wl,-rpath,/apps/rhel5/intel/mkl/10.2.2.025/lib/em64t/ -L/apps/rhel5/intel/mkl/10.2.2.025/lib/em64t/ -lmkl_lapack -lmkl_core -lmkl_sequential -lmkl_intel_lp64 -lgfortran

#DYLIB_EXT    = a
DYLIB_EXT    = so
CXXDYLIB     = -shared
# -c flag is inserted individually
# define NOUNDERSCORE if XL version of lapack/blas is used
#CXXFLAGS    += -g -O3 -q64 -Wall -Wno-deprecated -Wno-deprecated-declarations -DNOUNDERSCORE
CXXFLAGS    += -fPIC -g -O3 -Wall -Wno-deprecated -Wno-deprecated-declarations

# not used 
LINKER       = ar r
CXXLIBFLAGS  = 

