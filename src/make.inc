# MAKEFILE ANGEL
# KUBUNTU 10.04 VERSION

# ============================================================
# external include directories and libraries
# ============================================================

BASE_PATH    = /Users/seanzhou/Documents/angel-r41
LIB_PATH     = $(BASE_PATH)/lib

DEF_LIB      = 
GCC_LIB      = 

BOOST_INC    = -I /usr/local/Cellar/boost/1.58.0/include

#CBLAS_INC    = -I /usr/local/include
CBLAS_INC    = -I /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
CBLAS_LIB    = -L /usr/local/lib

#here /flens, not /flens/flens, otherwise"flens/flens.h: No such file or directory"
FLENS_INC    = -I./flens
FLENS_LIB    = -framework Accelerate -Wl,-rpath,$(PWD)/../lib -L $(PWD)/../lib -lflens 

LINALG_INC   = 
LINALG_LIB   =  -L /usr/local/lib

TDKP_INC     = 
TDKP_LIB     = 

ZLIB_INC     = 
ZLIB         = -lz

UMFPACK_INC  = -I./umfpack/UMFPACK -I./umfpack/AMD/Include -I./umfpack/UFconfig
UMFPACK_LIB  = -Wl,-rpath,$(PWD)/../lib -L $(PWD)/../lib -lumfpack 

DFISE_INC    = 
DFISE_LIB    = 

MPI_INC      = 
MPI_LIB      = 

PYTHON_H     = $(PYTHONPATH)/include/python2.5

RAPP_INC     = -I /home/steiger/src/rappture/20090818/include
RAPP_LIB     = -Bdynamic -Wl,-rpath,/home/steiger/src/rappture/20090818/lib -L /home/steiger/src/rappture/20090818/lib -lrappture

NEGF_VAR     = $(BASE_PATH)

# ============================================================
# compiler settings
# ============================================================
# S.Z. cannot put -DVECLIB here, otherwise conflict functions !
DEFINES      = -DDEBUG -DMPICH_IGNORE_CXX_SEEK -DNOACML -DPPC -DCSCS -DNOTDKP -DNODFISE

#should not use -std=c++11, since the original code by steiger was written before c++11 standard
CPP          = mpic++  

CPP_FLAGS    = -framework Accelerate -fPIC -Wall -Wno-deprecated -g $(DEFINES)
# -O1 gives weird memory error in inlined Logger->emit() routine

CPP_INC      = -Icommon/ -Igeometry/ -Iio/ -Ipoisson/ -Inegf/ \
			   -Inegf/selfenergies/ -Inegf/hamiltonian/ -Ikspace_Espace/ \
			   $(TDKP_INC) $(BOOST_INC) $(FLENS_INC) $(CBLAS_INC) $(MPI_INC)

CPP_ALL_INC  = $(CPP_INC) $(DFISE_INC) $(UMFPACK_INC) $(ZLIB_INC) $(LINALG_INC) 

CPP_LIBS     = $(TDKP_LIB) $(DFISE_LIB) $(UMFPACK_LIB) $(ZLIB) $(FLENS_LIB) $(CBLAS_LIB) $(LINALG_LIB) $(MPI_LIB) $(GCC_LIB) $(DEF_LIB) 

SWIG         = swig


# =========================================
# FLENS: used in src/flens/Makefile.common
# =========================================

FLENS_RM     = rm -f
FLENS_RMDIR  = rm -f -r
FLENS_MKDIR  = mkdir -p
FLENS_MAKE   = make
FLENS_SILENT = @

FLENS_BASEDIR = $(BASE_PATH)
FLENS_PREFIX  = $(FLENS_BASEDIR)/angel-nanohub/src/flens/


FFTWPATH     = /usr/include/
CBLASPATH    = -I/usr/local/include
# this is enough for resolving "cblas.h" error during compiling for libflens.so, no need for "-framework Accerlerate"
CBLASPATH   += -I/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers

#FLENS_LDFLAGS = -L /usr/lib -lfftw3 -llapack -lblas -lcblas   ->S.Z. -lblas/cblas should be the OSX shipped version VECLIB
FLENS_LDFLAGS = -framework Accelerate -L/usr/local/lib -lfftw3

# -c flag is inserted individually
# define NOUNDERSCORE if XL version of lapack/blas is used
#FLENS_CXXFLAGS = -fPIC -g -O1 -Wall -Wno-deprecated -Wno-deprecated-declarations -DNOUNDERSCORE
FLENS_CXXFLAGS  = -fPIC -g -O1 -Wall -Wno-deprecated -Wno-deprecated-declarations


# =========================================
# UMFPACK configuration file
# see src/umfpack/UFconfig/ directory
# =========================================
UF_MK = UFconfig.mk.kubuntu


