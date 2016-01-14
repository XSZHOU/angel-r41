# MAKEFILE ANGEL
# KUBUNTU 10.04 VERSION

# ============================================================
# external include directories and libraries
# ============================================================

BASE_PATH    = /galileo/home/userexternal/xzhou000/test_negf/negf-reference-master
LIB_PATH     = $(BASE_PATH)/lib

DEF_LIB      = 
GCC_LIB      = -lmkl_intel_lp64

BOOSTs_INC    = -I${BOOST_INC}

CBLAS_INC    = -I${MKL_INC} 
CBLAS_LIB    = -L${MKL_LIB} -lmkl_intel_lp64

#here /flens, not /flens/flens, otherwise"flens/flens.h: No such file or directory"
FLENS_INC    = -I./flens
FLENS_LIB    = -Wl,-rpath,$(PWD)/../lib -L $(PWD)/../lib -lflens 

#LINALG_INC   = 
#LINALG_LIB   =  -L /usr/local/lib

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

NEGF_VAR     = $(BASE_PATH)

# ============================================================
# compiler settings
# ============================================================
# S.Z. cannot put -DVECLIB here, otherwise conflict functions !
DEFINES      = -DMKL -DDEBUG -DMPICH_IGNORE_CXX_SEEK -DNOACML -DPPC -DCSCS -DNOTDKP -DNODFISE

#should not use -std=c++11, since the original code by steiger was written before c++11 standard
CPP          = mpiicpc 

CPP_FLAGS    = -DMKL -DMKL_ILP64 -fPIC -Wall -Wno-deprecated -g $(DEFINES)
# -O1 gives weird memory error in inlined Logger->emit() routine

CPP_INC      = -Icommon/ -Igeometry/ -Iio/ -Ipoisson/ -Inegf/ \
			   -Inegf/selfenergies/ -Inegf/hamiltonian/ -Ikspace_Espace/ \
			   $(TDKP_INC) $(BOOSTs_INC) $(FLENS_INC) $(CBLAS_INC) $(MPI_INC)

CPP_ALL_INC  = $(CPP_INC) $(DFISE_INC) $(UMFPACK_INC) $(ZLIB_INC) $(CBLAS_INC) 

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
FLENS_PREFIX  = $(FLENS_BASEDIR)/src/flens/


FFTWPATH     = -I${FFTW_INC}
CBLASPATH    = -I${MKL_INC} 

FLENS_LDFLAGS =  -L${MKL_LIB} -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

# -c flag is inserted individually
# define NOUNDERSCORE if XL version of lapack/blas is used
#FLENS_CXXFLAGS = -fPIC -g -O1 -Wall -Wno-deprecated -Wno-deprecated-declarations -DNOUNDERSCORE
FLENS_CXXFLAGS  = -DMKL -fPIC -g -O1 -Wall -Wno-deprecated


# =========================================
# UMFPACK configuration file
# see src/umfpack/UFconfig/ directory
# =========================================
# This is importrant, S.Z.
UF_MK = UFconfig.mk.kubuntu


