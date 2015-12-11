# MAKEFILE ANGEL
# COATES.RCAC.PURDUE.EDU VERSION

# ============================================================
# external include directories and libraries
# ============================================================

BASE_PATH    = /home/ba01/u107/steiger/src
LIB_PATH     = $(BASE_PATH)/lib

DEF_LIB      = 
GCC_LIB      = 

BOOST_INC    = $(BOOST_INCLUDE)

CBLAS_INC    = $(LAPACK_INCLUDE)
CBLAS_LIB    = $(LINK_LAPACK)

FLENS_INC    = -I ./flens
FLENS_LIB    = -Wl,-rpath,$(PWD)/../lib -L $(PWD)/../lib -lflens 

LINALG_INC   = 
LINALG_LIB   = $(LAPACK_INCLUDE)

TDKP_INC     = 
TDKP_LIB     = 

ZLIB_INC     = 
ZLIB         = -lz

UMFPACK_INC  = -I ./umfpack/UMFPACK -I ./umfpack/AMD/Include -I ./umfpack/UFconfig
UMFPACK_LIB  = -Wl,-rpath,$(PWD)../lib -L $(PWD)../lib -lumfpack 

DFISE_INC    =
DFISE_LIB    = 

MPI_INC      = 
MPI_LIB      = 

PYTHON_H     = $(PYTHONPATH)/include/python2.5

RAPP_INC     = -I $(BASE_PATH)/rappture/20090818/include
RAPP_LIB     = -Bdynamic -Wl,-rpath,$(BASE_PATH)/rappture/20090818/lib -L $(BASE_PATH)/rappture/20090818/lib -lrappture

NEGF_VAR     = $(BASE_PATH)/angel-nanohub

# ============================================================
# compiler settings
# ============================================================

DEFINES      = -DDEBUG -DMPICH_IGNORE_CXX_SEEK -DNOACML -DPPC -DCSCS -DNOTDKP -DNODFISE

CPP          = $(CXX)

CPP_FLAGS    = -fPIC -Wall -Wno-deprecated -g $(DEFINES)
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


FFTWPATH     = $(FFTW_HOME)
CBLASPATH    = $(LAPACK_INCLUDE) # is actually wrong.... we use MKL, not Atlas

FLENS_LDFLAGS = $(FFTW_LOADLIBES) -Wl,-rpath,/apps/rhel5/intel/mkl/10.2.2.025/lib/em64t/ -L/apps/rhel5/intel/mkl/10.2.2.025/lib/em64t/ -lmkl_lapack -lmkl_core -lmkl_sequential -lmkl_intel_lp64 -lgfortran

# -c flag is inserted individually
# define NOUNDERSCORE if XL version of lapack/blas is used
#FLENS_CXXFLAGS = -g -O3 -q64 -Wall -Wno-deprecated -Wno-deprecated-declarations -DNOUNDERSCORE
FLENS_CXXFLAGS  = -fPIC -g -O3 -Wall -Wno-deprecated -Wno-deprecated-declarations


# =========================================
# UMFPACK configuration file
# see src/umfpack/UFconfig/ directory
# =========================================
UF_MK = UFconfig.mk.coates.gnu

