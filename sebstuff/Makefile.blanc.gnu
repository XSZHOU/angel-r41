# MAKEFILE PROJECT NEGF
# CSCS-GNU VERSION

# ============================================================
# external include directories and libraries
# ============================================================

BASE_PATH    = /users/steiger/src/release/gnu
LIB_PATH     = $(BASE_PATH)/lib

DEF_LIB      = -lc -lrt
GCC_LIB      = -lstdc++ -L/opt/ibmhpc/ppe.poe/lib/libmpi64 -lmpi_ibm -L/opt/ibmhpc/ppe.poe/lib/libmpi64 -lpoe -L/usr/lib64 -llapi -lgfortranbegin -lgfortran -lm -shared-libgcc
# the latter was copy-pasted from "mpCC -compiler g++ -v" and "mpfort -compiler gfortran -v"

BOOST_INC    = -I/users/steiger/src/master/boost_1_35_0

CBLAS_INC    = -I$(BASE_PATH)/cblas/src
CBLAS_LIB    = -L$(LIB_PATH) -lcblas

FLENS_INC    = -I$(BASE_PATH)/FLENS-lite_0608
FLENS_LIB    = -L$(LIB_PATH) -lflens 

LINALG_INC   = 
#LINALG_LIB  = /apps/lapack/3.1.1/XL/lib/lapack_ppc64.a /apps/lapack/3.1.1/XL/lib/blas_ppc64.a 
LINALG_LIB   = /apps/lapack/3.1.1/GNU/lapack_ppc64.a /apps/lapack/3.1.1/GNU/blas_ppc64.a

TDKP_INC     = -I$(BASE_PATH)/tdkp
TDKP_LIB     = -L$(LIB_PATH) -ltdkp -ljdqz

ZLIB_INC     = -I/apps/zlib/1.2.3/XL/include
ZLIB         = -L/apps/zlib/1.2.3/XL/lib -lz

UMFPACK_INC  = -I$(BASE_PATH)/umfpack/UMFPACK -I$(BASE_PATH)/umfpack/AMD/Include -I$(BASE_PATH)/umfpack/UFconfig
UMFPACK_LIB  = -L$(LIB_PATH) -lumfpack -lamd

#DFISE_INC    = -I$(BASE_PATH)/DF-ISE
#DFISE_LIB    = -L$(LIB_PATH) -lDF-ISE++ -lDF-ISE
DFISE_INC    = -I$(BASE_PATH)/sebise/sebise
DFISE_LIB    = -L$(LIB_PATH) -lsebise

MPI_INC      = 
MPI_LIB      = 

# ============================================================
# compiler settings
# ============================================================

# defines
DEFINES      = -DDEBUG -DMPICH_IGNORE_CXX_SEEK -DNOACML -DPPC -DCSCS

CPP          = mpCC -compiler g++

CPP_FLAGS    = -q64 -Wall -Wno-deprecated -g -O1 $(DEFINES) -m64

CPP_INC      = -Isrc/common/ -Isrc/geometry/ -Isrc/io/ -Isrc/newton/ -Isrc/newton/solvers/ -Isrc/negf/ \
			   -Isrc/negf/selfenergies/ -Isrc/negf/hamiltonian/ -Isrc/kspace_Espace/ \
			   $(TDKP_INC) $(BOOST_INC) $(FLENS_INC) $(CBLAS_INC) $(MPI_INC)

CPP_ALL_INC  = $(CPP_INC) $(DFISE_INC) $(UMFPACK_INC) $(ZLIB_INC) $(LINALG_INC) 

CPP_LIBS     = -nodefaultlibs -Bstatic  $(TDKP_LIB) $(DFISE_LIB) $(UMFPACK_LIB) $(ZLIB) $(FLENS_LIB) $(CBLAS_LIB) $(LINALG_LIB) $(MPI_LIB) $(GCC_LIB) $(DEF_LIB) 

# ===============================================================
# from here on the makefile is platform- and compiler-independent
# ===============================================================

include definitions.negf
