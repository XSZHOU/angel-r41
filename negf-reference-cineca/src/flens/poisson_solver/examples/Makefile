ifndef FLENS_HOME
FLENS_HOME = $(PWD)/../..
endif

include $(FLENS_HOME)/Makefile.common

SOURCEFILES        = $(wildcard *.cc)
MPI_SOURCEFILES    = $(filter %_mpi.cc, $(SOURCEFILES))
NOMPI_SOURCEFILES  = $(filter-out %_mpi.cc,$(SOURCEFILES))

TARGETS      = $(patsubst %.cc,%,$(NOMPI_SOURCEFILES))
MPI_TARGETS  = $(patsubst %.cc,%,$(MPI_SOURCEFILES))
LINK         = -L/sw/lib -L$(FLENS_HOME) -lpoisson -lflens $(LDFLAGS) $(EX_LDFLAGS) -lfftw3
LINK_MPI     = -L/sw/lib -L$(FLENS_HOME) -lpoisson_mpi -lflens $(LDFLAGS) $(EX_LDFLAGS) -lfftw3
INCDIRS     += -I/sw/include


all: $(TARGETS) $(MPI_TARGETS)

clean:
	$(RM) *.o
	$(RM) $(TARGETS) $(MPI_TARGETS)

$(TARGETS): % : %.cc
	$(CXX) $(CXXFLAGS) $(INCDIRS) -I . -o $@ $< $(LINK)

$(MPI_TARGETS): % : %.cc
	$(MXX) -DUSE_MPI $(CXXFLAGS) $(INCDIRS) -I . -o $@ $< $(LINK_MPI)
