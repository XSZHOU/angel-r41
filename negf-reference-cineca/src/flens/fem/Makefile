ifndef FLENS_HOME
FLENS_HOME = $(PWD)/..
endif

include $(FLENS_HOME)/Makefile.common

SOURCEFILES        = $(filter-out fem.cc, $(wildcard *.cc))

TARGETS      = example create_mesh
LINK         = -L/sw/lib -L$(FLENS_HOME) -lflens $(LDFLAGS) $(EX_LDFLAGS)
INCDIRS     += -I/sw/include


all: $(TARGETS)

clean:
	$(RM) *.o
	$(RM) $(TARGETS)

$(TARGETS): % : %.o fem.o
	$(CXX) -o $@ fem.o $< $(LINK)

%.o : %.cc
	$(CXX) $(CXXFLAGS) $(INCDIRS) -I . -o $@ -c $<
