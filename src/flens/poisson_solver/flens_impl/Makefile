ifndef FLENS_HOME
FLENS_HOME = $(PWD)/../..
endif

include $(FLENS_HOME)/Makefile.common

DYLIB       = libpoisson.$(DYLIB_EXT)
DYLIB_MPI   = libpoisson_mpi.$(DYLIB_EXT)

OBJDIR  = .obj/
DEPDIR  = .dep/
ALL_SOURCES = $(wildcard *.cc)
MPI_SOURCES = $(wildcard mpi*.cc) $(wildcard *mpi.cc) $(wildcard distributed*.cc)
SOURCE      = $(filter-out $(MPI_SOURCES), $(ALL_SOURCES))

HEADERS     = $(wildcard *.h) $(wildcard *.tcc)
OBJECTS     = $(addprefix $(OBJDIR),$(SOURCE:.cc=.o))
MPI_OBJECTS = $(addprefix $(OBJDIR),$(MPI_SOURCES:.cc=.o))
DEPS        = $(addprefix $(DEPDIR),$(ALL_SOURCES:.cc=.d))

INCDIRS     += -I/sw/include

.PHONY : all
all: $(OBJDIR) $(DEPDIR) $(OBJECTS) $(MPI_OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) -o $(DYLIB) $(OBJECTS) $(LDFLAGS) -L$(FLENS_HOME) -lflens -L/sw/lib -l fftw3
	cp $(DYLIB) $(FLENS_HOME)
	$(MXX) $(MXXDYLIB) $(CXXFLAGS) -o $(DYLIB_MPI) $(OBJECTS) $(MPI_OBJECTS) $(LDFLAGS) -L$(FLENS_HOME) -lflens -L/sw/lib -l fftw3
	cp $(DYLIB_MPI) $(FLENS_HOME)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

.PHONY : install
install: $(OBJECTS)
	$(INSTALL) $(HEADERS) $(INCDIR)
	
$(OBJDIR):
	$(MKDIR) $(OBJDIR)
	
$(DEPDIR):
	$(MKDIR) $(DEPDIR)
	
$(OBJECTS): $(OBJDIR)%.o: %.cc
	$(CXX) $(INCDIRS) $(CXXFLAGS) -o $@ -c $<

$(MPI_OBJECTS): $(OBJDIR)%.o: %.cc
	$(MXX) $(INCDIRS) $(CXXFLAGS) -o $@ -c $<

$(DEPDIR)%.d: %.cc $(DEPDIR)
	$(CXX) -c -MM $(INCDIRS) $(CXXFLAGS) $< \
		| sed 's/$*.o/$(OBJDIR:/=\/)$*.o/g' \
	| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

.PHONY : clean
clean:
	$(RMDIR) .obj/ .dep/
	$(RM) $(DYLIB)
	
