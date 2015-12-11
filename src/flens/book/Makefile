include ../Makefile.common

CXXFLAGS    += -Wno-deprecated -Wno-deprecated-declarations
SOURCEFILES  = $(wildcard *.cc)
TARGETS      = $(patsubst %.cc,%,$(SOURCEFILES))
LINK         = -L.. -lflens -lwavelets $(LDFLAGS)

all: $(TARGETS)

clean:
	$(RM) *.o
	$(RM) $(TARGETS)
	
$(TARGETS): % : %.o
	$(CXX) -o $@ $< $(LINK)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCDIRS) -I . -o $@ -c $<
