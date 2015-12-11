include ../Makefile.common

SOURCEFILES  = $(wildcard *.cc)
TARGETS      = $(patsubst %.cc,%,$(SOURCEFILES))
LINK         = -L.. -lflens $(LDFLAGS)

all: $(TARGETS)

clean:
	$(RM) *.dat
	$(RM) *.o
	$(RM) $(TARGETS)

$(TARGETS): % : %.o
	$(MXX) -o $@ $< $(LINK)

%.o: %.cc
	$(MXX) $(MXXFLAGS) $(INCDIRS) -I . -o $@ -c $<


