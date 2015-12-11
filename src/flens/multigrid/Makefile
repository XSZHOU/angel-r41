include ../Makefile.common

INCDIRS += -I . 
OBJDIR  = .obj/
DEPDIR  = .dep/
SOURCES := $(wildcard *.cc)
HEADERS := $(wildcard *.h) $(wildcard *.tcc)
OBJECTS = $(addprefix $(OBJDIR),$(SOURCES:.cc=.o))
DEPS    = $(addprefix $(DEPDIR),$(SOURCES:.cc=.d))
LDFLAGS += -L.. -lflens

all: $(OBJDIR) $(DEPDIR) $(OBJECTS)
	
install: $(OBJECTS)
	$(SILENT) $(INSTALL) $(HEADERS) $(INCDIR)
	
$(OBJDIR):
	$(SILENT) $(MKDIR) $(OBJDIR)
	
$(DEPDIR):
	$(SILENT) $(MKDIR) $(DEPDIR)
	
$(OBJDIR)%.o:
	$(SILENT) $(MXX) $(INCDIRS) $(MXXFLAGS) -o $@ -c $<

$(DEPDIR)%.d: %.cc $(DEPDIR)
	$(SILENT) $(MXX) -MM $(INCDIRS) $(MXXFLAGS) $< \
				| sed 's/$*.o/$(OBJDIR:/=\/)$*.o/g' \
				| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

clean:
	$(SILENT) $(RMDIR) .obj/ .dep/
	
-include $(addprefix $(DEPDIR),$(SOURCES:.cc=.d))
