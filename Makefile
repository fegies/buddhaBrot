#The program Name (Name of the resulting executable)

PROG       = buddhaBrot
VPATH      = src include
ODIR       = ./bin
OBJDIR     = objs
SHAREFLAGS = -pipe -Wall -pedantic -O2
CPPCFLAGS  = $(SHAREFLAGS)
CCFLAGS    = $(SHAREFLAGS) -std=gnu11
LINKFLAGS  = $(SHAREFLAGS) -lpthread
CPPCOMPILER= $(CXX)
CCOMPILER  = $(CC)

#Name of subpaths inside Odir (Must be the same in ./include and ./src as well)
SUBPATHS   = .

#The Objects that are compiled
OBJS       = $(BASEOBS)

BASEOBS    = main.o qdbmp.o xorshift.o

OPROG = $(addprefix $(ODIR)/, $(PROG))
RUNFLAGS = -o buddha.bmp -p 8

DEPDIR := deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
COMPILE.c = $(CCOMPILER) $(DEPFLAGS) $(CCFLAGS) -c
COMPILE.cpp = $(CPPCOMPILER) $(DEPFLAGS) $(CPPCFLAGS) -c
POSTCOMPILE= mv -f $(DEPDIR)/$*.TD $(DEPDIR)/$*.d

run : all
	time $(OPROG) $(RUNFLAGS)

all : buildbin $(OPROG)

flagless : all
	$(OPROG)

clean:
	rm -rf ./bin $(DEPDIR)

.PHONY: run all clean flagless
	

#linking
$(OPROG): $(addprefix $(ODIR)/$(OBJDIR)/, $(OBJS))
	$(CPPCOMPILER) $(LINKFLAGS) -o $@ $^

#compiling
$(ODIR)/$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cpp) -o $@ $< -I./include

$(ODIR)/$(OBJDIR)/%.o : %.c $(DEPDIR)/%.d
	$(COMPILE.c) -o $@ $< -I./include

#Building the Directories if they don't exist
buildbin: | \
	$(addprefix $(ODIR)/$(OBJDIR)/,$(SUBPATHS)) \
	$(addprefix $(DEPDIR)/,$(SUBPATHS))

$(addprefix $(ODIR)/$(OBJDIR)/,$(SUBPATHS)) \
	$(addprefix $(DEPDIR)/,$(SUBPATHS)):
	mkdir -p $@

$(addprefix $(DEPDIR)/,$(OBJS:.o=.d)): ;
.PRECIOUS: $(DEPDIR)/%.d

-include $(addprefix $(DEPDIR)/,$(OBJS:.o=.d))
