PROGRAM = Analysis

SRCDIR    = src
OBJDIR    = obj

HDRS = $(SRCDIR)/McMassFitter.h\
			 $(SRCDIR)/MassFitter.h\
			 $(SRCDIR)/SPlotter.h\

OBJS = $(OBJDIR)/McMassFitter.o \
			 $(OBJDIR)/MassFitter.o \
			 $(OBJDIR)/SPlotter.o \
			 $(OBJDIR)/polarizationUtilities.o\
			 $(OBJDIR)/dict_$(PROGRAM).o


ObjSuf = o
SrcSuf = cpp
DllSuf = so
OutPutOpt = -o

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

# compiler and flags
CXX = g++
CXXFLAGS = -O2 -Wall -fPIC -g
LD = g++
LDFLAGS = -O2 -Wall -g
SOFLAGS = -shared

CXXFLAGS += $(ROOTCFLAGS) -I. -I$(SRCDIR) -I../../../headers
LIBS = $(ROOTLIBS)

# ===== compile rules =====
# .cpp to .o
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) -c $(CXXFLAGS) $< $(OutPutOpt) $@

# combine objects to .so
$(PROGRAM).so: $(OBJS)
	@echo "Linking $(PROGRAM).so ..."
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJS) $(ROOTLIBS) -lRooStats -o$(PROGRAM).so
	@echo "linking done"

# compile files out of src folder
obj/polarizationUtilities.o: ../../../headers/polarizationUtilities.cpp | $(OBJDIR)
	$(CXX) -c $(CXXFLAGS) $< -o $@
	

# ===== dictioranry rules =====
# dictionary source compile
$(OBJDIR)/dict_$(PROGRAM).o: $(OBJDIR)/dict_$(PROGRAM).cpp | $(OBJDIR)
	$(CXX) -c $(CXXFLAGS) $< $(OutPutOpt) $@

# make dictionary.cpp from Linkdef.h and headers
$(OBJDIR)/dict_$(PROGRAM).cpp: $(HDRS) Linkdef.h | $(OBJDIR)
	@echo "Generating dictionary..."
	rootcling -f $(OBJDIR)/dict_$(PROGRAM).cpp -c -p -I$(SRCDIR) $(HDRS) Linkdef.h
	@cp $(OBJDIR)/dict_$(PROGRAM)_rdict.pcm .

# make a obj folder if not obj
$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	@rm -rf $(OBJDIR) $(PROGRAM).so dict_$(PROGRAM)_rdict.pcm