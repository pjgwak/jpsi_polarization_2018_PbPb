PROGRAM = Analysis

SRCDIR    = src
OBJDIR    = obj

HDRS = $(SRCDIR)/ModelBuilder.h\
			 $(SRCDIR)/JpsiFitter.h\

OBJS = $(OBJDIR)/ModelBuilder.o\
			 $(OBJDIR)/JpsiFitter.o\
       $(OBJDIR)/dict_$(PROGRAM).o # need for ROOT and pyroot

#-------------------------------------------------------------

ObjSuf = o # suffix of obj
SrcSuf = cpp
ExeSuf = exe
DllSuf = so
OutPutOpt = -o

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

# directories
HEADRDIR = .

# debug option: -g
# no exceptions: -fno-exceptions

CXX = g++
CXXFLAGS = -O2 -Wall -fPIC -g# -std=c++17
LD = g++
LDFLAGS = -O2 -Wall -g
SOFLAGS = -shared
ARCHFLAGS =

CXXFLAGS += $(ROOTCFLAGS) -I$(HEADRDIR)
LIBS = $(ROOTLIBS)
GLIBS = $(ROOTGLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) -c $(CXXFLAGS) $< $(OutPutOpt) $@

$(PROGRAM).so: $(OBJS)
	@echo "Linking $(PROGRAM).so ..."
	$(LD) $(ARCHFLAGS) $(LDFLAGS) $(SOFLAGS) $(OBJS) $(ROOTLIBS) -lRooStats -o$(PROGRAM).so
	@echo "linking done"

$(OBJDIR)/dict_$(PROGRAM).o: $(OBJDIR)/dict_$(PROGRAM).cpp | $(OBJDIR)
	$(CXX) -c $(CXXFLAGS) $< $(OutPutOpt) $@

$(OBJDIR)/dict_$(PROGRAM).cpp: $(HDRS) | $(OBJDIR)
	@echo "Generating dictionary..."
	rootcling -f $(OBJDIR)/dict_$(PROGRAM).cpp -c -p -I$(HEADRDIR) $(HDRS) Linkdef.h
	@cp $(OBJDIR)/dict_$(PROGRAM)_rdict.pcm .

# make an obj folder if it doesn't exist
$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	@rm -f $(OBJS) $(PROGRAM).so $(OBJDIR)/dict_$(PROGRAM).cpp dict_$(PROGRAM)_rdict.pcm