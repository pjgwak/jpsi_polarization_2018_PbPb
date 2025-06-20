#-------------------------------------------------------------

PROGRAM = Analysis

HDRS = EventLoop.h\

OBJS = EventLoop.o\
			 dict_$(PROGRAM).o # need for ROOT and pyroot

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
CXXFLAGS = -O3 -Wall -fPIC -g -std=c++17
LD = g++
LDFLAGS = -O3 -Wall -g
SOFLAGS = -shared
ARCHFLAGS = 

CXXFLAGS += $(ROOTCFLAGS) -I$(HEADRDIR)
LIBS = $(ROOTLIBS)
GLIBS = $(ROOTGLIBS)

%/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< $(OutPutOpt) $@

$(PROGRAM).so: $(OBJS)
	@echo "Linking $(PROGRAM).so ..."
	$(LD) $(ARCHFLAGS) $(LDFLAGS) $(SOFLAGS) $(OBJS) $(ROOTLIBS) -o$(PROGRAM).so
	@echo "linking done"

dict_$(PROGRAM).o: dict_$(PROGRAM).cpp

dict_$(PROGRAM).cpp: $(HDRS)
	@echo "Generating dictionary..."
	rootcint -f dict_$(PROGRAM).cpp -c -p -I$(HEADRDIR) $(HDRS) Linkdef.h

clean:
	@rm -f $(OBJS) $(PROGRAM).so dict_*

#-------------------------------------------------------------