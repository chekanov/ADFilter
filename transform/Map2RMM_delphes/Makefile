# Makefile example for Promc+ROOT+fastjet
# S.Chekanov (ANL) 

ifndef ROOTSYS 
$(error ROOTSYS env variable is not set. Install ROOT first. See https://root.cern.ch/)
endif

include ${ROOTSYS}/etc/Makefile.arch
include ${PROMC}/etc/config.mk

# Root variables
ROOTCFLAGS    = $(shell root-config --nonew --cflags)
ROOTLIBS      = $(shell root-config --nonew --libs)
ROOTGTTLIBS   = $(shell root-config --nonew --glibs)
# Assign or add variables
CXXFLAGS     += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS)

LIBS += -L./map2rmm/lib -lmap2rmm
INCLUDE1= -I./inc -I./
INCLUDE2= -I./src -I./map2rmm/inc/


Tasks:     clean example

SOURCE_FILES := $(shell ls -1 example.cc)

# build object files 
objects       = $(patsubst %.cc,%.o,$(SOURCE_FILES))


%.o: %.cc
	$(CXX) $(OPT) -Wno-unused-variable  -Wunused-but-set-variable $(CXXFLAGS) $(INCLUDE1) $(INCLUDE2) $(INCLUDE3) -o $@ -c $<

LIBOBJS = $(patsubst %.cc,%.o,$(SOURCE_FILES))

example: $(objects)
	$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
clean:
	        @rm -f *.o example *~;  echo "Clear.." 
