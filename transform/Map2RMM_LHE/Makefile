# Makefile example for Promc+ROOT+fastjet
# S.Chekanov (ANL) 

ifndef ROOTSYS 
$(error ROOTSYS env variable is not set. Install ROOT first. See https://root.cern.ch/)
endif

ifndef PROMC 
  $(error ProMC file format to read events from the HepSim repository is not found.   PROMC env variable is not set. Install ProMC first. See https://atlaswww.hep.anl.gov/asc/promc/)
endif

ifndef FASTJET
  $(error FASTJET  env variable is not set. Install FASTJET first. See http://fastjet.fr/)
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
LIBS += -L${FASTJET}/lib -lfastjet -L./map2rmm/lib -lmap2rmm
LIBS += -L${PROMC}/lib -lpromc -lprotoc -lprotobuf -lprotobuf-lite -lcbook -lz


INCLUDE1= -I./inc -I./
INCLUDE2= -I./src -I./map2rmm/inc/
INCLUDE3= -I${PROMC}/include -I$(PROMC)/src -I${FASTJET}/include 


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
