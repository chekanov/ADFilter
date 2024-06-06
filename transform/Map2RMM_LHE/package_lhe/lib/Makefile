# S.Chekanov
# ANL analysis tutorial 
include ${ROOTSYS}/test/Makefile.arch


# Root variables
ROOTCFLAGS    = $(shell root-config --nonew --cflags)
ROOTLIBS      = $(shell root-config --nonew --libs)
ROOTGTTLIBS   = $(shell root-config --nonew --glibs)
# Assign or add variables
CXXFLAGS     += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS)

SOURCE_FILES := $(shell ls -1 *.cxx)
SOURCE_FILES += $(shell ls -1 src/*.cxx)

INCLUDE1= -I./inc
INCLUDE2= -I.

# build object files 
objects       = $(patsubst %.cxx,%.o,$(SOURCE_FILES))


%.o: %.cxx
	$(CXX) $(OPT) $(CXXFLAGS) $(INCLUDE1) $(INCLUDE2) -o $@ -c $<

Tasks:    mydict main


tar:
	rm -f ana.tgz; cd ../; tar -cvzf ana.tgz --exclude='*.o' --exclude='*.root' \
	--exclude='*.log'  --exclude='*.tup' --exclude='*.eps' --exclude='*.png' --exclude='*.log.gz'\
	--exclude='main' --exclude='*.root.log.gz' --exclude='*.so' --exclude='*~' --exclude='*.pyc' ana; \
        mv ana.tgz ana/; cd ana; 

mydict: inc/LParticle.h
	@rm -f src/CParticle_dict*
	@rm -f inc/CParticle_dict*
	@echo "Generating dictionary for CParticle"
	@rootcint src/CParticle_dict.cxx -c inc/CParticle.h
	@rm -f src/LParticle_dict*
	@rm -f inc/LParticle_dict*
	@echo "Generating dictionary for LParticle"
	@rootcint src/LParticle_dict.cxx -c inc/LParticle.h

lib:  $(objects)
	@echo "Making a library";  
	ar rcs libshapes.a $(objects) 

main: $(objects)
	$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
clean:
	        @rm -f *.o  main *~ src/*.o src/*~ inc/*~; echo "Clear.." 
