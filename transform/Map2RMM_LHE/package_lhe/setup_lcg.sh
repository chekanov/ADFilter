#!/bin/bash
echo "Setup ROOT, PyROOT tensorflow"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "views LCG_103 x86_64-centos9-gcc11-opt"
#bash -c "source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos9-gcc11-opt/setup.sh"
echo "Setup Fastjet"
# bash -c "source /cvmfs/sft.cern.ch/lcg/latest/fastjet/3.4.1-5af57/x86_64-el9-gcc11-opt/fastjet-env.sh"
export FASTJET=/cvmfs/sft.cern.ch/lcg/latest/fastjet/3.4.1-5af57/x86_64-el9-gcc11-opt/

#ls /cvmfs/sft.cern.ch/lcg/latest/

source /var/www/html/asc/adfilter/transform/Map2RMM/library/ProMCBin/promc/setup.sh


export PYTHIA8=/cvmfs/sft.cern.ch/lcg/latest/MCGenerators/pythia8/311-f5a58/x86_64-el9-gcc11-opt/
export PYTHIADIR=$PYTHIA8
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
# LHAPDF6 configuration.
export LHAPDF6_USE=true

export LHAPDF6=/cvmfs/sft.cern.ch/lcg/latest/MCGenerators/lhapdf/6.5.3-480be/x86_64-el9-gcc11-opt

export LHAPDF6_BIN=$LHAPDF6/bin/
export PATH=$LHAPDF6_BIN:$PATH
export LHAPDF6_INCLUDE=$LHAPDF6/include/
export LHAPDF6_LIB=$LHAPDF6/lib
export LD_LIBRARY_PATH=$LHAPDF6_LIB/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH

export HEPMC=$MSOFT/HEPMC
export LD_LIBRARY_PATH=$HEPMC/lib:$LD_LIBRARY_PATH



#source  /cvmfs/sft.cern.ch/lcg/latest/fastjet/3.4.1-5af57/x86_64-el9-gcc11-opt/fastjet-env.sh 

