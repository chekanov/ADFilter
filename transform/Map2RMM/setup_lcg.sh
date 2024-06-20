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

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


source $SCRIPT_DIR/library/ProMCBin/promc/setup.sh

# source /var/www/html/asc/adfilter/transform/Map2RMM/library/ProMCBin/promc/setup.sh
# source /afs/cern.ch/user/c/chekanov/adfilter/transform/Map2RMM/library/ProMCBin/promc/setup.sh 
#source  /cvmfs/sft.cern.ch/lcg/latest/fastjet/3.4.1-5af57/x86_64-el9-gcc11-opt/fastjet-env.sh 

