#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase # use your path
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
#lsetup "views LCG_102b x86_64-centos7-gcc11-opt"
#lsetup "views LCG_103 x86_64-centos7-gcc11-opt"
lsetup "views LCG_103 x86_64-centos9-gcc11-opt"
#lsetup "views LCG_102b x86_64-centos9-gcc11-opt"
#export PYTHIADIR=/cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/
#export PYTHIA8_DIR=$PYTHIADIR
#export FASTJET=/cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/
