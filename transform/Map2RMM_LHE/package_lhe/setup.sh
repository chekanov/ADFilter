#!/bin/bash


echo "Pythia8 setup"
# run it as: source setup.sh
echo "Set ROOT enviroment for Dijet+Lepton program"

HH=`hostname -A`
echo "HOST=$HH"


export MSOFT=/users/admin/share/sl7

if [[ $HH =~ .*hep.anl.* ]]; then
  echo "HEP ANL ROOT setup"
  MSOFT=/users/admin/share/sl7
  source $MSOFT/setup.sh
  #source $MSOFT/set_asc.sh
fi

if [[ $HH =~ .*lcrc.anl.* ]]; then
  echo "LCRC ANL ROOT setup"
  MSOFT=/soft/hep
  source $MSOFT/hep_setup.sh
fi


if [[ $HH =~ .*atlaslogin* ]]; then
  MSOFT=/users/admin/share/sl7
  source $MSOFT/setup.sh
fi

echo "MSOFT=$MSOFT"

export PYTHIA8=$MSOFT/pythia8
export PYTHIADIR=$PYTHIA8
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
# LHAPDF6 configuration.
export LHAPDF6_USE=true
export LHAPDF6_BIN=$MSOFT/bin/
export PATH=$LHAPDF6_BIN:$PATH
export LHAPDF6_INCLUDE=$MSOFT/lhapdf6/include/
export LHAPDF6_LIB=$MSOFT/lhapdf6/lib

export LD_LIBRARY_PATH=$LHAPDF6_LIB/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH

export HEPMC=$MSOFT/HEPMC
export LD_LIBRARY_PATH=$HEPMC/lib:$LD_LIBRARY_PATH

export PROMC=$MSOFT/promc
export LD_LIBRARY_PATH=$PROMC/lib:$LD_LIBRARY_PATH
export PATH=$PROMC/bin:$PATH
export PKG_CONFIG_PATH=${PROMC}/lib/pkgconfig:${PKG_CONFIG_PATH}
export PYTHONPATH=${PROMC}/python/lib/python2.4/site-packages:$PYTHONPATH
