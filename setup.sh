#!/bin/bash

echo "Setting up ROOT"
export SCRAM_ARCH=slc6_amd64_gcc530
export CMSSW_VERSION=CMSSW_9_2_4
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src
eval `scramv1 runtime -sh`
cd -

echo 'which root'
which root

#eof
