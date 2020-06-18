#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh   ## if a bash script, use .sh instead of .csh
setenv SCRAM_ARCH slc7_amd64_gcc630  
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/nbower/CMSSW.tgz ./CMSSW.tgz
tar -xf CMSSW.tgz
rm CMSSW.tgz
#setenv SCRAM_ARCH slc6_amd64_gcc530
cd CMSSW_9_4_13/src/
source /cvmfs/cms.cern.ch/cmsset_default.csh 
scramv1 b ProjectRename
cd ./BoostedDiTauID/
ls
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
#cd BoostedDiTau/BoostedDiTauReco/
echo "Arguments passed to this script are: for 1: $1"
python ${1} ${2} ${3} ${4}
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_9_4_13
