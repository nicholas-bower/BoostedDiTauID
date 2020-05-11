#!/bin/bash
 WORKDIR="/uscms_data/d3/nbower/FSU/EFilterTest/CMSSW_9_4_13/src/BoostedDiTau/BoostedDiTauReco/"
 export LSB_JOB_REPORT_MAIL=N
    export X509_USER_PROXY=/uscms_data/d3/nbower/x509up_u94443
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /uscms_data/d3/nbower/FSU/EFilterTest/CMSSW_9_4_13/src/BoostedDiTau/BoostedDiTauReco/ && eval `scramv1 runtime -sh`
    cd ${WORKDIR}
    cmsenv
    python /uscms_data/d3/nbower/FSU/EFilterTest/CMSSW_9_4_13/src/BoostedDiTau/BoostedDiTauReco/plotRecoIdEff.py ${1} ${2} ${3}  
