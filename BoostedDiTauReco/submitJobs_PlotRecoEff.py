from __future__ import print_function, division

import os
import sys
# import numpy as np
import argparse
import ROOT
inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs to produce efficiency plots.')
inputArgumentsParser.add_argument('--queue', default="8nh", help="queue to submit to for each set of events", type=str)
inputArgumentsParser.add_argument('--process', default = 'var', help="Type or process being carried out, var or eff", type=str)
inputArgumentsParser.add_argument('--sample', default = 'DYJetsToLL', help="The type of sample being run on (DY, QCD, TCP etc...)", type=str)
inputArgumentsParser.add_argument('--isDryRun', default =False, help="Do not submit the actual jobs; instead print what would have been submitted.", type = bool)

inputArguments = inputArgumentsParser.parse_args()

print(" >> Submitting jobs for running event selection...")

# Load input TTrees into TChain
sample = inputArguments.sample
process = inputArguments.process
inFileDir='/uscms_data/d3/nbower/FSU/EFilterTest/CMSSW_9_4_13/src/BoostedDiTau/BoostedDiTauReco/filelists/DYJetsToLL_94X/'
n = 1
conFile = open("condSub.jdl","w") 
conFile.write("request_memory = 4200\n")
conFile.write("Executable = run_condor.csh\n")  
conFile.write("Should_Transfer_Files = YES\n")  
conFile.write("WhenToTransferOutput = ON_EXIT\n")  
conFile.write("Transfer_Input_Files = run_condor.csh\n")  
conFile.write("Output = ./condorOut/condor_$ENV(JOBNUMBER).stdout\n")
conFile.write("/n")  

for subdir, dirs, files in os.walk(inFileDir):
    for file in files:
        if n !=1:continue
        n+=1
        commandToCall = "bsub -q {queue} job_PlotRecoIDEff.sh {file} {process} {sample}".format(queue=inputArguments.queue,  file=file, sample=sample, process = process)
        print ("Generated command: " + commandToCall)      
        if (inputArguments.isDryRun): 
            print("Not submitting due to dryRun flag.")
        else:
            os.system(commandToCall)   
            print ("Submitted.")     
print ("Final submitted, wish me luck.")
