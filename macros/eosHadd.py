import subprocess
import sys,string,math,os
import ConfigParser
import glob
import numpy as np


filesPerList=50000


def checkAndMakeDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def clearDir(dir):
    for fil in glob.glob(dir+"/*"):
        os.remove(fil)

prefix= "root://cmseos.fnal.gov./"
rootFileDir="/store/user/nbower/QCD_eff_plots/"
outFileName= './QCD_Eff_6_2.root'
searchString='/Neutrino_E-10_gun/RunIISummer19ULPrePremix-UL17_106X_mc2017_realistic_v6-v1/PREMIX'
os.system('xrdfsls' + rootFileDir)
#query = 'eosls '+ rootFileDir
query = 'xrdfs ' + prefix + " ls " + rootFileDir
os.system(query)
files=os.popen(query).read().split()
#out = open("runFi.csh", "w")
#out.write('hadd '+ outFileName + ' ')
runString = 'hadd -f -k '+ outFileName + ' '
for nf in range(1, len(files)+1):
    filelistIdx=int((nf-1))
    #out.write(prefix+files[nf-1]+" ")
    runString +=prefix+files[nf-1]+" "
print 'Running: \n'+runString
os.system(runString)
