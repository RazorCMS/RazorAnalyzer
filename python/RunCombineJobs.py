from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob
from GChiPairs import gchipairs
    
def writeBashScript(box,btag,model,mg,mchi,lumi,config,submitDir,isData,fit):
    
    massPoint = "%i_%i"%(mg, mchi)
    dataString = ''
    if isData:
        dataString = '--data'

    fitString = ''
    if fit:
        fitString = '--fit'
        
    # prepare the script to run
    outputname = submitDir+"/submit_"+model+"_"+massPoint+"_lumi-%.3f_"%(lumi)+btag+"_"+box+".src"
        
    ffDir = submitDir+"/logs_"+model+"_"+massPoint+"_"+btag+"_"+box
    user = os.environ['USER']
    pwd = os.environ['PWD']
    
    combineDir = "/afs/cern.ch/work/%s/%s/RAZORRUN2/CMSSW_7_1_5/src/RazorAnalyzer/cards_data/"%(user[0],user)

    script =  '#!/usr/bin/env bash -x\n'
    script += 'mkdir -p %s\n'
    script += 'echo $SHELL\n'
    script += 'pwd\n'
    script += 'cd /afs/cern.ch/work/%s/%s/RAZORRUN2/CMSSW_7_1_5/src/RazorAnalyzer \n'%(user[0],user)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc481\n"
    script += "export CMSSW_BASE=/afs/cern.ch/work/%s/%s/RAZORRUN2/CMSSW_7_1_5\n"%(user[0],user)
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'source setup.sh\n'
    script += 'python python/RunCombine.py --mGluino %i --mLSP %i %s -c %s --lumi-array %f -d %s -b %s %s'%(mg,mchi,dataString,config,lumi,submitDir,box,fitString)
    
    outputfile = open(outputname,'w')
    outputfile.write(script)
    outputfile.close

    return outputname,ffDir



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('--lumi',dest="lumi", default=0.210,type="float",
                  help="lumi in fb^-1, e.g.: 0.210")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--fit',dest="fit",default=False,action='store_true',
                  help="Turn on pre-fit")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes for data")
    parser.add_option('--no-sub',dest="noSub", default=False,action='store_true',
                  help="no submission")
    parser.add_option('-q','--queue',dest="queue",default="1nh",type="string",
                  help="queue: 1nh, 8nh, 1nd, etc.")


    (options,args) = parser.parse_args()


    btag = '0-3btag'

    for (mg, mchi) in gchipairs(options.model):
        
        outputname,ffDir = writeBashScript(options.box,btag,options.model,mg,mchi,options.lumi,options.config,options.outDir,options.isData,options.fit)
        
        pwd = os.environ['PWD']
        os.system("echo bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)
        if not options.noSub:
            time.sleep(3)
            os.system("bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)


