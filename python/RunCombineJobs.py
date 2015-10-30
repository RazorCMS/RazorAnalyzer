import time
from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob
from GChiPairs import gchipairs
    
def writeBashScript(box,btag,model,mg,mchi,lumi,config,submitDir,isData,fit,penalty,inputFitFile):
    
    massPoint = "%i_%i"%(mg, mchi)
    dataString = ''
    if isData:
        dataString = '--data'

    fitString = ''
    if fit:
        fitString = '--fit'
        
    penaltyString = ''
    if penalty:
        penaltyString = '--penalty'
        
    # prepare the script to run
    outputname = submitDir+"/submit_"+model+"_"+massPoint+"_lumi-%.3f_"%(lumi)+btag+"_"+box+".src"
        
    ffDir = submitDir+"/logs_"+model+"_"+massPoint+"_"+btag+"_"+box
    user = os.environ['USER']
    pwd = os.environ['PWD']
    
    combineDir = "/afs/cern.ch/work/%s/%s/RAZORRUN2/Limits/"%(user[0],user)

    script =  '#!/usr/bin/env bash -x\n'
    script += 'mkdir -p %s\n'%combineDir
    
    script += 'echo $SHELL\n'
    script += 'pwd\n'
    script += 'cd /afs/cern.ch/work/%s/%s/RAZORRUN2/CMSSW_7_1_5/src/RazorAnalyzer \n'%(user[0],user)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc481\n"
    script += "export CMSSW_BASE=/afs/cern.ch/work/%s/%s/RAZORRUN2/CMSSW_7_1_5\n"%(user[0],user)
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'cd - \n'
    script += "export TWD=${PWD}/%s_%s_lumi-%.3f_%s_%s\n"%(model,massPoint,lumi,btag,box)
    script += "mkdir -p $TWD\n"
    script += "cd $TWD\n"
    script += 'pwd\n'
    script += 'git clone git@github.com:RazorCMS/RazorAnalyzer\n'
    script += 'cd RazorAnalyzer\n'
    script += 'source setup.sh\n'
    script += 'make\n'
    script += 'mkdir -p Datasets\n'
    script += 'mkdir -p %s\n'%submitDir
    script += 'python python/RunCombine.py -i %s --mGluino %i --mLSP %i %s -c %s --lumi-array %f -d %s -b %s %s %s\n'%(inputFitFile,mg,mchi,dataString,config,lumi,submitDir,box,fitString,penaltyString)
    script += 'cp %s/higgsCombine* %s/\n'%(submitDir,combineDir) 
    script += 'cd ../..\n'
    script += 'rm -rf $TWD\n'
    
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
    parser.add_option('-l','--lumi',dest="lumi", default=0.210,type="float",
                  help="lumi in fb^-1, e.g.: 0.210")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--fit',dest="fit",default=False,action='store_true',
                  help="Turn on pre-fit")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes for data")
    parser.add_option('--penalty',dest="penalty",default=False,action='store_true',
                  help="penalty terms on background parameters")
    parser.add_option('--no-sub',dest="noSub", default=False,action='store_true',
                  help="no submission")
    parser.add_option('-q','--queue',dest="queue",default="1nh",type="string",
                  help="queue: 1nh, 8nh, 1nd, etc.")
    parser.add_option('--mg-geq',dest="mgMin",default=-1,type="float",
                  help="mgMin ")
    parser.add_option('--mg-lt',dest="mgMax",default=10000,type="float",
                  help="mgMax ")
    parser.add_option('--mchi-geq',dest="mchiMin",default=-1,type="float",
                  help="mchiMin ")
    parser.add_option('--mchi-lt',dest="mchiMax",default=10000,type="float",
                  help="mchiMax ")
    parser.add_option('--done-file',dest="doneFile",default=None,type="string",
                  help="file containing output files")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default='FitResults/BinnedFitResults.root',type="string",
                  help="input fit file")

    (options,args) = parser.parse_args()


    cfg = Config.Config(options.config)

    boxes = options.box.split('_')

    btag = ''
    for box in boxes:            
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
        btagMin = z[0]
        btagMax = z[-1]
        if btagMax-1>btagMin:          
            btag = '%i-%ibtag'%(btagMin,btagMax-1)
        else:
            btag = '%ibtag'%(btagMin)
                
    nJobs = 0
    donePairs = []
    if options.doneFile is not None:
        with open(options.doneFile,'r') as f:            
            allFiles = [ line.replace('\n','') for line in f.readlines()]
            for (mg, mchi) in gchipairs(options.model):
                outputname = 'higgsCombine%s_%i_%i_lumi-%.3f_%s_%s.Asymptotic.mH120.root'%(options.model,mg,mchi,options.lumi,btag,options.box)
                if outputname in allFiles: donePairs.append((mg,mchi))

                    
    for (mg, mchi) in gchipairs(options.model):
        if not (mg >= options.mgMin and mg < options.mgMax): continue
        if not (mchi >= options.mchiMin and mchi < options.mchiMax): continue
        if (mg, mchi) in donePairs: continue
        nJobs+=1
        outputname,ffDir = writeBashScript(options.box,btag,options.model,mg,mchi,options.lumi,options.config,options.outDir,options.isData,options.fit,options.penalty,options.inputFitFile)
        
        pwd = os.environ['PWD']
        os.system("mkdir -p "+pwd+"/"+ffDir)
        os.system("echo bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)        
        #os.system("source "+pwd+"/"+outputname)
        if not options.noSub:
            time.sleep(3)
            os.system("bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)

    print "nJobs = %i"%nJobs


