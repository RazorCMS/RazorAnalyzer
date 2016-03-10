import time
from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import glob

#local imports
from framework import Config
from GChiPairs import gchipairs
from WriteRazorMADDCard import LUMI
    
def writeBashScript(box,model,mg,mchi,submitDir,noSys,fitSys,signif=False,contamination=False):
    
    massPoint = "%i_%i"%(mg, mchi)

    sysString = ''
    if noSys:
        sysString += '--no-sys --no-stat'
    if fitSys:
        sysString += '--fit-sys'
    sigString = ''
    if signif:
        sigString = '--signif'
    contamString = ''
    if contamination:
        contamString = '--contamination'

    particleString = '--mGluino'
    if 'T2' in model:
        particleString = '--mStop'

    # prepare the script to run
    outputname = 'Limits/'+submitDir+"/submit_"+model+"_"+massPoint+"_lumi-%.3f_"%(LUMI*1.0/1000)+box+".src"
        
    ffDir = 'Limits/'+submitDir+"/logs_"+model+"_"+massPoint+"_"+box
    user = os.environ['USER']
    pwd = os.environ['PWD']
    cmsswBase = os.environ['CMSSW_BASE']
    
    combineDir = "/afs/cern.ch/work/%s/%s/RAZORRUN2/Limits/%s/%s/"%(user[0],user,submitDir,model) # directory where combine output files will be copied

    script =  '#!/usr/bin/env bash -x\n'
    script += 'mkdir -p %s\n'%combineDir
    
    script += 'echo $SHELL\n'
    script += 'pwd\n'
    script += 'cd %s/src/RazorAnalyzer \n'%(cmsswBase)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc481\n"
    script += "export CMSSW_BASE=%s\n"%(cmsswBase)
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'cd - \n'
    script += "export TWD=${PWD}/%s_%s_lumi-%.3f_%s\n"%(model,massPoint,LUMI*1.0/1000,box)
    script += "mkdir -p $TWD\n"
    script += "cd $TWD\n"
    script += 'pwd\n'
    script += 'git clone https://github.com/RazorCMS/RazorAnalyzer.git\n'
    script += 'cd RazorAnalyzer\n'
    script += 'git checkout -b Limits LimitsMADD20160310\n' 
    script += 'make lxplus\n'
    script += 'mkdir -p %s\n'%submitDir
    script += 'python python/WriteRazorMADDCard.py --model %s %s %i --mLSP %i --dir %s --box %s %s %s %s\n'%(model,particleString,mg,mchi,submitDir,box,sysString,sigString,contamString)
    script += 'cp %s/higgsCombine* %s/\n'%(submitDir,combineDir) 
    script += 'cd ../..\n'
    script += 'rm -rf $TWD\n'
    
    outputfile = open(outputname,'w')
    outputfile.write(script)
    outputfile.close

    return outputname,ffDir



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
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
    parser.add_option('--no-sys',dest="noSys",default=False,action='store_true',
                  help="no shape systematic uncertainties")
    parser.add_option('--fit-sys',dest="fitSys",default=False,action='store_true',
                  help="use fit vs MC systematic")
    parser.add_option('--signif',dest="signif",default=False,action='store_true',
                  help="Compute significance instead of limit")
    parser.add_option('--contamination',dest="contamination",default=False,action='store_true',
                  help="Propagate uncertainty on signal contamination")

    (options,args) = parser.parse_args()

    boxes = options.box.split('_')

    nJobs = 0
    donePairs = []
    if options.doneFile is not None:
        if options.signif:
            combineMethod = 'ProfileLikelihood'
        else:
            combineMethod = 'Asymptotic'
        with open(options.doneFile,'r') as f:            
            allFiles = [ line.replace('\n','') for line in f.readlines()]
            for (mg, mchi) in gchipairs(options.model):
                outputname = 'higgsCombineMADD_%s_SMS-%s_%i_%i.%s.mH120.root'%(options.box,options.model,mg,mchi,combineMethod)
                if outputname in allFiles: donePairs.append((mg,mchi))

    thyXsec = {}
    if "T1" in options.model or "T5" in options.model:
        xsecFile = 'data/gluino13TeV.txt'
    if "T2" in options.model:
        xsecFile = 'data/stop13TeV.txt'
        
    for line in open(xsecFile,'r'):
        for (mg, mchi) in gchipairs(options.model):
            if str(int(mg))==line.split(',')[0]:
                thyXsec[(mg,mchi)] = float(line.split(',')[1]) #pb

    for (mg, mchi) in gchipairs(options.model):
        if not (mg >= options.mgMin and mg < options.mgMax): continue
        if not (mchi >= options.mchiMin and mchi < options.mchiMax): continue
        if (mg, mchi) in donePairs: 
            print (mg,mchi),"is already done; skipping"
            continue
        nJobs+=1

        pwd = os.environ['PWD']
        os.system("mkdir -p "+pwd+"/Limits/"+options.outDir)
        outputname,ffDir = writeBashScript(options.box, options.model, mg, mchi, 
                options.outDir, options.noSys, options.fitSys, options.signif, options.contamination)
        
        os.system("mkdir -p "+pwd+"/"+ffDir)
        os.system("echo bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)        
        if not options.noSub:
            time.sleep(3)
            os.system("bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)

    print "nJobs = %i"%nJobs
