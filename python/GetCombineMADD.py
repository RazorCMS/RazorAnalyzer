from optparse import OptionParser
import ROOT as rt
import sys
import glob
from math import *
import os
from array import *

#local imports 
import rootTools
from framework import Config
from GChiPairs import gchipairs
from SidebandMacro import config
from GetCombine import writeXsecTree
from RunCombine import exec_me

if __name__ == '__main__':

    
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")
    parser.add_option('--signif',dest="doSignificance",default=False,action='store_true',
                  help="for significance instead of limit")
    parser.add_option('--toys',dest="doHybridNew",default=False,action='store_true',
                  help="for toys instead of asymptotic")

    (options,args) = parser.parse_args()

    #get options
    boxInput = options.box
    model = options.model
    directory = options.outDir
    doHybridNew = options.doHybridNew
    doSignificance = options.doSignificance

    #configuration
    cfg = Config.Config(config)
    boxes = boxInput.split('_')
    btag = ''

    haddOutputs = []
    thyXsec = {}
    thyXsecErr = {}
    refXsecFile =  "./data/gluino13TeV.txt" 
    if 'T2' in model:
        refXsecFile = './data/stop13TeV.txt'

    #get theory cross sections and errors
    for mg in range(100,2025,25):
        for line in open(refXsecFile,'r'):
            line = line.replace('\n','')
            if str(mg)==line.split(',')[0]:
                thyXsec[mg] = float(line.split(',')[1]) #pb
                thyXsecErr[mg] = 0.01*float(line.split(',')[2])

    #get combine results
    for mg, mchi in gchipairs(model):
        print "Looking for",mg,mchi
        refXsec = 1.e3*thyXsec[mg]
        modelName = 'SMS-'+('_'.join([model, str(mg), str(mchi)]))
        
        #open file if present
        combineName = 'MADD_'+options.box+'_'+modelName
        if doSignificance and doHybridNew:
            prefix = 'higgsCombineSignif'
            suffix = 'HybridNew'
        elif doHybridNew: 
            prefix = 'higgsCombineToys'
            suffix = 'HybridNew'
        elif doSignificance:
            prefix = 'higgsCombine'
            suffix = 'ProfileLikelihood'
        else:
            prefix = 'higgsCombine'
            suffix = 'Asymptotic'
        filename = directory+'/'+prefix+combineName+'.'+suffix+'.mH120.root'
        if not glob.glob(filename): 
            print "Didn't find file",filename
            continue
        tFile = rt.TFile.Open(filename)
        print "File:",filename;

        #check file and tree within
        try:
            if tFile.InheritsFrom("TFile") is False:
                print "Isn't a TFile"
                continue
        except:
            print "File not found!"
            continue
        limit = tFile.Get("limit")
        try:
            if limit.InheritsFrom("TTree") is False: 
                tFile.cd()
                tFile.Close()
                print "'limit' doesn't appear to be a TTree"
                continue
        except:
            tFile.cd()
            tFile.Close()
            print "Error getting tree from file"
            continue
        if doSignificance and limit.GetEntries() < 1: 
            tFile.cd()
            tFile.Close()
            print "Not enough entries in tree!"
            continue
        if (not doSignificance) and limit.GetEntries() < 6: 
            tFile.cd()
            tFile.Close()
            print "Not enough entries in tree!"
            continue

        #get limits out
        limit.Draw('>>elist','','entrylist')
        elist = rt.gDirectory.Get('elist')
        entry = elist.Next()
        limit.GetEntry(entry)
        limits = []
        while True:
            if entry == -1: break
            limit.GetEntry(entry)
            if doSignificance:
                limits.append(max(0.0,limit.limit))
            else:
                limits.append(refXsec*(1.e-3)*limit.limit)
            entry = elist.Next()
        tFile.cd()
        tFile.Close()
            
        limits.reverse()
        print limits
        
        haddOutput = writeXsecTree(boxInput, model, directory, mg, mchi, [limits[0]],[limits[1]],[limits[2]],[limits[3]],[limits[4]],[limits[5]])
        haddOutputs.append(haddOutput)

    if doHybridNew:
        os.system("hadd -f %s/xsecUL_HybridNew_%s.root %s"%(directory,boxInput," ".join(haddOutputs)))
    else:
        os.system("hadd -f %s/xsecUL_Asymptotic_%s.root %s"%(directory,boxInput," ".join(haddOutputs)))
