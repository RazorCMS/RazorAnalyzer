import ROOT as rt
import sys
import rootTools
import glob
from math import *
import os
from array import *
import random
from optparse import OptionParser
from framework import Config


def testWorkspace(w,outFile,singleBox):
     
    w.Print("v")
    CMS_th1x = w.var("CMS_th1x")
    CMS_th1x.Print("v")
    CMS_channel = w.cat("CMS_channel")
    
    CMS_set = rt.RooArgSet()
    CMS_set.add(CMS_channel)
    CMS_set.add(CMS_th1x)
    
    data_obs = w.data("data_obs")

    model_b = w.pdf("model_b")

    signalNuis = ["lumi"]
            
    for varName in signalNuis:
        print varName
        w.var(varName).setConstant()

    allParams = model_b.getParameters(data_obs)
    rt.RooStats.RemoveConstantParameters(allParams)

    nll = model_b.createNLL(data_obs,rt.RooFit.CloneData(False),rt.RooFit.Constrain(allParams))

    minim = rt.RooMinimizer(nll)

    strategy = rt.Math.MinimizerOptions.DefaultStrategy()

    minim.setStrategy(strategy)
    tol = rt.Math.MinimizerOptions.DefaultTolerance()
    tol = max(tol,1.0)
    minim.setEps(tol)

    minim.optimizeConst(2)

    minimizer = rt.Math.MinimizerOptions.DefaultMinimizerType()
    algorithm = rt.Math.MinimizerOptions.DefaultMinimizerAlgo()

    status = minim.minimize(minimizer,algorithm)
    
    fr_b = minim.save()
    nll_b = fr_b.minNll()
    fr_b.Print("v")
    
    print "background-only nll = %f"%nll_b
    minim.optimizeConst(False)

    data_obs.Print("v")
    
        
    argList = fr_b.floatParsFinal()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-o','--output',dest="output",default="./test.root",type="string",
                  help="Output file to store results")

    (options,args) = parser.parse_args()

     
    rt.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    
    for f in args:
    if f.lower().endswith('.root'):
        workspaceFileName = f
    workspaceFile = rt.TFile.Open(workspaceFileName,"READ")
    outFile = rt.TFile.Open(options.output,"RECREATE")
    w = workspaceFile.Get("w")
    singleBox = options.box

    print "SINGLE BOX %s"%singleBox 
    testWorkspace(w,outFile,singleBox)
