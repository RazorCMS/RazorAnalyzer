from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import glob
from RunCombine import exec_me
import datetime

if __name__ == '__main__':
    
    boxes = ['MultiJet','MuMultiJet','EleMultiJet','LooseLeptonMultiJet']
    fits = ['Full','Sideband']
    weights = ['weighted']
    configs = ['config/run2_sideband.config']
    
    lumi = 32000
    dryRun = True
    btag = '0-3btag'
    
    for box in boxes:
        for cfg in configs:                
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.3f_%s_%s.root'%(lumi/1000.,btag,box)
            exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ Backgrounds/*.root -l %f -w'%(cfg,box,lumi),dryRun)
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.3f_%s_%s.root'%(lumi/1000.,btag,box)
            exec_me('python python/RooDataSet2UnweightedDataset.py -c %s -b %s -d Datasets/ Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.3f_%s_%s.root'%(cfg,box,lumi/1000.,btag,box),dryRun)
            for weight in weights:
                backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_%s_lumi-%.3f_%s_%s.root'%(weight,lumi/1000.,btag,box)
                for fit in fits:                     
                    fitString = ''
                    if fit=='Sideband':
                        fitString = '--fit-region LowMR,LowRsq'

                    dateString = str(datetime.date.today()).replace("-","_")
                                       
                    outDir = "fits_%s/%s_%iifb/%s/"%(dateString,box,lumi/1000.,fit)

                    exec_me('python python/BinnedFit.py -c %s -d %s -b %s -l %f %s %s' %(cfg,outDir,box,lumi,backgroundDsName,fitString),dryRun)
                    exec_me('python python/RunToys.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t 10000' %(cfg,outDir,box,lumi,outDir,box),dryRun)
                    exec_me('python python/PlotFit.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root -t %s/toys_Bayes_%s.root %s' %(cfg,outDir,box,lumi,outDir,box,outDir,box,fitString),dryRun)
                    #exec_me('python python/PlotFit.py -c %s -d %s -b %s -l %f -i %s/BinnedFitResults_%s.root %s' %(cfg,outDir,box,lumi,outDir,box,fitString),dryRun)

                if len(fits)==2:
                    dateString = str(datetime.date.today()).replace("-","_")                    
                    outDir1 = "fits_%s/%s_%iifb/%s/"%(dateString,box,lumi/1000.,fits[0])
                    outDir2 = "fits_%s/%s_%iifb/%s/"%(dateString,box,lumi/1000.,fits[1])
                    outDir = "fits_%s/%s_%iifb/"%(dateString,box,lumi/1000.)
                    #exec_me('python python/CompareFits.py -b %s -c %s -l %f -1 %s/BinnedFitResults_%s.root -2 %s/BinnedFitResults_%s.root -d %s'%(box,cfg,lumi,outDir1,box,outDir2,box,outDir),dryRun)
            
            

