from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import glob
from RunCombine import exec_me
import datetime

if __name__ == '__main__':
    
    boxes = ['MultiJet','MuMultiJet','EleMultiJet','LooseLeptonMultiJet','DiJet']
    fits = ['Sideband','Full']
    configs = ['config/run2.config']
    
    lumi = 3000
    dryRun = False
    btag = '0-3btag'
    
    for box in boxes:
        for cfg in configs:
                
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root'%(lumi/1000.,btag,box)
            exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w -q Backgrounds/*.root'%(cfg,box),dryRun)
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.1f_%s_%s.root'%(lumi/1000.,btag,box)
            exec_me('python python/RooDataSet2UnweightedDataset.py -c %s -b %s -d Datasets/ Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root'%(cfg,box,lumi/1000.,btag,box),dryRun)
            lumiString = '%.1f'%(lumi/1000)
            lumiString = lumiString.replace('.','p')

            for weight in ["weighted","unweighted"]:
                backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_%s_lumi-%.1f_%s_%s.root'%(weight,lumi/1000.,btag,box)
                for fit in fits:                    
                    
                    fitString = ''
                    if fit=='Sideband':
                        fitString = '--fit-region LowMR,LowRsq'

                    dateString = str(datetime.date.today())
                    
                    outDir = "fits_%s/%s_lumi-%s/%s_%s_%s/"%(dateString,weight,lumiString,box,btag,fit)
                    exec_me('mkdir -p %s' %(outDir),dryRun)
                    exec_me('python python/BinnedFit.py -c %s -d %s/ -b %s -l %f %s %s' %(cfg,outDir,box,lumi,backgroundDsName,fitString),dryRun)
                    #exec_me('python python/RunToys.py -c %s -d %s/ -b %s -l %f -i %s/BinnedFitResults_%s.root -t 10000' %(cfg,outDir,box,lumi,outDir,box),dryRun)
                    #exec_me('python python/PlotFit.py -c %s -d %s/ -b %s -l %f -i %s/BinnedFitResults_%s.root -t %s/toys_Bayes_%s.root %s' %(cfg,outdir,box,lumi,outDir,box,outDir,box,fitString),dryRun)
                    exec_me('python python/PlotFit.py -c %s -d %s/ -b %s -l %f -i %s/BinnedFitResults_%s.root %s' %(cfg,outDir,box,lumi,outDir,box,fitString),dryRun)
                      
            
            

