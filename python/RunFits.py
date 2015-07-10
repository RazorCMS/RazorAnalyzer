from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import glob
from RunCombine import exec_me

if __name__ == '__main__':


    #boxes = ['MultiJet','FourJet','SixJet','DiJet']
    #configs = ['config/run2_rsqturnoff_0btag.config','config/run2_0btag.config','config/run2_rsqturnoff_1btag.config','config/run2_1btag.config']

    #boxes = ['MuEle','EleEle','MuMu','MuMultiJet','MuJet','EleMultiJet','EleJet','LooseLeptonMultiJet','MultiJet','DiJet']
    boxes = ['MuMultiJet','MuJet','EleMultiJet','EleJet','LooseLeptonMultiJet','MultiJet','DiJet']
    boxes.extend(['MuFourJet','MuSixJet','EleFourJet','EleSixJet','SixJet','FourJet','LooseLeptonSixJet','LooseLeptonFourJet'])

    
    boxes = ['MultiJet','MuMultiJet','EleMultiJet','LooseLeptonMultiJet','DiJet']
    configs = ['config/run2_1-3btag.config']
    
    
    lumi = 3000
    
    for box in boxes:
        for cfg in configs:
            if box in ['MuEle','EleEle','MuMu']:
                btag = '0-1btag'
            elif '1-3btag' in cfg:
                btag = '1-3btag'
            else:
                btag = '0-3btag'
                
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root'%(lumi/1000.,btag,box)
            if not glob.glob(backgroundDsName):                
                exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w -q Backgrounds/*.root'%(cfg,box),False)
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.1f_%s_%s.root'%(lumi/1000.,btag,box)
            if not glob.glob(backgroundDsName):
                exec_me('python python/RooDataSet2UnweightedDataset.py -c %s -b %s -d Datasets/ Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root'%(cfg,box,lumi/1000.,btag,box),False)
            
            exec_me('mkdir -p fits_07102015/weighted/%s_%s/' %(box,btag),False)
            exec_me('mkdir -p fits_07102015/unweighted/%s_%s/' %(box,btag),False)
            #  --fit-region LowMR,LowRsq 
            #exec_me('python python/Fit.py -c %s -d fits_06212015/unbinned/%s_%s/ --fit-region LowMR,LowRsq -b %s Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.1f_%s_%s.root' %(cfg,box,btag,box,lumi/1000,btag,box),False)

            #exec_me('python python/BinnedFit.py -c %s -d fits_07102015/weighted/%s_%s/ -b %s Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root' %(cfg,box,btag,box,lumi/1000,btag,box),False)
            #exec_me('python python/RunToys.py -c %s -d fits_07102015/weighted/%s_%s/ -b %s -i fits_07102015/weighted/%s_%s/BinnedFitResults_%s.root -t 10000' %(cfg,box,btag,box,box,btag,box),False)
            exec_me('python python/PlotFit.py -c %s -d fits_07102015/weighted/%s_%s/ -b %s -i fits_07102015/weighted/%s_%s/BinnedFitResults_%s.root -t fits_07102015/weighted/%s_%s/toys_Bayes_%s.root --print-errors' %(cfg,box,btag,box,box,btag,box,box,btag,box),False)
                        
            #exec_me('python python/BinnedFit.py -c %s -d fits_07102015/unweighted/%s_%s/ -b %s Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.1f_%s_%s.root' %(cfg,box,btag,box,lumi/1000,btag,box),False)
            #exec_me('python python/RunToys.py -c %s -d fits_07102015/unweighted/%s_%s/ -b %s -i fits_07102015/unweighted/%s_%s/BinnedFitResults_%s.root -t 10000' %(cfg,box,btag,box,box,btag,box),False)
            exec_me('python python/PlotFit.py -c %s -d fits_07102015/unweighted/%s_%s/ -b %s -i fits_07102015/unweighted/%s_%s/BinnedFitResults_%s.root -t fits_07102015/unweighted/%s_%s/toys_Bayes_%s.root --print-errors' %(cfg,box,btag,box,box,btag,box,box,btag,box),False)
                      
            
            

