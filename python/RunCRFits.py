from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import glob
from RunCombine import exec_me

if __name__ == '__main__':


    boxes = ['TTBarSingleLepton','WSingleLepton']
    fits = ['Sideband','Full']
    configs = ['config/controlsample.config']
    
    
    lumi = 23.8

    dryRun=False
    
    for box in boxes:
        for cfg in configs:            

            exec_me('python python/CRTuple2RooDataSet.py -c %s -b %s -d ControlSampleFits/ ControlSampleFits/RunTwoRazorControlRegions_OneLeptonFull_Run2015B_GoodLumiDCS_NoDuplicates.root -l %f'%(cfg,box,lumi),dryRun)
            
            for fit in fits:                
                if box in ['TTBarSingleLepton']:
                    btag = '1btag'
                else:
                    btag = '0btag'
                
                lumiString = '%.4f'%(lumi/1000)
                #lumiString = lumiString.replace('.','p')

                fitString = ''
                if fit=='Sideband':
                    fitString = '--fit-region LowMR,LowRsq'
                    
                dsName = 'ControlSampleFits/RunTwoRazorControlRegions_OneLeptonFull_Run2015B_GoodLumiDCS_NoDuplicates_lumi-%s_%s_%s.root'%(lumiString,btag,box)
                exec_me('python python/BinnedFit.py -c %s -b %s -l %f -d ControlSampleFits/%s/%s/ %s %s --data' %(cfg,box,lumi,box,fit,fitString,dsName),dryRun)
                #exec_me('python python/RunToys.py -c %s -b %s -i ControlSampleFits/%s/%s/BinnedFitResults_%s.root -d ControlSampleFits/%s/%s/ -t 3000'%(cfg,box,box,fit,box,box,fit),dryRun)
                #exec_me('python python/PlotFit.py -c %s -b %s -l %f -i ControlSampleFits/%s/%s/BinnedFitResults_%s.root -d ControlSampleFits/%s/%s/ -t ControlSampleFits/%s/%s/toys_Bayes_%s.root --data' %(cfg,box,lumi,box,fit,box,box,fit,box,fit,box),dryRun)
                exec_me('python python/PlotFit.py -c %s -b %s -l %f -i ControlSampleFits/%s/%s/BinnedFitResults_%s.root -d ControlSampleFits/%s/%s/ %s --data' %(cfg,box,lumi,box,fit,box,box,fit,fitString),dryRun)
                      
            
            
