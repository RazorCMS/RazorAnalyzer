from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob

def exec_me(command,dryRun=True):
    print command
    if not dryRun: os.system(command)

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('--mGluino',dest="mGluino", default=-1,type="float",
                  help="mass of gluino")
    parser.add_option('--mStop',dest="mStop", default=-1,type="float",
                  help="mass of stop")
    parser.add_option('--mLSP',dest="mLSP", default=-1,type="float",
                  help="mass of LSP")
    parser.add_option('--lumi-array',dest="lumi_array", default="0.2,0.5,1,3,4,7,10",type="string",
                  help="lumi array in fb^-1, e.g.: 0.2,0.5,1,3,4,7,10")
    parser.add_option('--signif',dest="signif",default=False,action='store_true',
                  help="calculate significance instead of limit")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--fit',dest="fit",default=False,action='store_true',
                  help="Turn on pre-fit")
    parser.add_option('--min-tol',dest="min_tol",default=0.001,type="float",
                  help="minimizer tolerance (default = 0.001)")
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")

    (options,args) = parser.parse_args()

    boxes = options.box.split('_')

    signif = options.signif
    
    lumiArray = array('d',[float(lumi) for lumi in options.lumi_array.split(',')])
    cfg = Config.Config(options.config)
    
    model = options.model
    if options.mGluino>-1:
        massPoint = 'mGl-%i_mLSP-%i'%(options.mGluino,options.mLSP)
    elif options.mStop>-1:
        massPoint = 'mStop-%i_mLSP-%i'%(options.mStop,options.mLSP)


    fit = ''
    if options.fit:
        fit = '--fit'
        
    for lumi in lumiArray:
        for box in boxes:            
            z = array('d', cfg.getBinning(box)[2]) # nBtag binning
            btagMin = z[0]
            btagMax = z[-1]
            if btagMax-1>btagMin:          
                btag = '%i-%ibtag'%(btagMin,btagMax-1)
            else:
                btag = '%ibtag'%(btagMin)    
            #signalDsName = 'Datasets/RazorInclusive_SMS-%s_2J_%s_weighted_lumi-%.1f_%s_%s.root'%(model,massPoint,lumi,btag,box)
            #backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.1f_%s_%s.root'%(lumi,btag,box)
            signalDsName = 'Datasets/RazorInclusive_SMS-%s_2J_%s_weighted_lumi-%.1f_%s_%s.root'%(model,massPoint,3.0,btag,box)
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root'%(3.0,btag,box)
            if not glob.glob(signalDsName):
                exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w Signals/RazorInclusive_SMS-%s_2J_%s_*.root'%(options.config,box,model,massPoint),options.dryRun)
            if not glob.glob(backgroundDsName):                
                exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w -q Backgrounds/*.root'%(options.config,box),options.dryRun)
                #exec_me('python python/DustinTuple2RooDataSet.py -c %s -b %s -l %f -d Datasets/ -w -q Backgrounds/*.root'%(options.config,box,1000.*lumi),options.dryRun)
                #exec_me('python python/RooDataSet2UnweightedDataSet.py -c %s -b %s -d Datasets/ Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root'%(options.config,box,lumi,btag,box),options.dryRun)
            exec_me('python python/WriteDataCard.py -l %f -c %s -b %s -d %s %s %s %s'%(1000*lumi,options.config,box,options.outDir,fit,signalDsName,backgroundDsName),options.dryRun)
            
            if signif:
                exec_me('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 --toysFreq %s/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.1f_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)
            else:
                exec_me('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt -n %s_%s_lumi-%.1f_%s_%s --minimizerTolerance %f'%(options.outDir,model,massPoint,lumi,btag,box,model,massPoint,lumi,btag,box,options.min_tol),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.1f_%s_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,btag,box,options.outDir),options.dryRun)
        if len(boxes)>1:
            for box in boxes: exec_me('cp %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt .'%(options.outDir,model,massPoint,lumi,btag,box),options.dryRun)
            cmds = ['%s=razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(box,model,massPoint,lumi,btag,box) for box in boxes]
            exec_me('combineCards.py %s > %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(' '.join(cmds),options.outDir,model,massPoint,lumi,btag,options.box),options.dryRun)
            exec_me('cat %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(options.outDir,model,massPoint,lumi,btag,options.box),options.dryRun)
            if signif:
                exec_me('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 --toysFreq %s/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s_%s'%(options.outDir,model,massPoint,lumi,options.box,model,massPoint,lumi,btag,options.box),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.1f_%s_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,btag,options.box,options.outDir),options.dryRun)
            else:
                exec_me('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt -n %s_%s_lumi-%.1f_%s_%s --minimizerTolerance %f'%(options.outDir,model,massPoint,lumi,btag,options.box,model,massPoint,lumi,btag,options.box,options.min_tol),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%.1f_%s_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,btag,options.box,options.outDir),options.dryRun)
            for box in boxes: exec_me('rm razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(model,massPoint,lumi,btag,box),options.dryRun)
 
