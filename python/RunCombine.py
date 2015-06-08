from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use"),
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
    parser.add_option('--lumi-array',dest="lumi_array", default="0.2,3,4,10",type="string",
                  help="lumi array in fb^-1, e.g.: 0.2,3,4,10")
    parser.add_option('--signif',dest="signif",default=False,action='store_true',
                  help="calculate significance instead of limit")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--fit',dest="fit",default=False,action='store_true',
                  help="Turn on pre-fit")


    (options,args) = parser.parse_args()

    boxes = options.box.split('_')

    signif = options.signif
    
    lumiArray = array('d',[float(lumi) for lumi in options.lumi_array.split(',')])
    cfg = options.config
    
    model = options.model
    if options.mGluino>-1:
        massPoint = 'mGl-%i_mLSP-%i'%(options.mGluino,options.mLSP)
    elif options.mStop>-1:
        massPoint = 'mStop-%i_mLSP-%i'%(options.mStop,options.mLSP)

    btag = '0-3btag'

    fit = ''
    if options.fit:
        fit = '--fit'
        
    for lumi in lumiArray:
        for box in boxes:
            os.system('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w -q Backgrounds/*.root'%(cfg,box))
            os.system('python python/DustinTuple2RooDataSet.py -c %s -b %s -d Datasets/ -w Signals/RazorInclusive_SMS-%s_2J_%s_*.root'%(cfg,box,model,massPoint))
            os.system('python python/RooDataSet2UnweightedDataset.py -c %s -b %s -d Datasets/ Datasets/RazorInclusive_SMCocktail_weighted_lumi-%.1f_%s_%s.root'%(cfg,box,lumi,btag,box))
            signalDsName = 'Datasets/RazorInclusive_SMS-%s_2J_%s_weighted_lumi-%.1f_%s_%s.root'%(model,massPoint,lumi,btag,box)
            backgroundDsName = 'Datasets/RazorInclusive_SMCocktail_unweighted_lumi-%.1f_%s_%s.root'%(lumi,btag,box)
            os.system('python python/WriteDataCard.py -l %f -c %s -b %s -d %s %s %s %s'%(1000*lumi,cfg,box,options.outDir,fit,signalDsName,backgroundDsName))
            
            if signif:
                os.system('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 --toysFreq %s/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box))
                os.system('mv higgsCombine%s_%s_lumi-%.1f_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir))
            else:
                os.system('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s_%s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,btag,box))
                os.system('mv higgsCombine%s_%s_lumi-%.1f_%s_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,btag,box,options.outDir))
        if len(boxes)>1:
            for box in boxes: os.system('cp %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt .'%(options.outDir,model,massPoint,lumi,btag,box))
            cmds = ['%s=razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(box,model,massPoint,lumi,btag,box) for box in boxes]
            os.system('combineCards.py %s > %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(' '.join(cmds),options.outDir,model,massPoint,lumi,btag,options.box))
            os.system('cat %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(options.outDir,model,massPoint,lumi,btag_options.box))
            if signif:
                os.system('combine -M ProfileLikelihood --signif --expectSignal=1 -t -1 --toysFreq %s/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s_%s'%(options.outDir,model,massPoint,lumi,options.box,model,massPoint,lumi,btag,options.box))
                os.system('mv higgsCombine%s_%s_lumi-%.1f_%s_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,lumi,btag,options.box,options.outDir))
            else:
                os.system('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.1f_%s_%s.txt -n %s_%s_lumi-%.1f_%s'%(options.outDir,model,massPoint,lumi,options.box,model,massPoint,lumi,btag,options.box))
                os.system('mv higgsCombine%s_%s_lumi-%.1f_%s_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,btag,options.box,options.outDir))
            for box in boxes: os.system('rm razor_combine_%s_%s_lumi-%.1f_%s_%s.txt'%(model,massPoint,lumi,btag,box))
 
