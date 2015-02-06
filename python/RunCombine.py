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
    parser.add_option('--lumi-array',dest="lumi_array", default="0.2,4,10",type="string",
                  help="lumi array in fb^-1, e.g.: 0.2,4,10")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (use MC directly)")


    (options,args) = parser.parse_args()

    boxes = options.box.split('_')
    
    lumiArray = array('d',[float(lumi) for lumi in options.lumi_array.split(',')])
    cfg = options.config
    
    model = options.model
    if options.mGluino>-1:
        massPoint = 'mGl-%i_mLSP-%i'%(options.mGluino,options.mLSP)
    elif options.mStop>-1:
        massPoint = 'mStop-%i_mLSP-%i'%(options.mStop,options.mLSP)
        
    lumi_in = 4000.


    noFit = ''
    if options.noFit:
        noFit = '--no-fit'
        
    for lumi in lumiArray:
        for box in boxes:
            signalDsName = 'Datasets/RazorAnalysis_SMS-%s_2J_%s_25ns_weighted_lumi-%.1f_%s.root'%(model,massPoint,lumi_in/1000,box)
            backgroundDsName = 'Datasets/RazorAnalysis_SMCocktail_weighted_lumi-%.1f_%s.root'%(lumi_in/1000,box)
            os.system('python python/WriteDataCard.py -l %f -c %s -b %s -d %s %s %s %s'%(1000*lumi,cfg,box,options.outDir,noFit,signalDsName,backgroundDsName))
            os.system('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s'%(options.outDir,model,massPoint,lumi,box,model,massPoint,lumi,box))
            os.system('mv higgsCombine%s_%s_lumi-%.1f_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,box,options.outDir))
        if len(boxes)>1:
            for box in boxes: os.system('cp %s/razor_combine_%s_%s_lumi-%.1f_%s.txt .'%(options.outDir,model,massPoint,lumi,box))
            cmds = ['%s=razor_combine_%s_%s_lumi-%.1f_%s.txt'%(box,model,massPoint,lumi,box) for box in boxes]
            os.system('combineCards.py %s > %s/razor_combine_%s_%s_lumi-%.1f_%s.txt'%(' '.join(cmds),options.outDir,model,massPoint,lumi,options.box))
            os.system('cat %s/razor_combine_%s_%s_lumi-%.1f_%s.txt'%(options.outDir,model,massPoint,lumi,options.box))
            os.system('combine -M Asymptotic %s/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s'%(options.outDir,model,massPoint,lumi,options.box,model,massPoint,lumi,options.box))
            os.system('mv higgsCombine%s_%s_lumi-%.1f_%s.Asymptotic.mH120.root %s/'%(model,massPoint,lumi,options.box,options.outDir))
            for box in boxes: os.system('rm razor_combine_%s_%s_lumi-%.1f_%s.txt'%(model,massPoint,lumi,box))
 
