from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('--lumi-array',dest="lumi_array", default="0.2,4,10",type="string",
                  help="lumi array in fb^-1, e.g.: 0.2,4,10")


    (options,args) = parser.parse_args()

    lumiArray = array('d',[float(lumi) for lumi in options.lumi_array.split(',')])
    box = options.box
    cfg = options.config
    

    backgroundFileNames = []
    signalFileName = ''
    model = ''
    massPoint = ''
        
    for f in args:
        if f.lower().endswith('.root'):
            if f.lower().find('sms')!=-1:
                model = f.split('-')[1].split('_')[0]
                massPoint = '_'.join(f.split('_')[3:5])
                signalFileName = f
            else:
                backgroundFileNames.append(f)
    
    for lumi in lumiArray:
        os.system('python python/DustinTuple2RooDataSet.py -l %f -w -c %s -b %s -d Datasets/ %s'%(1000*lumi,cfg,box,' '.join(backgroundFileNames)))
        os.system('python python/DustinTuple2RooDataSet.py -l %f -w -c %s -b %s -d Datasets/ %s'%(1000*lumi,cfg,box,signalFileName))
        signalDsName = signalFileName.replace('.root','_lumi-%.1f_%s.root'%(lumi,box)).replace('Signals/','Datasets/')
        backgroundDsName = 'Datasets/RazorAnalysis_SMCocktail_weighted_lumi-%.1f_%s.root'%(lumi,box)
        os.system('python python/WriteDataCard.py -l %f -c %s -b %s -d cards/ %s %s'%(1000*lumi,cfg,box,signalDsName,backgroundDsName))
        os.system('combine -M Asymptotic cards/razor_combine_%s_%s_lumi-%.1f_%s.txt -n %s_%s_lumi-%.1f_%s'%(model,massPoint,lumi,box,model,massPoint,lumi,box))
        os.system('mv higgsCombine%s_%s_lumi-%.1f_%s.Asymptotic.mH120.root cards/'%(model,massPoint,lumi,box))
