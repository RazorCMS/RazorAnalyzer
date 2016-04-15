import sys,os
import argparse
import ROOT as rt
import glob

#local imports
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import makeRazorBinEvidencePlots
from RunCombine import exec_me
import macro.macro as macro
from SidebandMacro import SAMPLES, LUMI, config
from framework import Config

SIGNAL_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined"
NOPATHOLOGIES_SIGNAL_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_RemovedPathologicalJets20160414/combined"
BACKGROUND_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorMADD2015"

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--dir', help="output directory", default="SignalRegionPlots", dest='outDir')
    parser.add_argument('-m','--model', default="T1bbbb", help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=int, help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=int, help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=int, help="mass of LSP")
    parser.add_argument('--no-pathologies', dest="noPathologies", action='store_true', 
            help='Remove pathological fastsim events')

    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    outDir = args.outDir

    if args.box is None:
        print "Please choose an analysis box with --box"
        sys.exit()
    box = args.box
    boxList = box.split('_') #interpret box1_box2_... as a list of boxes to combine

    dirToUse = SIGNAL_DIR
    if args.noPathologies:
        dirToUse = NOPATHOLOGIES_SIGNAL_DIR

    #make output directory
    os.system('mkdir -p '+outDir)

    for curBox in boxList:
        #get configuration for this box
        samples = SAMPLES[curBox]
        unrollBins = [(xbinsSignal[curBox][str(btags)+'B'], colsSignal[curBox][str(btags)+'B']) for btags in range(4)]

        #call SMS template maker
        if 'T1' in args.model or 'T5' in args.model:
            modelName = 'SMS-'+args.model+'_'+str(args.mGluino)+'_'+str(args.mLSP)
        else:
            modelName = 'SMS-'+args.model+'_'+str(args.mStop)+'_'+str(args.mLSP)
        signalFilename=dirToUse+'/'+modelName+'.root'
        brString = ""
        if 'T1x' in args.model:
            xBR = float(args.model[args.model.find('x')+1:args.model.find('y')].replace('p','.'))
            yBR = float(args.model[args.model.find('y')+1:].replace('p','.'))
            brString = '--xBR %.2f --yBR %.2f'%(xBR,yBR)
            signalFilename = dirToUse+'/SMS-T1ttbb_'+str(args.mGluino)+'_'+str(args.mLSP)+'.root'
        signalHistFilename = '%s/%s_lumi-%.3f_0-3btag_%s.root'%(outDir,modelName,LUMI*1.0/1000,curBox)
        if not glob.glob(signalHistFilename):
            exec_me('python python/SMSTemplates.py --merge-bins -c %s -d %s --lumi %d --box %s %s %s --no-signal-sys'%(config, outDir, LUMI, curBox, signalFilename, brString), False) 
        #load SMS template histograms
        signalHists = macro.importHists(signalHistFilename)
        #update with correct names
        for x,h in signalHists.items():
            h.SetName(h.GetName().replace(curBox+'_'+args.model,modelName))
            signalHists[h.GetName()] = signalHists.pop(x)
        signalHist = signalHists[modelName]

        #make combined unrolled histograms for background
        makeRazorBinEvidencePlots(curBox, samples=samples, inDir=BACKGROUND_DIR, outDir=outDir, 
                signalHist=signalHist, unrollBins=unrollBins, debugLevel=debugLevel)

