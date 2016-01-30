import sys,os
import argparse
import ROOT as rt

#local imports
from framework import Config
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import plotControlSampleHists
import macro.macro as macro

LUMI = 2185 #in /pb
MCLUMI = 1 

SAMPLES_HADRONIC = ["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets2L", "TTJets1L"]
SAMPLES_LEPTONIC = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"]
SAMPLES = { "MultiJet":SAMPLES_HADRONIC, "MuMultiJet":SAMPLES_LEPTONIC, "EleMultiJet":SAMPLES_LEPTONIC }
BOXES = ["MultiJet", "MuMultiJet", "EleMultiJet"]

config = "config/run2_20151108_Preapproval_2b3b_data.config"
cfg = Config.Config(config)
binsMRHad = cfg.getBinning("MultiJet")[0]
binsRsqHad = cfg.getBinning("MultiJet")[1]
hadronicBinning = { "MR":binsMRHad, "Rsq":binsRsqHad, ("MR","Rsq"):[] }
binsMRLep = cfg.getBinning("MuMultiJet")[0]
binsRsqLep = cfg.getBinning("MuMultiJet")[1]
leptonicBinning = { "MR":binsMRLep, "Rsq":binsRsqLep, ("MR","Rsq"):[] }
binning = { "MultiJet":hadronicBinning, "MuMultiJet":leptonicBinning, "EleMultiJet":leptonicBinning}

blindBins = {b:[(x,y) for x in range(2,len(binning[b]["MR"])+1) for y in range(2,len(binning[b]["Rsq"])+1)] for b in binning}

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--unblind", help="do not blind signal sensitive region", action='store_true')
    parser.add_argument("--full", help="label plot as full fit", action='store_true')
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--btags', help="choose number of b-tags", type=int)
    parser.add_argument('--dir', help="output directory (should contain the ROOT files with the razor histograms)", default="SignalRegionPlots", dest='dirName')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    dirName = args.dirName

    plotOpts = {"ymin":1e-3}
    doSideband=(not args.full)
    if not doSideband:
        plotOpts['sideband'] = False
    else:
        plotOpts['sideband'] = True

    boxesToUse = BOXES
    if args.box is not None:
        boxesToUse = [args.box]

    #make output directory
    os.system('mkdir -p '+dirName)

    #estimate yields in signal region
    for boxName in boxesToUse:

        #apply options
        blindBinsToUse = blindBins[boxName]
        if args.unblind: blindBinsToUse = None
        samplesToUse = SAMPLES[boxName]

        #loop over btag bins
        if args.btags is not None:
            btaglist = [args.btags]
        else:
            btaglist = [0,1,2,3]
        for btags in btaglist:
            print "\n---",boxName,"Box,",btags,"B-tags ---"

            extboxName = boxName+str(btags)+"BTag"
            unrollBins = (xbinsSignal[boxName][str(btags)+'B'], colsSignal[boxName][str(btags)+'B'])
            inFile = dirName+'/razorHistograms'+extboxName+'.root'

            plotControlSampleHists(extboxName, inFile, samples=samplesToUse, plotOpts=plotOpts, lumiMC=MCLUMI, lumiData=LUMI, boxName=boxName, btags=btags, blindBins=blindBinsToUse, debugLevel=debugLevel, printdir=dirName, unrollBins=unrollBins)
