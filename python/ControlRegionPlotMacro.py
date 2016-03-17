import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import plotControlSampleHists
import macro.macro as macro
from SidebandMacro import LUMI, MCLUMI
from ComputeScaleFactorsMacro import SAMPLES_TTJ1L, SAMPLES_WJ1L, SAMPLES_WJ1L_INV, ControlRegionBinning

SAMPLES = { "WJetControlRegion":SAMPLES_WJ1L, "WJetInvControlRegion":SAMPLES_WJ1L_INV, "TTJetsSingleLeptonControlRegion":SAMPLES_TTJ1L }

boxMapping = { "WJetControlRegion":"WJetsSingleLepton", "WJetInvControlRegion":"WJetsSingleLeptonInv", "TTJetsSingleLeptonControlRegion":"TTJetsSingleLepton" }

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--dir', help="output directory (should contain the ROOT files with the razor histograms)", default="ControlSamplePlots", dest='dirName')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    dirName = args.dirName

    plotOpts = { "ymin":1e-3, 'comment':False }

    boxesToUse = ["WJetControlRegion","TTJetsSingleLeptonControlRegion","WJetInvControlRegion"]
    if args.box is not None:
        boxesToUse = [args.box]

    #make output directory
    os.system('mkdir -p '+dirName)

    #estimate yields in signal region
    for boxName in boxesToUse:

        #apply options
        samplesToUse = SAMPLES[boxName]

        print "\n---",boxName,"---"

        unrollBins = (xbinsSignal[boxName]['0B'], colsSignal[boxName]['0B'])
        inFile = dirName+'/controlHistograms'+boxMapping[boxName]+'.root'

        plotControlSampleHists(boxName, inFile, samples=samplesToUse, plotOpts=plotOpts, lumiMC=MCLUMI, lumiData=LUMI, boxName=boxName, debugLevel=debugLevel, printdir=dirName, unrollBins=unrollBins)
