import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import plotControlSampleHists
import macro.macro as macro
from SidebandMacro import SAMPLES, LUMI, MCLUMI, blindBins, shapes

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
    parser.add_argument('--no-sfs', help='Uncertainties from scale factor cross checks will be ignored', action='store_true',dest='noSFs')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    dirName = args.dirName

    plotOpts = {"SUS15004":True}
    doSideband=(not args.full)
    if not doSideband:
        plotOpts['sideband'] = False
    else:
        plotOpts['sideband'] = True

    boxesToUse = ["MultiJet","MuMultiJet","EleMultiJet"]
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

        shapesToUse = copy.copy(shapes[boxName])
        #option to disable scale factors
        if args.noSFs:
            print "Ignoring all uncertainties from scale factor cross checks."
            toRemove = ['btaginvcrosscheck','btagcrosscheckrsq','btagcrosscheckmr','sfsyszinv','ttcrosscheck','zllcrosscheck','sfsysttjets','sfsyswjets','vetolepptcrosscheck','vetotauptcrosscheck','vetolepetacrosscheck','vetotauetacrosscheck']
            #remove scale factor cross check uncertainties
            shapesToUse = [s for s in shapesToUse if s not in toRemove]
            #this removes scale factor uncertainties that are listed as tuples
            shapesToUse = [s for s in shapesToUse if not (hasattr(s, '__getitem__') and s[0] in toRemove)] 

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

            plotControlSampleHists(extboxName, inFile, samples=samplesToUse, plotOpts=plotOpts, lumiMC=MCLUMI, lumiData=LUMI, boxName=boxName, btags=btags, blindBins=blindBinsToUse, debugLevel=debugLevel, printdir=dirName, unrollBins=unrollBins, shapeErrors=shapesToUse)
