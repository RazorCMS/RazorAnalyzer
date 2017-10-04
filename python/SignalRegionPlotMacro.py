import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
from macro.razorAnalysis import Analysis
from macro.razorMacros import plotControlSampleHists
import macro.macro as macro
import SignalRegionMacro as sig

def adjustFilenameForOptions(fName, args):
    if args.noSFs:
        fName = fName.replace('.root','NoSFs.root')
    if args.noSys:
        fName = fName.replace('.root','NoSys.root')
    if args.nloZInv:
        fName = fName.replace('.root','NLOZInv.root')
    if args.fineGrained:
        fName = fName.replace('.root', 'FineGrained.root')
    return fName

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = sig.makeSignalRegionParser()
    parser.add_argument("--dir", help="Specify input/output directory")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    box = args.box
    btags = args.btags
    dirName = 'Plots/%s/%s%dB'%(tag,box,btags)+sig.getDirSuffix(args)

    analysis = Analysis(box, tag=tag, nbMin=btags, nbMax=btags)
    plotOpts = sig.getPlotOpts(args, analysis)

    blindBins = [(x,y) for x in range(1,len(analysis.binning["MR"])+1) 
            for y in range(1,len(analysis.binning["Rsq"])+1)]
    if args.unblind: blindBins = None

    shapesToUse = copy.copy(sig.shapes[box])

    if args.noSFs:
        print "Ignoring all uncertainties from scale factor cross checks."
        sig.removeSFShapes(shapesToUse)
    if args.noSys:
        print "Ignoring systematic uncertainties."
        shapesToUse = []
    if args.fineGrained:
        analysis, _, _, shapesToUse, plotOpts = sig.adjustForFineGrainedMCPred(
                analysis, {}, {}, shapesToUse, plotOpts)
    
    if args.dir is not None:
        dirName = args.dir

    print "\n---",box,"Box,",btags,"B-tags ---"

    extbox = box+str(btags)+"B"
    unrollBins = analysis.unrollBins
    lumi = analysis.lumi
    samples = analysis.samples
    inFile = dirName+'/razorHistograms'+extbox+'.root'
    inFile = adjustFilenameForOptions(inFile, args)
    print "Input file: {}".format(inFile)

    plotControlSampleHists(box, inFile, samples=samples, plotOpts=plotOpts, boxName=box, 
            btags=btags, blindBins=blindBins, debugLevel=debugLevel, printdir=dirName, lumiData=lumi,
            unrollBins=unrollBins, shapeErrors=shapesToUse, 
            doEmptyBinErrs=args.fineGrained)
