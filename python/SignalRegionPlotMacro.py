import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
from macro.razorAnalysis import Analysis
from macro.razorMacros import plotControlSampleHists
import macro.macro as macro
from SignalRegionMacro import shapes

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
    parser.add_argument('--box', help="choose a box",required=True)
    parser.add_argument('--btags', help="choose number of b-tags", type=int, required=True)
    parser.add_argument('--no-sfs', help='Uncertainties from scale factor cross checks will be ignored', action='store_true',dest='noSFs')
    parser.add_argument('--no-sys', help='Do not propagate systematic uncertainties', action='store_true',dest='noSys')
    parser.add_argument("--tag", dest="tag", required=True,
            help="Analysis tag, e.g. Razor2015")
    parser.add_argument("--dir", help="Specify input/output directory")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    box = args.box
    btags = args.btags
    dirName = 'Plots/%s/%s%dB'%(tag,box,btags)

    plotOpts = {"SUS15004":True}
    plotOpts["combineBackgrounds"] = { 
            "Other":["SingleTop","DYJets","Other"], 
            "TTJets":["TTJets1L","TTJets2L","TTJets"] }

    analysis = Analysis(box, tag=tag, nbMin=btags, nbMax=btags)

    #apply options
    blindBins = [(x,y) for x in range(2,len(analysis.binning["MR"])+1) 
            for y in range(2,len(analysis.binning["Rsq"])+1)]
    if args.unblind: blindBins = None
    samples = analysis.samples
    plotOpts["combineSamples"] = analysis.samplesReduced

    shapesToUse = copy.copy(shapes[box])

    #option to disable scale factors
    if args.noSFs:
        print "Ignoring all uncertainties from scale factor cross checks."
        toRemove = ['btaginvcrosscheck','btagcrosscheckrsq','btagcrosscheckmr','sfsyszinv','ttcrosscheck','zllcrosscheck','sfsysttjets','sfsyswjets','vetolepptcrosscheck','vetotauptcrosscheck','vetolepetacrosscheck','vetotauetacrosscheck']
        #remove scale factor cross check uncertainties
        shapesToUse = [s for s in shapesToUse if s not in toRemove]
        #this removes scale factor uncertainties that are listed as tuples
        shapesToUse = [s for s in shapesToUse if not (hasattr(s, '__getitem__') and s[0] in toRemove)] 
        dirName += 'NoSFs'
    ###### TEMPORARY UNTIL MISTAG WEIGHTS ARE FIXED #######
    if 'bmistag' in shapesToUse:
        shapesToUse.remove('bmistag')
    #######################################################
    #option to disable systematics
    if args.noSys:
        print "Ignoring systematic uncertainties."
        shapesToUse = []
        dirName += 'NoSys'

    print "\n---",box,"Box,",btags,"B-tags ---"

    if args.dir is not None:
        dirName = args.dir

    extbox = box+str(btags)+"B"
    unrollBins = analysis.unrollBins
    inFile = dirName+'/razorHistograms'+extbox+'.root'
    if args.noSFs:
        inFile = inFile.replace('.root','NoSFs.root')
    if args.noSys:
        inFile = inFile.replace('.root','NoSys.root')
    lumi = analysis.lumi

    plotControlSampleHists(box, inFile, samples=samples, plotOpts=plotOpts, boxName=box, 
            btags=btags, blindBins=blindBins, debugLevel=debugLevel, printdir=dirName, lumiData=lumi,
            unrollBins=unrollBins, shapeErrors=shapesToUse)
