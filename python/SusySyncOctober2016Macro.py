import sys, os, argparse
import ROOT as rt

from macro import macro
from macro.razorAnalysis import Analysis
from macro.razorMacros import makeControlSampleHistsForAnalysis

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = "Razor2016"

    #initialize
    plotOpts = { 'comment':False }

    regionsOrder = ["0BSusySync", "1BSusySync", "2BSusySync"]
    regions = {
            "0BSusySync":Analysis("SusySync", tag=tag, nbMax=0),
            "1BSusySync":Analysis("SusySync", tag=tag, nbMin=1),
            "2BSusySync":Analysis("SusySync", tag=tag, nbMin=2),
            }

    for region in regionsOrder:
        analysis = regions[region]
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                printdir=outdir, debugLevel=debugLevel ) 
        #export histograms
        macro.exportHists(hists, outFileName='controlHistograms'+region+'.root', 
            outDir=outdir, debugLevel=debugLevel)
