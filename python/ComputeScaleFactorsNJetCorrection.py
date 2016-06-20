import sys, os, argparse
import ROOT as rt

from macro import macro
from macro.razorAnalysis import Analysis
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--tag", dest="tag", default="Razor2015",
                                help="Analysis tag, e.g. Razor2015")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    if tag not in ["Razor2015","Razor2016"]:
        sys.exit("Error: tag "+tag+" not supported!")

    #initialize
    plotOpts = { "comment":False }
    regions = {
            "OneLeptonForNJets":Analysis("SingleLepton",tag=tag),
            "OneLeptonInvForNJets":Analysis("SingleLeptonInv",tag=tag)
            }
    sfVars = {
            "OneLeptonForNJets":("MR","Rsq"),
            "OneLeptonInvForNJets":{ "WJetsInv":("MR_NoW","Rsq_NoW"), "TTJets":("MR","Rsq") }
            }
    sfHists = { region:macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions[region].samples, debugLevel=debugLevel) for region in regions }
    sfNames = { "OneLeptonForNJets":"NJetsCorrection", "OneLeptonInvForNJets":"NJetsNoWCorrection" }
    njetsNames = { "OneLeptonForNJets":"NJets40", "OneLeptonInvForNJets":"NJets_NoW" }
    outfile = rt.TFile("data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag), 
            "RECREATE")

    for region,analysis in regions.iteritems():
        #make output directory
        outdir = 'Plots/'+region
        os.system('mkdir -p '+outdir)
        #set up analysis
        process = sfNames[region]
        (xbins,cols) = analysis.unrollBins
        sfHistsToUse = sfHists[region]
        sfVarsToUse = sfVars[region]
        njetsName = njetsNames[region]
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHistsToUse,
                sfVars=sfVarsToUse, printdir=outdir, debugLevel=debugLevel )
        #compute scale factors
        appendScaleFactors( process, hists, sfHistsToUse, lumiData=analysis.lumi, 
                debugLevel=debugLevel, var=njetsName, printdir=outdir )
        #export histograms
        macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                outDir=outdir, debugLevel=debugLevel )
        #write out scale factors
        print "Writing scale factor histogram",sfHistsToUse[process].GetName(),"to file"
        outfile.cd()
        sfHistsToUse[process].Write( sfHistsToUse[process].GetName() )

    outfile.Close()
