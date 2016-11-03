### This script is like ComputeScaleFactorsNJetsCorrection.py except that it 
### computes the NJets correction separately for TT+jets and W+jets
import sys, os, argparse
import ROOT as rt

from macro import macro, razorWeights
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
    parser.add_argument("--tag", help="Analysis tag, e.g. Razor2015", default="Razor2016")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag

    # 1) compute the TTJets scale factor for the 2-3 and >=4 jet bins
    # 2) compute the WJets scale factors for all bins, applying the TTJets ones
    # 3) compute the TTJets scale factor for the 1-jet bin, applying the WJets one
    plotOpts = { 'comment':False, "SUS15004CR":True }
    regionsOrder = ["TTJetsForNJets", "WJetsForNJets", "TTJetsForNJetsCorrected"]
    regions = {
            "TTJetsForNJets":Analysis("TTJetsSingleLepton",tag=tag,nbMin=1),
            "TTJetsForNJetsCorrected":Analysis("TTJetsSingleLepton",tag=tag,nbMin=1),
            "WJetsForNJets":Analysis("WJetsSingleLepton",tag=tag,nbMax=0),
            }
    sfVars = ("MR","Rsq")
    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions["TTJetsForNJets"].samples, 
            debugLevel=debugLevel)
    sfNames = { 
            "TTJetsForNJets":"TTJets", 
            "TTJetsForNJetsCorrected":"TTJets", 
            "WJetsForNJets":"WJets", 
            }
    njetsName = "NJets40"
    outfile = rt.TFile(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag), 
            "UPDATE")

    for region in regionsOrder:
        analysis = regions[region]
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #set up analysis
        process = sfNames[region]
        (xbins,cols) = analysis.unrollBins
        #use proper set of NJets correction factors.
        #do not correct WJets until WJets scale factors have been derived
        auxSFs = razorWeights.getNJetsSFs(analysis)
        if region == "TTJetsForNJets" or region == "WJetsForNJets":
            auxSFs["WJets"] = {}
        if region == "TTJetsForNJets" or region == "TTJetsForNJetsCorrected":
            auxSFs["TTJets"] = {}
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                sfHists=sfHists, sfVars=sfVars, printdir=outdir, 
                auxSFs=auxSFs, debugLevel=debugLevel )
        #compute scale factors
        sfHistsCopy = sfHists.copy()
        appendScaleFactors( process, hists, sfHistsCopy, lumiData=analysis.lumi, 
                debugLevel=debugLevel, var=njetsName, printdir=outdir )
        if region == "TTJetsForNJets":
            sfHists["NJetsTTJets"] = sfHistsCopy["TTJets"]
            print "Setting 1-jet scale factor for TTJets to 1.0"
            sfHists["NJetsTTJets"].SetBinContent(1, 1.0)
            #store the tt+jets histogram for later
            tmpTTJets = sfHists["NJetsTTJets"].Clone()
        elif region == "WJetsForNJets":
            sfHists["NJetsWJets"] = sfHistsCopy["WJets"]
        elif region == "TTJetsForNJetsCorrected":
            print "Recovering 2-3 and 4+ jet scale factors from earlier..."
            for i in [2,3]:
                sfHistsCopy["TTJets"].SetBinContent(i, tmpTTJets.GetBinContent(i))
                sfHistsCopy["TTJets"].SetBinError(i, tmpTTJets.GetBinError(i))
        #export histograms
        macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                outDir=outdir, debugLevel=debugLevel )
        #write out scale factors
        if region != "TTJetsForNJets":
            print "Writing scale factor histogram",sfHistsCopy[process].GetName(),"to file"
            outfile.cd()
            sfHistsCopy[process].Write( sfHistsCopy[process].GetName() )

    outfile.Close()
