import sys, os, argparse, copy
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
    parser.add_argument("--tag", dest="tag", default="Razor2015",
                                help="Analysis tag, e.g. Razor2015")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    if tag not in ["Razor2015","Razor2016"]:
        sys.exit("Error: tag "+tag+" not supported!")

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regions = {
            "DYJetsDileptonInvDiJet":Analysis("DYJetsDileptonInv",tag=tag,njetsMin=2,njetsMax=3),
            "DYJetsDileptonInvMultiJet":Analysis("DYJetsDileptonInvMultiJet",tag=tag,njetsMin=4),
            }
    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions["DYJetsDileptonInvDiJet"].samples, 
            scaleFactorNames={ "DYJetsInv":"GJetsInv" }, debugLevel=debugLevel)
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJets'] = sfNJetsFile.Get("NJetsCorrectionScaleFactors")
    sfHists['NJetsInv'] = sfNJetsFile.Get("GJetsScaleFactorVsNJets")
    sfVars = { "WJets":("MR","Rsq"), "TTJets":("MR","Rsq"), "DYJetsInv":("MR_NoZ","Rsq_NoZ") }
    outfile = rt.TFile(
        "data/ScaleFactors/RazorMADD2015/RazorDYJetsDileptonInvCrossCheck_%s.root"%(tag), "RECREATE")
    
    for region,analysis in regions.iteritems():
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #prepare analysis
        auxSFs = razorWeights.getNJetsSFs(analysis,jetName='NJets_NoZ')
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
            sfVars = sfVars, printdir=outdir, auxSFs=auxSFs, debugLevel=debugLevel )
        #record discrepancies > 1 sigma
        tmpSFHists = copy.copy(sfHists)
        del tmpSFHists["DYJetsInv"]
        appendScaleFactors("DYJetsInv", hists, tmpSFHists, lumiData=analysis.lumi, 
            debugLevel=debugLevel, var=sfVars["DYJetsInv"], signifThreshold=1.0, printdir=outdir)
        #export histograms
        macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                outDir=outdir, debugLevel=debugLevel )
        #write out scale factors
        print "Writing histogram",tmpSFHists["DYJetsInv"].GetName(),"to file"
        outfile.cd()
        tmpSFHists["DYJetsInv"].Write(region+"ScaleFactors")

    outfile.Close()

