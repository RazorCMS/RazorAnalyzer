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
    parser.add_argument("--tag", dest="tag", default="Razor2016",
                                help="Analysis tag, e.g. Razor2015")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regionsOrder = ["TTJetsDilepton","TTJetsDileptonMultiJet","TTJetsDileptonDiJet"]
    regions = {
            "TTJetsDilepton":Analysis("TTJetsDilepton",tag=tag),
            "TTJetsDileptonDiJet":Analysis("TTJetsDilepton",tag=tag,njetsMin=2,njetsMax=3),
            "TTJetsDileptonMultiJet":Analysis("TTJetsDileptonMultiJet",tag=tag,njetsMin=4),
            }
    sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
    sfHists = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["TTJetsDilepton"].samples, debugLevel=debugLevel )
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
    outfile = rt.TFile(
        "data/ScaleFactors/RazorMADD2015/RazorTTJetsDileptonCrossCheck_%s.root"%(tag), "RECREATE")
    
    for region in regionsOrder:
        analysis = regions[region]
        analysis.weightOpts.append('ttbardileptonmt')
        analysis.dataWeightOpts.append('ttbardileptonmt')
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #prepare analysis
        auxSFs = razorWeights.getNJetsSFs(analysis,jetName='NJets40')
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
            printdir=outdir, auxSFs=auxSFs, debugLevel=debugLevel )
        #record discrepancies > 1 sigma
        tmpSFHists = copy.copy(sfHists)
        if 'TTJets2L' in tmpSFHists: del tmpSFHists["TTJets2L"]
        appendScaleFactors("TTJets2L", hists, tmpSFHists, lumiData=analysis.lumi, 
            debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir)
        #write out scale factors
        print "Writing histogram",tmpSFHists["TTJets2L"].GetName(),"to file"
        outfile.cd()
        tmpSFHists["TTJets2L"].Write(region+"ScaleFactors")
        #export histograms
        macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                outDir=outdir, debugLevel=debugLevel )

    outfile.Close()

