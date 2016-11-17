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
    parser.add_argument("--tag", dest="tag", default="Razor2016",
                                help="Analysis tag, e.g. Razor2015")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag

    #initialize
    plotOpts = { 'comment':False, "SUS15004CR":True }
    regions = {}
    #define all tests
    for name,jets in {"DiJet":(2,3),"MultiJet":(4,-1)}.iteritems():
        regionName = "GJetsInv"+name+"ClosureTest"
        regions[regionName] = Analysis("GJetsInv",tag=tag,
                njetsMin=jets[0], njetsMax=jets[1])

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions.itervalues().next().samples, debugLevel=debugLevel)
    sfVars = ("MR_NoPho","Rsq_NoPho")
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsInv'] = sfNJetsFile.Get("GJetsInvScaleFactors")
    razorWeights.loadPhotonPurityHists(sfHists, tag, debugLevel)
    for region,analysis in regions.iteritems():
        print "\nRegion:",region,"\n"
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #set up analysis
        (xbins,cols) = analysis.unrollBins
        auxSFs = razorWeights.getNJetsSFs(analysis, jetName=analysis.jetVar)
        razorWeights.getPhotonPuritySFs(auxSFs)
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, btags=analysis.nbMin,
                dataDrivenQCD=True, debugLevel=debugLevel )
        #compute scale factors
        appendScaleFactors( region+"NBJets", hists, sfHists, lumiData=analysis.lumi, 
                var="NBJetsMedium", debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir )
        #export histograms
        macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                outDir=outdir, debugLevel=debugLevel )
        #write out scale factors
        outfile = rt.TFile(
                "data/ScaleFactors/RazorMADD2015/RazorGJetsBTagClosureTests_%s.root"%(tag),
                "UPDATE")
        print "Writing scale factor histogram",sfHists[region+"NBJets"].GetName(),"to file"
        outfile.cd()
        sfHists[region+"NBJets"].Write( sfHists[region+"NBJets"].GetName() )
        outfile.Close()
