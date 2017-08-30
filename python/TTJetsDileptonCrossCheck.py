import sys, os, copy
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regionsOrder = ["TTJetsDilepton","TTJetsDileptonMultiJet","TTJetsDileptonDiJet"]
    regions = {
            "TTJetsDilepton":Analysis("TTJetsDilepton",tag=tag,boostCuts=boostCuts),
            "TTJetsDileptonDiJet":Analysis("TTJetsDilepton",tag=tag,njetsMin=2,njetsMax=3,boostCuts=boostCuts),
            "TTJetsDileptonMultiJet":Analysis("TTJetsDileptonMultiJet",tag=tag,njetsMin=4,boostCuts=boostCuts),
            }
    sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
    sfHists = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["TTJetsDilepton"].samples, debugLevel=debugLevel )
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
    if not args.noSave:
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
            printdir=outdir, auxSFs=auxSFs, debugLevel=debugLevel, noFill=args.noFill )
        #record discrepancies > 1 sigma
        tmpSFHists = copy.copy(sfHists)
        if 'TTJets2L' in tmpSFHists: del tmpSFHists["TTJets2L"]
        appendScaleFactors("TTJets2L", hists, tmpSFHists, lumiData=analysis.lumi, 
            debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir)
        #write out scale factors
        print "Writing histogram",tmpSFHists["TTJets2L"].GetName(),"to file"
        if not args.noSave:
            outfile.cd()
            tmpSFHists["TTJets2L"].Write(region+"ScaleFactors")
            #export histograms
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )

    if not args.noSave:
        outfile.Close()

