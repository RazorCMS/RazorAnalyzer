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
    parser.add_argument("--closure", action="store_true", help="include uncertainties from scale factor cross check")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    closure = args.closure

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    #Process inclusive sample twice; the first pass will compute the overall normalization 
    #and the second pass will be a rerun with the corrected normalization
    if closure:
        regionsOrder = ["DYJetsDileptonInv", "DYJetsDileptonInvDiJet", "DYJetsDileptonInvMultiJet"]
    else:
        regionsOrder = ["DYJetsDileptonInvUncorr", "DYJetsDileptonInv", 
                "DYJetsDileptonInvDiJet", "DYJetsDileptonInvMultiJet",
                "DYJetsDileptonInvDiJetWJetsCorr", "DYJetsDileptonInvMultiJetWJetsCorr",
                "DYJetsDileptonInvNoSFs"]
    regions = {
            "DYJetsDileptonInvUncorr":Analysis("DYJetsDileptonInv",tag=tag),
            "DYJetsDileptonInv":Analysis("DYJetsDileptonInv",tag=tag),
            "DYJetsDileptonInvDiJet":Analysis("DYJetsDileptonInv",tag=tag,njetsMin=2,njetsMax=3),
            "DYJetsDileptonInvMultiJet":Analysis("DYJetsDileptonInvMultiJet",tag=tag,njetsMin=4),
            "DYJetsDileptonInvDiJetWJetsCorr":Analysis("DYJetsDileptonInv",tag=tag,njetsMin=2,njetsMax=3),
            "DYJetsDileptonInvMultiJetWJetsCorr":Analysis("DYJetsDileptonInvMultiJet",tag=tag,njetsMin=4),
            "DYJetsDileptonInvNoSFs":Analysis("DYJetsDileptonInv",tag=tag),
            }
    sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
    #make two dictionaries of scale factor histograms, one with GJets and one with WJets corrections
    sfHists = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["DYJetsDileptonInvDiJet"].samples, 
            scaleFactorNames={ "DYJetsInv":"GJetsInv" }, debugLevel=debugLevel )
    sfHistsForWCorr = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["DYJetsDileptonInvDiJet"].samples, 
            scaleFactorNames={ "DYJetsInv":"WJetsInv" }, debugLevel=debugLevel )
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    for d in [sfHists, sfHistsForWCorr]:
        d['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
        d['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
        d['NJetsInv'] = sfNJetsFile.Get("GJetsInvScaleFactors")
        d['NJetsWJetsInv'] = sfNJetsFile.Get("WJetsInvScaleFactors")
    sfVars = { "WJets":("MR","Rsq"), "TTJets":("MR","Rsq"), "DYJetsInv":("MR_NoZ","Rsq_NoZ") }
    outfile = rt.TFile(
        "data/ScaleFactors/RazorMADD2015/RazorDYJetsDileptonInvCrossCheck_%s.root"%(tag), "RECREATE")
    #optionally inflate scale factor uncertainties to cover difference between G+jets and W+jets SFs
    if args.closure:
        sfFile = rt.TFile.Open(sfFilename)
        downHist = sfFile.Get("GJetsInvScaleFactors_Down")
        for bn in range(sfHists["DYJetsInv"].GetNumberOfBins()+1):
            err = sfHists["DYJetsInv"].GetBinError(bn)
            sysErr = sfHists["DYJetsInv"].GetBinContent(bn) - downHist.GetBinContent(bn)
            newErr = ( err*err + sysErr*sysErr )**(0.5)
            sfHists["DYJetsInv"].SetBinError(bn-1, newErr) # adjust bin number by 1 to account for bug in ROOT
            print "Increasing error on DYJets scale factor bin",bn,"from",err,"to",sfHists["DYJetsInv"].GetBinError(bn)
    
    updateNorm = (not closure)
    for region in regionsOrder:
        analysis = regions[region]
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        if closure:
            outdir += '_Closure'
        os.system('mkdir -p '+outdir)
        #prepare analysis
        auxSFs = razorWeights.getNJetsSFs(analysis,jetName='NJets_NoZ')
        #use the correct set of scale factors
        if 'WJetsCorr' in region:
            sfHistsToUse = sfHistsForWCorr
            auxSFs['DYJetsInv'] = {'NJetsWJetsInv': ('NJets_NoZ', '1')}
        elif 'NoSFs' in region:
            sfHistsToUse = {}
            auxSFs = { proc:{} for proc in auxSFs }
        else:
            sfHistsToUse = sfHists
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHistsToUse,
            sfVars = sfVars, printdir=outdir, auxSFs=auxSFs, debugLevel=debugLevel )
        #record discrepancies > 1 sigma
        tmpSFHists = copy.copy(sfHists)
        del tmpSFHists["DYJetsInv"]
        appendScaleFactors("DYJetsInv", hists, tmpSFHists, lumiData=analysis.lumi, 
            debugLevel=debugLevel, var=sfVars["DYJetsInv"], signifThreshold=1.0, printdir=outdir)
        #write out scale factors
        print "Writing histogram",tmpSFHists["DYJetsInv"].GetName(),"to file"
        outfile.cd()
        tmpSFHists["DYJetsInv"].Write(region+"ScaleFactors")

        #in the first pass, update the normalization of the G+jets scale factor histogram
        if updateNorm:
            dataNorm = hists["Data"]["1"].GetBinContent(1)
            dataNormErr = hists["Data"]["1"].GetBinError(1)
            mcNorm = 0
            mcNormErr = 0
            for proc in analysis.samples:
                mcNorm += hists[proc]["1"].GetBinContent(1)
                mcNormErr += (hists[proc]["1"].GetBinError(1))**2
            mcNormErr = mcNormErr**(0.5)
            normUpdate = dataNorm / mcNorm
            normUpdateErr = ( (dataNormErr/mcNorm)**2 + (dataNorm*mcNormErr/(mcNorm*mcNorm))**2 )**(0.5)

            sfFilenameUncorr="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s_Uncorr.root"%(tag)
            print "Writing old G+jets scale factor histogram to",sfFilenameUncorr
            sfFileUncorr = rt.TFile.Open(sfFilenameUncorr, "RECREATE")
            sfHistUncorr = sfHists["DYJetsInv"].Clone()
            sfFileUncorr.WriteTObject(sfHistUncorr, "GJetsInvScaleFactors", "WriteDelete")
            sfFileUncorr.Close()
            print "Scaling G+jets scale factor histogram by",normUpdate,"( = %.3f / %.3f )"%(dataNorm,mcNorm)
            print "Propagating error on scale factor:",normUpdateErr
            for bn in range(sfHists["DYJetsInv"].GetNumberOfBins()+1):
                sfHists["DYJetsInv"].SetBinContent( bn, sfHists["DYJetsInv"].GetBinContent(bn)*normUpdate )
                sfHists["DYJetsInv"].SetBinError( bn-1, sfHists["DYJetsInv"].GetBinError(bn)*normUpdate )
                sfHists["DYJetsInv"].SetBinError( bn-1, 
                        ( ( sfHists["DYJetsInv"].GetBinError(bn))**2 + normUpdateErr*normUpdateErr )**(0.5) )
            sfHistCorr = sfHists["DYJetsInv"].Clone()
            #Write the corrected G+jets scale factor file over the old one
            sfFile = rt.TFile.Open(sfFilename, "UPDATE")
            sfFile.WriteTObject(sfHistCorr, "GJetsInvScaleFactors", "WriteDelete")
            sfFile.Close()
            updateNorm = False

        #export histograms
        macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                outDir=outdir, debugLevel=debugLevel )

    outfile.Close()

