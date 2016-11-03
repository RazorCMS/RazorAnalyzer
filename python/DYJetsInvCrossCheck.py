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
    #Process inclusive sample twice; the first pass will compute the overall normalization 
    #and the second pass will be a rerun with the corrected normalization
    regionsOrder = ["DYJetsDileptonInvUncorr", "DYJetsDileptonInv", 
            "DYJetsDileptonInvDiJet", "DYJetsDileptonInvMultiJet"]
    regions = {
            "DYJetsDileptonInvUncorr":Analysis("DYJetsDileptonInv",tag=tag),
            "DYJetsDileptonInv":Analysis("DYJetsDileptonInv",tag=tag),
            "DYJetsDileptonInvDiJet":Analysis("DYJetsDileptonInv",tag=tag,njetsMin=2,njetsMax=3),
            "DYJetsDileptonInvMultiJet":Analysis("DYJetsDileptonInvMultiJet",tag=tag,njetsMin=4),
            }
    sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
    sfHists = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["DYJetsDileptonInvDiJet"].samples, 
            scaleFactorNames={ "DYJetsInv":"GJetsInv" }, debugLevel=debugLevel )
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
    sfHists['NJetsInv'] = sfNJetsFile.Get("NJetsNoPhoCorrectionScaleFactors")
    sfVars = { "WJets":("MR","Rsq"), "TTJets":("MR","Rsq"), "DYJetsInv":("MR_NoZ","Rsq_NoZ") }
    outfile = rt.TFile(
        "data/ScaleFactors/RazorMADD2015/RazorDYJetsDileptonInvCrossCheck_%s.root"%(tag), "RECREATE")
    
    updateNorm = True
    for region in regionsOrder:
        analysis = regions[region]
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
                sfHists["DYJetsInv"].SetBinError( bn, sfHists["DYJetsInv"].GetBinError(bn)*normUpdate )
                sfHists["DYJetsInv"].SetBinError( bn, 
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

