import sys, os, argparse, copy
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors, makeVetoLeptonCorrectionHist

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

    #load the MT cut efficiency as a function of lepton pt
    mtHists = {}
    mtLepFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoLeptonMTCutEfficiency.root")
    mtHists["VetoLeptonPt"] = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsPt")
    mtHists["VetoLeptonEta"] = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsEta")
    #load the MT cut efficiency as a function of tau pt
    mtTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoTauMTCutEfficiency.root")
    mtHists["VetoTauPt"] = mtTauFile.Get("VetoTauMTCutEfficiencyVsPt")
    mtHists["VetoTauEta"] = mtTauFile.Get("VetoTauMTCutEfficiencyVsEta")
    #load the dPhi cut efficiency as a function of lepton pt
    dphiHists = {}
    dphiLepFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostLepton.root")
    dphiHists["VetoLeptonPt"] = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsPt")
    dphiHists["VetoLeptonEta"] = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsEta")
    #load the dPhi cut efficiency as a function of tau pt
    dphiTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostTau.root")
    dphiHists["VetoTauPt"] = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsPt")
    dphiHists["VetoTauEta"] = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsEta")

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regionsOrder = []
    regions = {}
    regionsCorrespondence = {} #get the control region corresponding to each signal region
    regionMtHists = {}
    regionDphiHists = {}
    for ltype in ["VetoLepton","VetoTau"]:
        for jtype,jets in {"DiJet":(2,3),"MultiJet":(4,-1)}.iteritems():
            #veto lepton/tau control region
            regionsOrder.append(ltype+jtype) 
            regions[ltype+jtype] = Analysis(ltype+"ControlRegion",tag=tag,
                    njetsMin=jets[0],njetsMax=jets[1])
            #corresponding signal region 
            sigR = jtype+"For"+ltype
            regionsOrder.append(sigR) 
            regions[jtype+"For"+ltype] = Analysis(sigR+"ControlRegion",
                    tag=tag,njetsMin=jets[0],njetsMax=jets[1])
            regionsCorrespondence[sigR] = ltype+jtype
            regionMtHists[sigR] = mtHists[ltype+"Pt"]
            regionDphiHists[sigR] = dphiHists[ltype+"Pt"]
            #veto lepton/tau control region
            regionsOrder.append(ltype+jtype+"PtCorr") 
            regions[ltype+jtype+"PtCorr"] = Analysis(ltype+"ControlRegion",tag=tag,
                    njetsMin=jets[0],njetsMax=jets[1])
            #corresponding signal region 
            sigRPtCorr = jtype+"For"+ltype+"PtCorr"
            regionsOrder.append(sigRPtCorr) 
            regions[sigRPtCorr] = Analysis(sigR+"ControlRegion",
                    tag=tag,njetsMin=jets[0],njetsMax=jets[1])
            regionsCorrespondence[sigRPtCorr] = ltype+jtype+"PtCorr"
            regionMtHists[sigRPtCorr] = mtHists[ltype+"Eta"]
            regionDphiHists[sigRPtCorr] = dphiHists[ltype+"Eta"]

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions["VetoLeptonDiJet"].samples, scaleFactorNames={ "ZInv":"GJetsInv" },
            debugLevel=debugLevel)
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJets'] = sfNJetsFile.Get("NJetsCorrectionScaleFactors")
    sfHists['NJetsInv'] = sfNJetsFile.Get("GJetsScaleFactorVsNJets")
    sfVars = ("MR","Rsq")
    outfile = rt.TFile("data/ScaleFactors/RazorMADD2015/RazorVetoLeptonClosureTests_%s.root"%(tag), "RECREATE")
    outfile.Close()

    hists = {}
    for region in regionsOrder:
        analysis = regions[region]
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #set up analysis
        if region.startswith('Veto'):
            auxSFsToUse = razorWeights.getNJetsSFs(analysis,jetName='NJets40')
            treeName = "ControlSampleEvent"
        else:
            auxSFsToUse = razorWeights.getNJetsSFs(analysis,jetName='nSelectedJets')
            treeName = "RazorInclusive"
        #set up lepton pt correction
        varForCorrection = "lep1.Pt()"
        sigVarForCorrection = "leadingGenLeptonPt"
        if "PtCorr" in region: 
            varForCorrection = "abs(lep1.Eta())"
            sigVarForCorrection = "abs(leadingGenLeptonEta)"
            if region.startswith('VetoLepton'):
                auxSFsToUse[region.replace('PtCorr','')] = ("lep1.Pt()", 
                        "abs(lep1Type) == 11 || abs(lep1Type) == 13")
            elif region.startswith('VetoTau'):
                auxSFsToUse[region.replace('PtCorr','')] = ("lep1.Pt()", "abs(lep1Type) == 15")
            elif 'VetoLepton' in region:
                auxSFsToUse[region.replace('PtCorr','')] = ("leadingGenLeptonPt", 
                        "abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13")
            elif 'VetoTau' in region:
                auxSFsToUse[region.replace('PtCorr','')] = ("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 15")
        #sanity check
        print "\nRegion:",region
        print "Tree name:",treeName
        print "Aux SFs to use:",auxSFsToUse
        print "Variable for correction:",varForCorrection
        print "Signal region variable for correction:",sigVarForCorrection,"\n"
        #perform analysis
        hists[region] = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                sfHists=sfHists, sfVars=sfVars, printdir=outdir, auxSFs=auxSFsToUse, 
                treeName=treeName, debugLevel=debugLevel )
        #export histograms
        macro.exportHists(hists[region], outFileName='controlHistograms'+region+'.root', 
                outDir=outdir, debugLevel=debugLevel, delete=False)
        #compute correction factors
        if region.startswith('Veto'):
            #make control region scale factors
            sfHists[region] = makeVetoLeptonCorrectionHist(hists[region], lumiData=analysis.lumi, 
                    debugLevel=debugLevel, var=varForCorrection, signifThreshold=1.0, 
                    regionName=region, doDataOverMC=True, sfHists=sfHists, 
                    printdir=outdir)
        else:
            #make signal region scale factors
            controlRegionHists = hists[regionsCorrespondence[region]]
            mtHistToUse = regionMtHists[region]
            dphiHistToUse = regionDphiHists[region]
            sfHists[region] = makeVetoLeptonCorrectionHist(controlRegionHists, 
                    lumiData=analysis.lumi, debugLevel=debugLevel, var=varForCorrection, 
                    signifThreshold=1.0, regionName=region, doDataOverMC=False, 
                    histsToCorrect=hists[region], signalRegionVar=sigVarForCorrection, 
                    mtEfficiencyHist=mtHistToUse, dPhiEfficiencyHist=dphiHistToUse, 
                    sfHists=sfHists, printdir=outdir)
            #write out to file
            sfHistClone = sfHists[region].Clone()
            print "Writing scale factor histogram",sfHistClone.GetName(),"to file"
            outfile = rt.TFile("data/ScaleFactors/RazorMADD2015/RazorVetoLeptonClosureTests_%s.root"%(tag), "UPDATE")
            outfile.Append(sfHistClone)
            outfile.Write()
            outfile.Close()

