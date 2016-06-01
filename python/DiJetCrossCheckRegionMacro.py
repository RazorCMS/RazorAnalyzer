import sys
import argparse
import copy
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import wjetsSingleLeptonCutsMC, wjetsSingleLeptonCutsData, ttjetsSingleLeptonCutsMC, ttjetsSingleLeptonCutsData, cutsMultiJetForVetoLepton, cutsMultiJetForVetoTau, cutsDiJetForVetoLepton, cutsDiJetForVetoTau, razorCuts
from macro.razorMacros import *
from SidebandMacro import LUMI as LUMI_DATA
from CrossCheckRegionMacro import MCLUMI, SAMPLES_VetoLepton, SAMPLES_MultiJet, SAMPLES_VetoTau, SAMPLES_DYJ2L_INV, ScaleFactorNames_DYJ2L_INV, ScaleFactorVars_DYJ2L_INV, FILENAMES_VetoLepton, FILENAMES_VetoTau, FILENAMES_2L_INV, FILENAMES_MULTIJET, weightOpts, VetoLeptonControlRegionBinning, MultiJetControlRegionBinning, VetoTauControlRegionBinning, MultiJetTauControlRegionBinning, ZNuNu_2L_ControlRegionBinning

printdir="DiJetCrossCheckRegionPlots"

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

    #initialize
    weightHists = {}
    sfHists = {}
    plotOpts = { "comment":False, 'SUS15004CR':True } 

    #make output directory
    os.system('mkdir -p '+printdir)

    sfHists = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root", processNames=SAMPLES_VetoLepton, scaleFactorNames={ "ZInv":"GJetsInv" }, debugLevel=debugLevel)
    sfNJetsFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors.root")
    sfHists['NJets'] = sfNJetsFile.Get("NJetsCorrectionScaleFactors")
    sfVars = ("MR","Rsq")
    sfVarsDYJetsDileptonInv = ("MR_NoZ", "Rsq_NoZ")
    auxSFs = {"NJets":("NJets40","1")} 
    auxSFsSR = {"NJets":("nSelectedJets","1")} 

    # #########################################################
    # #Veto Lepton cross-check region
    # #########################################################

    #load the MT cut efficiency as a function of lepton pt
    #mtLepFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoLeptonMTCutEfficiency.root")
    #assert mtLepFile
    #mtLepPtHist = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsPt")
    #assert mtLepPtHist
    #mtLepEtaHist = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsEta")
    #assert mtLepEtaHist

    ##load the dPhi cut efficiency as a function of lepton pt
    #dphiLepFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostLepton.root")
    #assert dphiLepFile
    #dphiLepPtHist = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsPt")
    #assert dphiLepPtHist
    #dphiLepEtaHist = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsEta")
    #assert dphiLepEtaHist

    #use these histograms to derive the additive veto lepton correction
    #vetoLeptonDiJetHists = makeControlSampleHists("VetoLeptonDiJetControlRegion", 
    #           filenames=FILENAMES_VetoLepton, samples=SAMPLES_VetoLepton, 
    #           cutsMC=vetoLeptonDiJetControlRegionCutsMC, cutsData=vetoLeptonDiJetControlRegionCutsData, 
    #           bins=VetoLeptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, auxSFs=auxSFs,
    #           weightHists=weightHists, plotDensity=True, sfHists=sfHists, weightOpts=weightOpts, 
    #           printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)
    #macro.exportHists(vetoLeptonDiJetHists, outFileName='controlHistogramsVetoLeptonDiJetControlRegion.root', 
    #        outDir=printdir, debugLevel=debugLevel)

    ##use these histograms to convert the additive veto lepton correction into a multiplicative one
    #dijetHistsForVetoLeptonCorrection = makeControlSampleHists("DiJetForVetoLeptonCorrection", 
    #        filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #        cutsMC=cutsDiJetForVetoLepton, cutsData=cutsDiJetForVetoLepton, 
    #        bins=MultiJetControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, auxSFs=auxSFsSR,
    #        weightHists=weightHists, plotDensity=False, sfHists=sfHists, treeName="RazorInclusive", 
    #        weightOpts=weightOpts, debugLevel=debugLevel, plotOpts=plotOpts, printdir=printdir)
    #macro.exportHists(dijetHistsForVetoLeptonCorrection, outDir=printdir, outFileName='controlHistogramsDiJetForVetoLeptonCorrection.root')

    #Make pt correction (control region data/MC)
    #vetoLeptonDiJetHists = macro.importHists(printdir+'/controlHistogramsVetoLeptonDiJetControlRegion.root', debugLevel)
    #dijetHistsForVetoLeptonCorrection = macro.importHists(printdir+'/controlHistogramsDiJetForVetoLeptonCorrection.root', debugLevel)
    #sfHistsCRPtCorr = sfHists.copy()
    #sfHistsCRPtCorr["VetoLepton"] = makeVetoLeptonCorrectionHist(vetoLeptonDiJetHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Lepton DiJet Pt CR", doDataOverMC=True, sfHists=sfHists, printdir=printdir)

    ##Make pt correction (signal region data/MC)
    #sfHistsSRPtCorr = sfHists.copy()
    #sfHistsSRPtCorr["VetoLepton"] = makeVetoLeptonCorrectionHist(vetoLeptonDiJetHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Lepton DiJet Pt", doDataOverMC=False, histsToCorrect=dijetHistsForVetoLeptonCorrection, signalRegionVar="leadingGenLeptonPt", mtEfficiencyHist=mtLepPtHist, dPhiEfficiencyHist=dphiLepPtHist, sfHists=sfHists, printdir=printdir)

    ##for applying pt correction
    #auxSFsCRForVetoLeptonPtCorr = auxSFs.copy()
    #auxSFsCRForVetoLeptonPtCorr['VetoLepton'] = ("lep1.Pt()", "abs(lep1Type) == 11 || abs(lep1Type) == 13")
    #auxSFsSRForVetoLeptonPtCorr = auxSFsSR.copy()
    #auxSFsSRForVetoLeptonPtCorr['VetoLepton'] = ("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13")

    #apply the pt correction and examine the residual eta agreement (control region)
    #vetoLeptonDiJetHistsPtCorr = makeControlSampleHists("VetoLeptonDiJetControlRegionForEtaCorrection", 
    #           filenames=FILENAMES_VetoLepton, samples=SAMPLES_VetoLepton, 
    #           cutsMC=vetoLeptonDiJetControlRegionCutsMC, cutsData=vetoLeptonDiJetControlRegionCutsData, 
    #           bins=VetoLeptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, plotDensity=True, sfHists=sfHistsCRPtCorr, weightOpts=weightOpts, 
    #           auxSFs=auxSFsCRForVetoLeptonPtCorr, printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)
    #macro.exportHists(vetoLeptonDiJetHistsPtCorr, outDir=printdir, outFileName='controlHistogramsVetoLeptonDiJetPtCorrControlRegion.root', debugLevel=debugLevel)

    #apply the pt correction and examine the residual eta agreement (signal region)
    #dijetHistsPtCorrForVetoLeptonCorrection = makeControlSampleHists("DiJetForVetoLeptonEtaCorrection", 
    #        filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #        cutsMC=cutsDiJetForVetoLepton, cutsData=cutsDiJetForVetoLepton, 
    #        bins=MultiJetControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #        weightHists=weightHists, plotDensity=False, sfHists=sfHistsSRPtCorr, treeName="RazorInclusive", 
    #        auxSFs=auxSFsSRForVetoLeptonPtCorr, weightOpts=weightOpts, debugLevel=debugLevel, 
    #        plotOpts=plotOpts, printdir=printdir)
    #macro.exportHists(dijetHistsPtCorrForVetoLeptonCorrection, outDir=printdir, outFileName='controlHistogramsDiJetPtCorrForVetoLeptonCorrection.root', debugLevel=debugLevel)

    #Make residual eta correction (control region data/MC)
    #vetoLeptonDiJetHistsPtCorr = macro.importHists(printdir+'/controlHistogramsVetoLeptonDiJetPtCorrControlRegion.root', debugLevel)
    #makeVetoLeptonCorrectionHist(vetoLeptonDiJetHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Lepton DiJet Eta CR", doDataOverMC=True, sfHists=sfHists, printdir=printdir)

    ##Make residual eta correction (signal region data/MC)
    #dijetHistsPtCorrForVetoLeptonCorrection = macro.importHists(printdir+'/controlHistogramsDiJetPtCorrForVetoLeptonCorrection.root')
    #makeVetoLeptonCorrectionHist(vetoLeptonDiJetHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Lepton DiJet Eta", doDataOverMC=False, histsToCorrect=dijetHistsPtCorrForVetoLeptonCorrection, signalRegionVar="abs(leadingGenLeptonEta)", mtEfficiencyHist=mtLepEtaHist, dPhiEfficiencyHist=dphiLepEtaHist, sfHists=sfHists, printdir=printdir)

    ##########################################################
    ##Veto Tau cross-check region
    ##########################################################

    ##load the MT cut efficiency as a function of tau pt
    #mtTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoTauMTCutEfficiency.root")
    #assert mtTauFile
    #mtTauPtHist = mtTauFile.Get("VetoTauMTCutEfficiencyVsPt")
    #assert mtTauPtHist
    #mtTauEtaHist = mtTauFile.Get("VetoTauMTCutEfficiencyVsEta")
    #assert mtTauEtaHist

    ##load the dPhi cut efficiency as a function of tau pt
    #dphiTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostTau.root")
    #assert dphiTauFile
    #dphiTauPtHist = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsPt")
    #assert dphiTauPtHist
    #dphiTauEtaHist = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsEta")
    #assert dphiTauEtaHist

    ##use these histograms to derive the additive veto tau correction
    #vetoTauDiJetHists = makeControlSampleHists("VetoTauDiJetControlRegion", 
    #           filenames=FILENAMES_VetoTau, samples=SAMPLES_VetoTau, 
    #           cutsMC=vetoTauDiJetControlRegionCutsMC, cutsData=vetoTauDiJetControlRegionCutsData, 
    #           bins=VetoTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, auxSFs=auxSFs,
    #           weightHists=weightHists, plotDensity=True, sfHists=sfHists, weightOpts=weightOpts, 
    #           printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)
    #macro.exportHists(vetoTauDiJetHists, outFileName='controlHistogramsVetoTauDiJetControlRegion.root', 
    #        outDir=printdir, debugLevel=debugLevel)

    ##use these histograms to convert the additive veto tau correction into a multiplicative one
    #dijetHistsForVetoTauCorrection = makeControlSampleHists("DiJetForVetoTauCorrection", 
    #        filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #        cutsMC=cutsDiJetForVetoTau, cutsData=cutsDiJetForVetoTau, 
    #        bins=MultiJetTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, auxSFs=auxSFsSR,
    #        weightHists=weightHists, plotDensity=False, sfHists=sfHists, treeName="RazorInclusive", 
    #        weightOpts=weightOpts, debugLevel=debugLevel, plotOpts=plotOpts, printdir=printdir)
    #macro.exportHists(dijetHistsForVetoTauCorrection, outDir=printdir, outFileName='controlHistogramsDiJetForVetoTauCorrection.root')

    ##Make pt correction (control region data/MC)
    #vetoTauDiJetHists = macro.importHists(printdir+'/controlHistogramsVetoTauDiJetControlRegion.root', debugLevel)
    #dijetHistsForVetoTauCorrection = macro.importHists(printdir+'/controlHistogramsDiJetForVetoTauCorrection.root', debugLevel)
    #sfHistsCRPtCorr = sfHists.copy()
    #sfHistsCRPtCorr["VetoTau"] = makeVetoLeptonCorrectionHist(vetoTauDiJetHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Tau DiJet Pt CR", doDataOverMC=True, sfHists=sfHists, printdir=printdir)

    ##Make pt correction (signal region data/MC)
    #sfHistsSRPtCorr = sfHists.copy()
    #sfHistsSRPtCorr["VetoTau"] = makeVetoLeptonCorrectionHist(vetoTauDiJetHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Tau DiJet Pt", doDataOverMC=False, histsToCorrect=dijetHistsForVetoTauCorrection, signalRegionVar="leadingGenLeptonPt", mtEfficiencyHist=mtTauPtHist, dPhiEfficiencyHist=dphiTauPtHist, sfHists=sfHists, printdir=printdir)

    ##for applying pt correction
    #auxSFsCRForVetoTauPtCorr = auxSFs.copy()
    #auxSFsCRForVetoTauPtCorr['VetoTau'] = ("lep1.Pt()", "abs(lep1Type) == 15")
    #auxSFsSRForVetoTauPtCorr = auxSFsSR.copy()
    #auxSFsSRForVetoTauPtCorr['VetoTau'] = ("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 15")

    ##apply the pt correction and examine the residual eta agreement (control region)
    #vetoTauDiJetHistsPtCorr = makeControlSampleHists("VetoTauDiJetControlRegionForEtaCorrection", 
    #           filenames=FILENAMES_VetoTau, samples=SAMPLES_VetoTau, 
    #           cutsMC=vetoTauDiJetControlRegionCutsMC, cutsData=vetoTauDiJetControlRegionCutsData, 
    #           bins=VetoTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, plotDensity=True, sfHists=sfHistsCRPtCorr, weightOpts=weightOpts, 
    #           auxSFs=auxSFsCRForVetoTauPtCorr, printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)
    #macro.exportHists(vetoTauDiJetHistsPtCorr, outDir=printdir, outFileName='controlHistogramsVetoTauDiJetPtCorrControlRegion.root', debugLevel=debugLevel)

    ##apply the pt correction and examine the residual eta agreement (signal region)
    #dijetHistsPtCorrForVetoTauCorrection = makeControlSampleHists("DiJetForVetoTauEtaCorrection", 
    #        filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #        cutsMC=cutsDiJetForVetoTau, cutsData=cutsDiJetForVetoTau, 
    #        bins=MultiJetTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #        weightHists=weightHists, plotDensity=False, sfHists=sfHistsSRPtCorr, treeName="RazorInclusive", 
    #        auxSFs=auxSFsSRForVetoTauPtCorr, weightOpts=weightOpts, debugLevel=debugLevel, 
    #        plotOpts=plotOpts, printdir=printdir)
    #macro.exportHists(dijetHistsPtCorrForVetoTauCorrection, outDir=printdir, outFileName='controlHistogramsDiJetPtCorrForVetoTauCorrection.root', debugLevel=debugLevel)

    ##Make residual eta correction (control region data/MC)
    #vetoTauDiJetHistsPtCorr = macro.importHists(printdir+'/controlHistogramsVetoTauDiJetPtCorrControlRegion.root', debugLevel)
    #makeVetoLeptonCorrectionHist(vetoTauDiJetHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Tau DiJet Eta CR", doDataOverMC=True, sfHists=sfHists, printdir=printdir)

    ##Make residual eta correction (signal region data/MC)
    #dijetHistsPtCorrForVetoTauCorrection = macro.importHists(printdir+'/controlHistogramsDiJetPtCorrForVetoTauCorrection.root')
    #makeVetoLeptonCorrectionHist(vetoTauDiJetHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Tau DiJet Eta", doDataOverMC=False, histsToCorrect=dijetHistsPtCorrForVetoTauCorrection, signalRegionVar="abs(leadingGenLeptonEta)", mtEfficiencyHist=mtTauEtaHist, dPhiEfficiencyHist=dphiTauEtaHist, sfHists=sfHists, printdir=printdir)

    ##########################################################
    #Z->LL dilepton control sample
    ##########################################################
    sfHistsDileptonInv = loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root", 
            processNames=SAMPLES_DYJ2L_INV, scaleFactorNames=ScaleFactorNames_DYJ2L_INV, debugLevel=debugLevel)
    sfNJetsFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors.root")
    sfHistsDileptonInv['NJets'] = sfNJetsFile.Get("NJetsCorrectionScaleFactors")
    auxSFsDileptonInv = {"NJets":("NJets_NoZ","1")} 
    dyjetsDileptonInvDiJetHists = makeControlSampleHists("DYJetsDileptonInvDiJet", 
                filenames=FILENAMES_2L_INV, samples=SAMPLES_DYJ2L_INV, 
                cutsMC=dyjetsDileptonInvDiJetCutsMC, cutsData=dyjetsDileptonInvDiJetCutsData, 
                bins=ZNuNu_2L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
                weightHists=weightHists, plotDensity=True, sfHists=sfHistsDileptonInv, 
                weightOpts=weightOpts, sfVars=ScaleFactorVars_DYJ2L_INV, auxSFs=auxSFsDileptonInv,
                printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)

    #Record discrepancies > 1 sigma
    tmpSFHists = copy.copy(sfHistsDileptonInv)
    del tmpSFHists["DYJetsInv"]
    appendScaleFactors("DYJetsInv", dyjetsDileptonInvDiJetHists, tmpSFHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var=sfVarsDYJetsDileptonInv, signifThreshold=1.0, printdir=printdir)

    #write DYJetsInv cross check scale factors
    dyjetsDileptonInvOutfile = rt.TFile("RazorDYJetsDileptonInvDiJetCrossCheck.root", "RECREATE")
    print "Writing histogram",tmpSFHists["DYJetsInv"].GetName(),"to file"
    tmpSFHists["DYJetsInv"].Write("DYJetsDileptonInvDiJetCrossCheckScaleFactors")
    dyjetsDileptonInvOutfile.Close()

