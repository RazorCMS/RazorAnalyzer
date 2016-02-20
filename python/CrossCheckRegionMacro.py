import sys
import argparse
import copy
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import wjetsSingleLeptonCutsMC, wjetsSingleLeptonCutsData, ttjetsSingleLeptonCutsMC, ttjetsSingleLeptonCutsData, cutsMultiJetForVetoLepton, cutsMultiJetForVetoTau, razorCuts
from macro.razorMacros import *
from SidebandMacro import LUMI as LUMI_DATA

MCLUMI = 1 

SAMPLES_TTJ1L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["Other", "DYJets", "SingleTop", "TTJets", "WJets"]
SAMPLES_TTJ2L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_DYJ2L = ["Other", "SingleTop", "WJets", "TTJets", "DYJets"]
SAMPLES_VetoLepton = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_MultiJet = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets1L", "TTJets2L"]
SAMPLES_VetoTau = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_DYJ2L_INV = ["Other", "SingleTop", "WJets", "TTJets", "DYJetsInv"]
ScaleFactorNames_DYJ2L_INV = {"Other"     : "Other", 
                              "SingleTop" : "SingleTop", 
                              "WJets"     : "WJets", 
                              "TTJets"    :"TTJets" , 
                              "DYJetsInv"    :"WJetsInv"
                              }
#apply correct version of MR when applying scale factors
ScaleFactorVars_DYJ2L_INV = { "WJets":("MR","Rsq"),
                              "TTJets":("MR","Rsq"),
                              "DYJetsInv":("MR_NoZ","Rsq_NoZ"),
                              }

DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonFull_1p23_2015Final/RazorSkim/"
PREFIX_1L = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
FILENAMES_1L = {
            "TTJets"   : DIR_1L+"/"+PREFIX_1L+"_TTJets_1pb_weighted_RazorSkim.root",
            "WJets"    : DIR_1L+"/"+PREFIX_1L+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_1L+"/"+PREFIX_1L+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_1L+"/"+PREFIX_1L+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_1L+"/"+PREFIX_1L+"_Other_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_1L+"/"+PREFIX_1L+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }

DIR_2L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RazorSkim/"
PREFIX_2L = "RunTwoRazorControlRegions_DileptonFull_DileptonSkim"
FILENAMES_2L = {
            "TTJets"   : DIR_2L+"/"+PREFIX_2L+"_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root",
            "WJets"    : DIR_2L+"/"+PREFIX_2L+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_2L+"/"+PREFIX_2L+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_2L+"/"+PREFIX_2L+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_2L+"/"+PREFIX_2L+"_VV_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_2L+"/"+PREFIX_2L+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }


DIR_VetoLepton = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLeptonFull_1p23_2015Final/RazorNJets80Skim/"
PREFIX_VetoLepton = "RunTwoRazorControlRegions_VetoLeptonFull"
FILENAMES_VetoLepton = {
            "TTJets"   : DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_TTJets_1pb_weighted_RazorSkim.root",
            "WJets"    : DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "QCD"      : DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",
            "ZInv"     : DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",            
            "Other"    : DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_Other_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_VetoLepton+"/"+PREFIX_VetoLepton+"_HTMHT_Run2015D_GoodLumiGolden_RazorSkim.root"
            }

DIR_VetoTau = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoTauFull_1p23_2015Final/"
PREFIX_VetoTau = "RunTwoRazorControlRegions_VetoTauFull_RazorSkim"
FILENAMES_VetoTau = {
            "TTJets"   : DIR_VetoTau+"/"+PREFIX_VetoTau+"_TTJets_1pb_weighted.root",
            "WJets"    : DIR_VetoTau+"/"+PREFIX_VetoTau+"_WJetsToLNu_HTBinned_1pb_weighted.root",
            "SingleTop": DIR_VetoTau+"/"+PREFIX_VetoTau+"_SingleTop_1pb_weighted.root",
            "DYJets"   : DIR_VetoTau+"/"+PREFIX_VetoTau+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted.root",
            "QCD"      : DIR_VetoTau+"/"+PREFIX_VetoTau+"_QCD_HTBinned_1pb_weighted.root",
            "ZInv"     : DIR_VetoTau+"/"+PREFIX_VetoTau+"_ZJetsToNuNu_HTBinned_1pb_weighted.root",            
            "Other"    : DIR_VetoTau+"/"+PREFIX_VetoTau+"_Other_1pb_weighted.root",
            "Data"     : DIR_VetoTau+"/"+PREFIX_VetoTau+"_HTMHT_Run2015D_GoodLumiGolden.root"
            }

DIR_2L_INV = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFullAddToMET_1p23_2015Final/RazorNJets80Skim/"
PREFIX_2L_INV = "RunTwoRazorControlRegions_DileptonAddToMetFull_DileptonSkim"
FILENAMES_2L_INV = {
            "TTJets"   : DIR_2L_INV+"/"+PREFIX_2L_INV+"_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root",
            "WJets"    : DIR_2L_INV+"/"+PREFIX_2L_INV+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_2L_INV+"/"+PREFIX_2L_INV+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJetsInv"   : DIR_2L_INV+"/"+PREFIX_2L_INV+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_2L_INV+"/"+PREFIX_2L_INV+"_Other_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_2L_INV+"/"+PREFIX_2L_INV+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }

DIR_MULTIJET = "Backgrounds"
#DIR_MULTIJET = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108/"
FILENAMES_MULTIJET = {
        "TTJets1L"    : DIR_MULTIJET+"/"+"FullRazorInclusive_TTJets1L_1pb_weighted.root",
        "TTJets2L"    : DIR_MULTIJET+"/"+"FullRazorInclusive_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
        "WJets"     : DIR_MULTIJET+"/"+"FullRazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop" : DIR_MULTIJET+"/"+"FullRazorInclusive_SingleTop_1pb_weighted.root",
        "Other" : DIR_MULTIJET+"/"+"FullRazorInclusive_Other_1pb_weighted.root",
        "DYJets"     : DIR_MULTIJET+"/"+"FullRazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted.root",
        "ZInv"     : DIR_MULTIJET+"/"+"FullRazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root",
        "QCD"      : DIR_MULTIJET+'/'+'FullRazorInclusive_QCD_HTBinned_1pb_weighted.root',
        }

weightOpts = []

config = "config/run2_20151229_ControlRegion.config"
cfg = Config.Config(config)

VetoLeptonBinsMRLep = cfg.getBinning("VetoLeptonControlRegion")[0]
VetoLeptonBinsRsqLep = cfg.getBinning("VetoLeptonControlRegion")[1]
VetoLeptonBinsLepPt = [5, 10, 15, 20.,30.,40.,100,1000]
VetoLeptonBinsLepEta = [0, 0.5, 1.0, 1.5, 2.0, 2.5]
VetoLeptonControlRegionBinning = { "MR":VetoLeptonBinsMRLep, "Rsq":VetoLeptonBinsRsqLep, "lep1.Pt()":VetoLeptonBinsLepPt , "abs(lep1.Eta())":VetoLeptonBinsLepEta, ("MR","Rsq"):[], ("abs(lep1.Eta())","lep1.Pt()"):[]}
MultiJetControlRegionBinning = { "MR":VetoLeptonBinsMRLep, "Rsq":VetoLeptonBinsRsqLep, "leadingGenLeptonPt":VetoLeptonBinsLepPt, "abs(leadingGenLeptonEta)":VetoLeptonBinsLepEta, ("MR","Rsq"):[] }

VetoLepton7JetBinsMRLep = [400, 700, 900, 4000]
VetoLepton7JetBinsRsqLep =  [0.25,0.30,0.41,1.5]
VetoLepton7JetBinsLepPt = [5, 15, 40.,1000]
VetoLepton7JetBinsLepEta = [0, 1.5, 2.5]
VetoLeptonControlRegion7JetBinning = { "MR":VetoLepton7JetBinsMRLep, "Rsq":VetoLepton7JetBinsRsqLep, "lep1.Pt()":VetoLepton7JetBinsLepPt , "abs(lep1.Eta())":VetoLepton7JetBinsLepEta, ("MR","Rsq"):[], ("abs(lep1.Eta())","lep1.Pt()"):[]}


VetoTauBinsLepPt = [20,30,40,100,1000]
VetoTauBinsLepEta = [0, 0.5, 1.0, 1.5, 2.0, 2.5]
VetoTauBinsMR =  [400, 500, 600, 700, 900, 4000]
VetoTauBinsRsq = [0.25,0.30,0.41,1.5]
VetoTauControlRegionBinning = { "MR":VetoTauBinsMR, "Rsq":VetoTauBinsRsq, "lep1.Pt()":VetoTauBinsLepPt , "abs(lep1.Eta())":VetoTauBinsLepEta, ("MR","Rsq"):[], ("abs(lep1.Eta())","lep1.Pt()"):[]}
MultiJetTauControlRegionBinning = { "MR":VetoTauBinsMR, "Rsq":VetoTauBinsRsq, "leadingGenLeptonPt":VetoTauBinsLepPt, "abs(leadingGenLeptonEta)":VetoTauBinsLepEta, ("MR","Rsq"):[] }

VetoTau7JetBinsMRLep = [400, 700, 900, 4000]
VetoTau7JetBinsRsqLep =  [0.25,0.30,0.41,1.5]
VetoTau7JetBinsLepPt = [5, 10, 15, 20.,30.,40.,100,1000]
VetoTau7JetBinsLepEta = [0, 1.5, 2.5]
VetoTauControlRegion7JetBinning = { "MR":VetoTau7JetBinsMRLep, "Rsq":VetoTau7JetBinsRsqLep, "lep1.Pt()":VetoTau7JetBinsLepPt , "abs(lep1.Eta())":VetoTau7JetBinsLepEta, ("MR","Rsq"):[], ("abs(lep1.Eta())","lep1.Pt()"):[]}

TTJetsDileptonBinsMRLep = cfg.getBinning("TTJetsDileptonControlRegion")[0]
TTJetsDileptonBinsRsqLep = cfg.getBinning("TTJetsDileptonControlRegion")[1]
TTJetsDileptonBinsNBTags = [0.,1.,2.,3.,4.]
TTJetsDileptonBinsNJets80 = [0.,1.,2.,3.,4.]
TTJetsDileptonBinsNJets = [0.,1.,2.,3.,4.,5.,6.,7.]
TTJetsDileptonBinsLepPt = [20.,25.,30.,35.,40.,45.,50.,70.,100]
TTJetsDileptonControlRegionBinning = { "MR":TTJetsDileptonBinsMRLep, "Rsq":TTJetsDileptonBinsRsqLep, "NBJetsMedium":TTJetsDileptonBinsNBTags, "NJets80":TTJetsDileptonBinsNJets80, "NJets40":TTJetsDileptonBinsNJets, ("MR","Rsq"):[]}

ZNuNu_2L_4Jet_BinsMRLep = [400, 500, 600, 700, 900, 4000]
ZNuNu_2L_4Jet_BinsRsqLep = [0.25,0.30,0.41,1.5]
ZNuNu_2L_BinsMRLep = [400, 500, 600, 700, 900, 1200, 4000]
ZNuNu_2L_BinsRsqLep = [0.25,0.30,0.41,0.52,0.64,1.5]
ZNuNu_2L_BinsNBTags = [0,1,2,3,4]
ZNuNu_2L_BinsNJets80 = [0,1,2,3,4]
ZNuNu_2L_BinsNJets = [0,1,2,3,4,5,6,7,8]
ZNuNu_2L_ControlRegionBinning = { "MR_NoZ":ZNuNu_2L_BinsMRLep, "Rsq_NoZ":ZNuNu_2L_BinsRsqLep, "NBJetsMedium":ZNuNu_2L_BinsNBTags, "NJets80":ZNuNu_2L_BinsNJets80, "NJets40":ZNuNu_2L_BinsNJets, ("MR_NoZ","Rsq_NoZ"):[] }
ZNuNu_2L_4Jet_ControlRegionBinning = { "MR_NoZ":ZNuNu_2L_4Jet_BinsMRLep, "Rsq_NoZ":ZNuNu_2L_4Jet_BinsRsqLep, "NBJetsMedium":ZNuNu_2L_BinsNBTags, "NJets80_NoZ":ZNuNu_2L_BinsNJets80, "NJets_NoZ":ZNuNu_2L_BinsNJets, ("MR_NoZ","Rsq_NoZ"):[] }

printdir="CrossCheckRegionPlots"

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
    plotOpts = { "comment":False } 

    #make output directory
    os.system('mkdir -p '+printdir)

    sfHists = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root", processNames=SAMPLES_VetoLepton, scaleFactorNames={ "ZInv":"GJetsInv" }, debugLevel=debugLevel)
    sfVars = ("MR","Rsq")
    sfVarsDYJetsDileptonInv = ("MR_NoZ", "Rsq_NoZ")

    ##########################################################
    #TTJets dilepton control sample
    ##########################################################
    #ttjetsDileptonHists = makeControlSampleHists("TTJetsDilepton", 
    #           filenames=FILENAMES_2L, samples=SAMPLES_TTJ2L, 
    #           cutsMC=ttjetsDileptonCutsMC, cutsData=ttjetsDileptonCutsData, 
    #           bins=TTJetsDileptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, plotDensity=True, sfHists=sfHists, weightOpts=weightOpts, 
    #           printdir=printdir, debugLevel=debugLevel)

    ##Record discrepancies > 1 sigma
    #tmpSFHists = copy.copy(sfHists)
    #del tmpSFHists["TTJets"]
    #appendScaleFactors("TTJets", ttjetsDileptonHists, tmpSFHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var=sfVars, signifThreshold=1.0, printdir=printdir)

    ##write TTJets cross check scale factors
    #ttjetsDileptonOutfile = rt.TFile("RazorTTJetsDileptonCrossCheck.root", "RECREATE")
    #print "Writing histogram",tmpSFHists["TTJets"].GetName(),"to file"
    #tmpSFHists["TTJets"].Write("TTJetsDileptonCrossCheckScaleFactors")
    #ttjetsDileptonOutfile.Close()

    # #########################################################
    # #Veto Lepton cross-check region
    # #########################################################

    # #load the MT cut efficiency as a function of lepton pt
    # mtLepFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoLeptonMTCutEfficiency.root")
    # assert mtLepFile
    # mtLepPtHist = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsPt")
    # assert mtLepPtHist
    # mtLepEtaHist = mtLepFile.Get("VetoLeptonMTCutEfficiencyVsEta")
    # assert mtLepEtaHist

    # #load the dPhi cut efficiency as a function of lepton pt
    # dphiLepFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostLepton.root")
    # assert dphiLepFile
    # dphiLepPtHist = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsPt")
    # assert dphiLepPtHist
    # dphiLepEtaHist = dphiLepFile.Get("VetoLeptonDPhiCutEfficiencyVsEta")
    # assert dphiLepEtaHist

    # #use these histograms to derive the additive veto lepton correction
    # vetoLeptonHists = makeControlSampleHists("VetoLeptonControlRegion", 
    #            filenames=FILENAMES_VetoLepton, samples=SAMPLES_VetoLepton, 
    #            cutsMC=vetoLeptonControlRegionCutsMC, cutsData=vetoLeptonControlRegionCutsData, 
    #            bins=VetoLeptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, plotDensity=False, sfHists=sfHists, weightOpts=weightOpts, 
    #            printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)

    # #use these histograms to convert the additive veto lepton correction into a multiplicative one
    # multijetHistsForVetoLeptonCorrection = makeControlSampleHists("MultiJetForVetoLeptonCorrection", 
    #         filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #         cutsMC=cutsMultiJetForVetoLepton, cutsData=cutsMultiJetForVetoLepton, 
    #         bins=MultiJetControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #         weightHists=weightHists, plotDensity=False, sfHists=sfHists, treeName="RazorInclusive", 
    #         weightOpts=weightOpts, debugLevel=debugLevel, plotOpts=plotOpts, printdir=printdir)

    # #Make pt correction (control region data/MC)
    # sfHistsCRPtCorr = sfHists.copy()
    # sfHistsCRPtCorr["VetoLepton"] = makeVetoLeptonCorrectionHist(vetoLeptonHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Lepton Pt CR", doDataOverMC=True, printdir=printdir)

    # #Make pt correction (signal region data/MC)
    # sfHistsSRPtCorr = sfHists.copy()
    # sfHistsSRPtCorr["VetoLepton"] = makeVetoLeptonCorrectionHist(vetoLeptonHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Lepton Pt", doDataOverMC=False, histsToCorrect=multijetHistsForVetoLeptonCorrection, signalRegionVar="leadingGenLeptonPt", mtEfficiencyHist=mtLepPtHist, dPhiEfficiencyHist=dphiLepPtHist, printdir=printdir)

    # #for applying pt correction
    # auxSFsCRForVetoLeptonPtCorr = { 
    #     "VetoLepton":("lep1.Pt()", "abs(lep1Type) == 11 || abs(lep1Type) == 13"), 
    #     }
    # auxSFsSRForVetoLeptonPtCorr = { 
    #     "VetoLepton":("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13"), 
    #     }

    # #apply the pt correction and examine the residual eta agreement (control region)
    # vetoLeptonHistsPtCorr = makeControlSampleHists("VetoLeptonControlRegionForEtaCorrection", 
    #            filenames=FILENAMES_VetoLepton, samples=SAMPLES_VetoLepton, 
    #            cutsMC=vetoLeptonControlRegionCutsMC, cutsData=vetoLeptonControlRegionCutsData, 
    #            bins=VetoLeptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, plotDensity=False, sfHists=sfHistsCRPtCorr, weightOpts=weightOpts, 
    #            auxSFs=auxSFsCRForVetoLeptonPtCorr, printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)

    # #apply the pt correction and examine the residual eta agreement (signal region)
    # multijetHistsPtCorrForVetoLeptonCorrection = makeControlSampleHists("MultiJetForVetoLeptonEtaCorrection", 
    #         filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #         cutsMC=cutsMultiJetForVetoLepton, cutsData=cutsMultiJetForVetoLepton, 
    #         bins=MultiJetControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #         weightHists=weightHists, plotDensity=False, sfHists=sfHistsSRPtCorr, treeName="RazorInclusive", 
    #         auxSFs=auxSFsSRForVetoLeptonPtCorr, weightOpts=weightOpts, debugLevel=debugLevel, 
    #         plotOpts=plotOpts, printdir=printdir)

    # #Make residual eta correction (control region data/MC)
    # sfHistsCRPtEtaCorr = sfHists.copy()
    # sfHistsCRPtEtaCorr["VetoLeptonPt"] = sfHistsCRPtCorr["VetoLepton"]
    # sfHistsCRPtEtaCorr["VetoLeptonEta"] = makeVetoLeptonCorrectionHist(vetoLeptonHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Lepton Eta CR", doDataOverMC=True, printdir=printdir)

    # #Make residual eta correction (signal region data/MC)
    # makeVetoLeptonCorrectionHist(vetoLeptonHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Lepton Eta", doDataOverMC=False, histsToCorrect=multijetHistsPtCorrForVetoLeptonCorrection, signalRegionVar="abs(leadingGenLeptonEta)", mtEfficiencyHist=mtLepEtaHist, dPhiEfficiencyHist=dphiLepEtaHist, printdir=printdir)

    # auxSFsCRForVetoLeptonPtEtaCorr = { 
    #     "VetoLeptonPt":("lep1.Pt()", "abs(lep1Type) == 11 || abs(lep1Type) == 13"), 
    #     "VetoLeptonEta":("abs(lep1.Eta())", "abs(lep1Type) == 11 || abs(lep1Type) == 13"), 
    #     }

    # #Check control region after both corrections
    # vetoLeptonHistsPtEtaCorr = makeControlSampleHists("VetoLeptonControlRegionPtEtaCorrected", 
    #            filenames=FILENAMES_VetoLepton, samples=SAMPLES_VetoLepton, 
    #            cutsMC=vetoLeptonControlRegionCutsMC, cutsData=vetoLeptonControlRegionCutsData, 
    #            bins=VetoLeptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, plotDensity=True, sfHists=sfHistsCRPtEtaCorr, weightOpts=weightOpts, 
    #            auxSFs=auxSFsCRForVetoLeptonPtEtaCorr, printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)

    # ##########################################################
    # ##Veto Tau cross-check region
    # ##########################################################

    # #load the MT cut efficiency as a function of tau pt
    # mtTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/VetoTauMTCutEfficiency.root")
    # assert mtTauFile
    # mtTauPtHist = mtTauFile.Get("VetoTauMTCutEfficiencyVsPt")
    # assert mtTauPtHist
    # mtTauEtaHist = mtTauFile.Get("VetoTauMTCutEfficiencyVsEta")
    # assert mtTauEtaHist

    # #load the dPhi cut efficiency as a function of tau pt
    # dphiTauFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/DPhiCutEfficiencyForLostTau.root")
    # assert dphiTauFile
    # dphiTauPtHist = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsPt")
    # assert dphiTauPtHist
    # dphiTauEtaHist = dphiTauFile.Get("VetoTauDPhiCutEfficiencyVsEta")
    # assert dphiTauEtaHist

    # #use these histograms to derive the additive veto tau correction
    # vetoTauHists = makeControlSampleHists("VetoTauControlRegion", 
    #            filenames=FILENAMES_VetoTau, samples=SAMPLES_VetoTau, 
    #            cutsMC=vetoTauControlRegionCutsMC, cutsData=vetoTauControlRegionCutsData, 
    #            bins=VetoTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, plotDensity=False, sfHists=sfHists, weightOpts=weightOpts, 
    #            printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)

    # #use these histograms to convert the additive veto tau correction into a multiplicative one
    # multijetHistsForVetoTauCorrection = makeControlSampleHists("MultiJetForVetoTauCorrection", 
    #         filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #         cutsMC=cutsMultiJetForVetoTau, cutsData=cutsMultiJetForVetoTau, 
    #         bins=MultiJetTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #         weightHists=weightHists, plotDensity=False, sfHists=sfHists, treeName="RazorInclusive", 
    #         weightOpts=weightOpts, debugLevel=debugLevel, plotOpts=plotOpts, printdir=printdir)

    # #Make pt correction (control region data/MC)
    # sfHistsCRPtCorr = sfHists.copy()
    # sfHistsCRPtCorr["VetoTau"] = makeVetoLeptonCorrectionHist(vetoTauHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Tau Pt CR", doDataOverMC=True, printdir=printdir)

    # #Make pt correction (signal region data/MC)
    # sfHistsSRPtCorr = sfHists.copy()
    # sfHistsSRPtCorr["VetoTau"] = makeVetoLeptonCorrectionHist(vetoTauHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var="lep1.Pt()", signifThreshold=1.0, regionName="Veto Tau Pt", doDataOverMC=False, histsToCorrect=multijetHistsForVetoTauCorrection, signalRegionVar="leadingGenLeptonPt", mtEfficiencyHist=mtTauPtHist, dPhiEfficiencyHist=dphiTauPtHist, printdir=printdir)

    # #for applying pt correction
    # auxSFsCRForVetoTauPtCorr = { 
    #     "VetoTau":("lep1.Pt()", "abs(lep1Type) == 15"), 
    #     }
    # auxSFsSRForVetoTauPtCorr = { 
    #     "VetoTau":("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 15"), 
    #     }

    # #apply the pt correction and examine the residual eta agreement (control region)
    # vetoTauHistsPtCorr = makeControlSampleHists("VetoTauControlRegionForEtaCorrection", 
    #            filenames=FILENAMES_VetoTau, samples=SAMPLES_VetoTau, 
    #            cutsMC=vetoTauControlRegionCutsMC, cutsData=vetoTauControlRegionCutsData, 
    #            bins=VetoTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, plotDensity=True, sfHists=sfHistsCRPtCorr, weightOpts=weightOpts, 
    #            auxSFs=auxSFsCRForVetoTauPtCorr, printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)

    # #apply the pt correction and examine the residual eta agreement (signal region)
    # multijetHistsPtCorrForVetoTauCorrection = makeControlSampleHists("MultiJetForVetoTauEtaCorrection", 
    #         filenames=FILENAMES_MULTIJET, samples=SAMPLES_MultiJet, 
    #         cutsMC=cutsMultiJetForVetoTau, cutsData=cutsMultiJetForVetoTau, 
    #         bins=MultiJetTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #         weightHists=weightHists, plotDensity=False, sfHists=sfHistsSRPtCorr, treeName="RazorInclusive", 
    #         auxSFs=auxSFsSRForVetoTauPtCorr, weightOpts=weightOpts, debugLevel=debugLevel, plotOpts=plotOpts, 
    #         printdir=printdir)

    # #Make residual eta correction (control region data/MC)
    # sfHistsCRPtEtaCorr = sfHists.copy()
    # sfHistsCRPtEtaCorr["VetoTauPt"] = sfHistsCRPtCorr["VetoTau"]
    # sfHistsCRPtEtaCorr["VetoTauEta"] = makeVetoLeptonCorrectionHist(vetoTauHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Tau Eta CR", doDataOverMC=True, printdir=printdir)

    # #Make residual eta correction (signal region data/MC)
    # makeVetoLeptonCorrectionHist(vetoTauHistsPtCorr, lumiData=LUMI_DATA, debugLevel=debugLevel, var="abs(lep1.Eta())", signifThreshold=1.0, regionName="Veto Tau Eta", doDataOverMC=False, histsToCorrect=multijetHistsPtCorrForVetoTauCorrection, signalRegionVar="abs(leadingGenLeptonEta)", mtEfficiencyHist=mtTauEtaHist, dPhiEfficiencyHist=dphiTauEtaHist, printdir=printdir)

    # auxSFsCRForVetoTauPtEtaCorr = { 
    #     "VetoTauPt":("lep1.Pt()", "abs(lep1Type) == 15"), 
    #     "VetoTauEta":("abs(lep1.Eta())", "abs(lep1Type) == 15"), 
    #     }

    # #Check control region after both corrections
    # vetoTauHistsPtEtaCorr = makeControlSampleHists("VetoTauControlRegionPtEtaCorrected", 
    #            filenames=FILENAMES_VetoTau, samples=SAMPLES_VetoTau, 
    #            cutsMC=vetoTauControlRegionCutsMC, cutsData=vetoTauControlRegionCutsData, 
    #            bins=VetoTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, plotDensity=False, sfHists=sfHistsCRPtEtaCorr, weightOpts=weightOpts, 
    #            auxSFs=auxSFsCRForVetoTauPtEtaCorr, printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)



    #########################################################
    #Veto Lepton 7Jet cross-check region
    #########################################################
    sfHists7Jet = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedTo7Jet.root", processNames=SAMPLES_VetoLepton, scaleFactorNames={ "ZInv":"GJetsInv" }, debugLevel=debugLevel)
    vetoLeptonHists = makeControlSampleHists("VetoLeptonControlRegion7Jet", 
               filenames=FILENAMES_VetoLepton, samples=SAMPLES_VetoLepton, 
               cutsMC=vetoLeptonControlRegion7JetCutsMC, cutsData=vetoLeptonControlRegion7JetCutsData, 
               bins=VetoLeptonControlRegion7JetBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
               weightHists=weightHists, plotDensity=False, sfHists=sfHists7Jet, weightOpts=weightOpts, 
               printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)


    ##########################################################
    ##Veto Tau 7Jet cross-check region
    ##########################################################
    vetoTauHists = makeControlSampleHists("VetoTauControlRegion7Jet", 
               filenames=FILENAMES_VetoTau, samples=SAMPLES_VetoTau, 
               cutsMC=vetoTauControlRegion7JetCutsMC, cutsData=vetoTauControlRegion7JetCutsData, 
               bins=VetoTauControlRegion7JetBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
               weightHists=weightHists, plotDensity=False, sfHists=sfHists7Jet, weightOpts=weightOpts, 
               printdir=printdir, plotOpts=plotOpts, debugLevel=debugLevel)


    ##########################################################
    #Z->LL dilepton control sample
    ##########################################################
    #sfHistsDileptonInv = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root", processNames=SAMPLES_DYJ2L_INV, scaleFactorNames=ScaleFactorNames_DYJ2L_INV, debugLevel=debugLevel)
    #dyjetsDileptonInvHists = makeControlSampleHists("DYJetsDileptonInv", 
    #            filenames=FILENAMES_2L_INV, samples=SAMPLES_DYJ2L_INV, 
    #            cutsMC=dyjetsDileptonInvCutsMC, cutsData=dyjetsDileptonInvCutsData, 
    #            bins=ZNuNu_2L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, plotDensity=True, sfHists=sfHistsDileptonInv, 
    #            weightOpts=weightOpts, sfVars=ScaleFactorVars_DYJ2L_INV,
    #            printdir=printdir, debugLevel=debugLevel)
    # dyjetsDileptonInv4JetHists = makeControlSampleHists("DYJetsDileptonInv4Jet", 
    #           filenames=FILENAMES_2L_INV, samples=SAMPLES_DYJ2L_INV, 
    #           cutsMC=dyjetsDileptonInv4JetCutsMC, cutsData=dyjetsDileptonInv4JetCutsData, 
    #           bins=ZNuNu_2L_4Jet_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, plotDensity=True, sfHists=sfHistsDileptonInv, 
    #           weightOpts=weightOpts, sfVars=ScaleFactorVars_DYJ2L_INV,
    #           printdir=printdir, debugLevel=debugLevel)

    ##Record discrepancies > 1 sigma
    #tmpSFHists = copy.copy(sfHistsDileptonInv)
    #del tmpSFHists["DYJetsInv"]
    #appendScaleFactors("DYJetsInv", dyjetsDileptonInvHists, tmpSFHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var=sfVarsDYJetsDileptonInv, signifThreshold=1.0, printdir=printdir)

    ##write DYJetsInv cross check scale factors
    #dyjetsDileptonInvOutfile = rt.TFile("RazorDYJetsDileptonInvCrossCheck.root", "RECREATE")
    #print "Writing histogram",tmpSFHists["DYJetsInv"].GetName(),"to file"
    #tmpSFHists["DYJetsInv"].Write("DYJetsDileptonInvCrossCheckScaleFactors")
    #dyjetsDileptonInvOutfile.Close()

