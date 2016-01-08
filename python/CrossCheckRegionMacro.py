import sys
import argparse
import copy
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import wjetsSingleLeptonCutsMC, wjetsSingleLeptonCutsData, ttjetsSingleLeptonCutsMC, ttjetsSingleLeptonCutsData
from macro.razorWeights import *
from macro.razorMacros import *

LUMI_DATA = 2185 #in /pb
MCLUMI = 1 

SAMPLES_TTJ1L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["Other", "DYJets", "SingleTop", "TTJets", "WJets"]
SAMPLES_TTJ2L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_DYJ2L = ["Other", "SingleTop", "WJets", "TTJets", "DYJets"]
SAMPLES_VetoLepton = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_VetoTau = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_DYJ2L_INV = ["Other", "SingleTop", "WJets", "TTJets", "DYJets"]
ScaleFactorNames_DYJ2L_INV = {"Other"     : "Other", 
                              "SingleTop" : "SingleTop", 
                              "WJets"     : "WJets", 
                              "TTJets"    :"TTJets" , 
                              "DYJets"    :"WJetsInv"
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
            "DYJets"   : DIR_2L_INV+"/"+PREFIX_2L_INV+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_2L_INV+"/"+PREFIX_2L_INV+"_Other_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_2L_INV+"/"+PREFIX_2L_INV+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }




weightOpts = []

config = "config/run2_20151229_ControlRegion.config"
cfg = Config.Config(config)
VetoLeptonBinsMRLep = cfg.getBinning("VetoLeptonControlRegion")[0]
VetoLeptonBinsRsqLep = cfg.getBinning("VetoLeptonControlRegion")[1]
VetoLeptonBinsLepPt = [5, 10, 15, 20.,30.,40.,100]
VetoLeptonBinsLepEta = [0, 0.5, 1.0, 1.5, 2.0, 2.5]
VetoLeptonControlRegionBinning = { "MR":VetoLeptonBinsMRLep, "Rsq":VetoLeptonBinsRsqLep, "lep1.Pt()":VetoLeptonBinsLepPt , "abs(lep1.Eta())":VetoLeptonBinsLepEta, ("MR","Rsq"):[], ("abs(lep1.Eta())","lep1.Pt()"):[]}
VetoTauControlRegionBinning = { "MR":VetoLeptonBinsMRLep, "Rsq":VetoLeptonBinsRsqLep, ("MR","Rsq"):[]}
TTJetsDileptonBinsMRLep = cfg.getBinning("TTJetsDileptonControlRegion")[0]
TTJetsDileptonBinsRsqLep = cfg.getBinning("TTJetsDileptonControlRegion")[1]
TTJetsDileptonBinsNBTags = [0.,1.,2.,3.,4.]
TTJetsDileptonBinsNJets80 = [0.,1.,2.,3.,4.]
TTJetsDileptonBinsNJets = [0.,1.,2.,3.,4.,5.,6.,7.]
TTJetsDileptonBinsLepPt = [20.,25.,30.,35.,40.,45.,50.,70.,100]
TTJetsDileptonControlRegionBinning = { "MR":TTJetsDileptonBinsMRLep, "Rsq":TTJetsDileptonBinsRsqLep, "NBJetsMedium":TTJetsDileptonBinsNBTags, "NJets80":TTJetsDileptonBinsNJets80, "NJets40":TTJetsDileptonBinsNJets, ("MR","Rsq"):[]}
ZNuNu_2L_BinsMRLep = cfg.getBinning("WJetControlRegion")[0]
ZNuNu_2L_BinsRsqLep = cfg.getBinning("WJetControlRegion")[1]
ZNuNu_2L_BinsNBTags = [0,1,2,3,4]
ZNuNu_2L_BinsNJets80 = [0,1,2,3,4]
ZNuNu_2L_BinsNJets = [0,1,2,3,4,5,6,7,8]
ZNuNu_2L_ControlRegionBinning = { "MR_NoZ":ZNuNu_2L_BinsMRLep, "Rsq_NoZ":ZNuNu_2L_BinsRsqLep, "NBJetsMedium":ZNuNu_2L_BinsNBTags, "NJets80":ZNuNu_2L_BinsNJets80, "NJets40":ZNuNu_2L_BinsNJets, ("MR_NoZ","Rsq_NoZ"):[] }


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

    #make output directory
    os.system('mkdir -p '+printdir)

    sfHists = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorScaleFactors.root", processNames=SAMPLES_VetoLepton, debugLevel=debugLevel)
    sfVars = ("MR","Rsq")


    ##########################################################
    #TTJets dilepton control sample
    ##########################################################
    #ttjetsDileptonHists = makeControlSampleHists("TTJetsDilepton", 
    #            filenames=FILENAMES_2L, samples=SAMPLES_TTJ2L, 
    #            cutsMC=ttjetsDileptonCutsMC, cutsData=ttjetsDileptonCutsData, 
    #            bins=TTJetsDileptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, 
    #            printdir=printdir, debugLevel=debugLevel)

    ##Record discrepancies > 1 sigma
    #tmpSFHists = copy.copy(sfHists)
    #del tmpSFHists["TTJets"]
    #appendScaleFactors("TTJets", ttjetsDileptonHists, tmpSFHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var=sfVars, signifThreshold=1.0, printdir=printdir)

    ##write TTJets cross check scale factors
    #ttjetsDileptonOutfile = rt.TFile("RazorTTJetsDileptonCrossCheck.root", "RECREATE")
    #print "Writing histogram",tmpSFHists["TTJets"].GetName(),"to file"
    #tmpSFHists["TTJets"].Write("TTJetsDileptonCrossCheckScaleFactors")
    #ttjetsDileptonOutfile.Close()

    ##########################################################
    #Veto Lepton cross-check region
    ##########################################################
    vetoLeptonHists = makeControlSampleHists("VetoLeptonControlRegion", 
                filenames=FILENAMES_VetoLepton, samples=SAMPLES_VetoLepton, 
                cutsMC=vetoLeptonControlRegionCutsMC, cutsData=vetoLeptonControlRegionCutsData, 
                bins=VetoLeptonControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
                weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, 
                printdir=printdir, debugLevel=debugLevel)

    #Record discrepancies > 1 sigma
    vetoLeptonSFs = makeVetoLeptonCorrectionHist(vetoLeptonHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var=sfVars, signifThreshold=1.0, regionName="Veto Lepton", printdir=printdir)

    #write veto lepton scale factors
    vetoLeptonOutfile = rt.TFile("RazorVetoLeptonCrossCheck.root", "RECREATE")
    print "Writing histogram",vetoLeptonSFs.GetName(),"to file"
    vetoLeptonSFs.Write("VetoLeptonScaleFactors")
    vetoLeptonOutfile.Close()

    ##########################################################
    ##Veto Tau cross-check region
    ##########################################################
    #vetoTauHists = makeControlSampleHists("VetoTauControlRegion", 
    #             filenames=FILENAMES_VetoTau, samples=SAMPLES_VetoTau, 
    #             cutsMC=vetoTauControlRegionCutsMC, cutsData=vetoTauControlRegionCutsData, 
    #             bins=VetoTauControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #             weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, 
    #             printdir=printdir, debugLevel=debugLevel)

    ##Record discrepancies > 1 sigma
    #vetoTauSFs = makeVetoLeptonCorrectionHist(vetoTauHists, lumiData=LUMI_DATA, debugLevel=debugLevel, var=sfVars, signifThreshold=1.0, regionName="Veto Tau", printdir=printdir)

    ##write veto tau correction factors
    #vetoTauOutfile = rt.TFile("RazorVetoTauCrossCheck.root", "RECREATE")
    #print "Writing histogram",vetoTauSFs.GetName(),"to file"
    #vetoTauSFs.Write("VetoTauScaleFactors")
    #vetoTauOutfile.Close()



    # ##########################################################
    # #Z->LL dilepton control sample
    # ##########################################################
    # sfHistsDileptonInv = loadScaleFactorHists(sfFilename="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/ScaleFactors/RazorScaleFactors.root", processNames=SAMPLES_DYJ2L_INV, scaleFactorNames=ScaleFactorNames_DYJ2L_INV, debugLevel=debugLevel)

    # dyjetsDileptonInvHists = makeControlSampleHists("DYJetsDileptonInv", 
    #             filenames=FILENAMES_2L_INV, samples=SAMPLES_DYJ2L_INV, 
    #             cutsMC=dyjetsDileptonInvCutsMC, cutsData=dyjetsDileptonInvCutsData, 
    #             bins=ZNuNu_2L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #             weightHists=weightHists, sfHists=sfHistsDileptonInv, weightOpts=weightOpts, 
    #             printdir=printdir, debugLevel=debugLevel)
 
