import sys, os
import argparse
import ROOT as rt

#local imports
from framework import Config
from macro import macro
from macro.razorAnalysis import OneLeptonForNJetsCutsData, OneLeptonForNJetsCutsMC, OneLeptonInvForNJetsCutsData, OneLeptonInvForNJetsCutsMC 
from macro.razorWeights import loadScaleFactorHists
from macro.razorMacros import makeControlSampleHists, appendScaleFactors
from ComputeScaleFactorsMacro import xbinsWJETS1L, colsWJETS1L, xbinsTTJETS1L, colsTTJETS1L, xbinsWJETS1LINV, colsWJETS1LINV
from SidebandMacro import LUMI as LUMI_DATA

MCLUMI = 1 

SAMPLES_TTJ1L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L_INV = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "TTJets", "WJetsInv"]
#apply correct version of MR when applying scale factors
ScaleFactorVars_WJETS1L_INV = { 
                              "WJetsInv":("MR_NoW","Rsq_NoW"),
                              "TTJets":("MR","Rsq"),
                              }

DIR_1L = "ControlRegions"
#DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonFull_1p23_2015Final/RazorSkim/"
PREFIX_1L = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
FILENAMES_1L = {
            "TTJets"   : DIR_1L+"/"+PREFIX_1L+"_TTJets_1pb_weighted_RazorSkim.root",
            "WJets"    : DIR_1L+"/"+PREFIX_1L+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_1L+"/"+PREFIX_1L+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_1L+"/"+PREFIX_1L+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_1L+"/"+PREFIX_1L+"_Other_1pb_weighted_RazorSkim.root",
            "ZInv"     : DIR_1L+"/"+PREFIX_1L+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
            "QCD"      : DIR_1L+"/"+PREFIX_1L+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",       
            "Data"     : DIR_1L+"/"+PREFIX_1L+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }

DIR_1L_INV = "ControlRegions"
#DIR_1L_INV = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonAddToMETFull_1p23_2015Final/RazorSkim/"
PREFIX_1L_INV = "RunTwoRazorControlRegions_OneLeptonAddToMetFull_SingleLeptonSkim"
FILENAMES_1L_INV = {
            "TTJets"   : DIR_1L_INV+"/"+PREFIX_1L_INV+"_TTJets_1pb_weighted_RazorSkim.root",
            "WJetsInv" : DIR_1L_INV+"/"+PREFIX_1L_INV+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_1L_INV+"/"+PREFIX_1L_INV+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_1L_INV+"/"+PREFIX_1L_INV+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_1L_INV+"/"+PREFIX_1L_INV+"_Other_1pb_weighted_RazorSkim.root",
            "ZInv"     : DIR_1L_INV+"/"+PREFIX_1L_INV+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
            "QCD"      : DIR_1L_INV+"/"+PREFIX_1L_INV+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_1L_INV+"/"+PREFIX_1L_INV+"_SingleLepton_Run2015D_RazorSkim_GoodLumiGolden_NoDuplicates.root"
            }

weightOpts = []

config = "config/run2_20151229_ControlRegion.config"
cfg = Config.Config(config)
binsMRLep = cfg.getBinning("WJetControlRegion")[0]
binsRsqLep = cfg.getBinning("WJetControlRegion")[1]
binsNBTags = [0.,1.,2.,3.,4.]
binsNJets80 = [0.,1.,2.,9.]
binsNJets = [1.,2.,4.,20.]
binsLepPt = [20.,25.,30.,35.,40.,45.,50.,70.,100]
ControlRegionBinning = { "MR":binsMRLep, "Rsq":binsRsqLep, "NJets40":binsNJets, ("MR","Rsq"):[]}
ControlRegionBinning_Inv = { "MR_NoW":binsMRLep, "Rsq_NoW":binsRsqLep, "NJets40_NoW":binsNJets, ("MR_NoW","Rsq_NoW"):[] }

printdir="ScaleFactorsNJetPlots"

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

    #make output directory
    os.system('mkdir -p '+printdir)

    sfVars = ("MR","Rsq")
    sfVars_NoW = ("MR_NoW", "Rsq_NoW")

    plotOpts = { "comment":False }

    ##########################################################
    #1L W+Jets & TTJets control sample
    ##########################################################    
    sfHists_OneLeptonForNJets = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root", processNames=SAMPLES_TTJ1L, debugLevel=debugLevel)
    OneLeptonForNJetsHists = makeControlSampleHists("OneLeptonForNJets", 
                  filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                  cutsMC=OneLeptonForNJetsCutsMC, cutsData=OneLeptonForNJetsCutsData, 
                  bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
                  weightHists=weightHists, sfHists=sfHists_OneLeptonForNJets, weightOpts=weightOpts, 
                  printdir=printdir, btags=-1, plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
                  plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
    appendScaleFactors("NJetsCorrection", OneLeptonForNJetsHists, sfHists_OneLeptonForNJets, var="NJets40", 
                  lumiData=LUMI_DATA, debugLevel=debugLevel, printdir=printdir)

    # ##########################################################
    # #W+Jets Add lepton to MET control sample
    # ##########################################################
    #sfHists_OneLeptonInvForNJets = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root", processNames=SAMPLES_WJ1L_INV, debugLevel=debugLevel)
    #OneLeptonInvForNJetsHists = makeControlSampleHists("OneLeptonInvForNJets", 
    #              filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
    #              cutsMC=OneLeptonInvForNJetsCutsMC, cutsData=OneLeptonInvForNJetsCutsData, 
    #              bins=ControlRegionBinning_Inv, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #              weightHists=weightHists, sfHists=sfHists_OneLeptonInvForNJets, weightOpts=weightOpts, 
    #              printdir=printdir, plotDensity=True, sfVars=ScaleFactorVars_WJETS1L_INV , debugLevel=debugLevel)
    #appendScaleFactors("NJetsNoWCorrection", OneLeptonInvForNJetsHists, sfHists_OneLeptonInvForNJets, 
    #              var="NJets40_NoW", lumiData=LUMI_DATA, debugLevel=debugLevel, printdir=printdir)

    #write scale factors
    outfile = rt.TFile("RazorNJetsScaleFactors.root", "RECREATE")
    print "Writing scale factor histogram",sfHists_OneLeptonForNJets["NJetsCorrection"].GetName(),"to file"
    sfHists_OneLeptonForNJets["NJetsCorrection"].Write(sfHists_OneLeptonForNJets["NJetsCorrection"].GetName())
    #print "Writing scale factor histogram",sfHists_OneLeptonInvForNJets["NJetsNoWCorrection"].GetName(),"to file"
    #sfHists_OneLeptonInvForNJets["NJetsNoWCorrection"].Write(sfHists_OneLeptonInvForNJets["NJetsNoWCorrection"].GetName())
    outfile.Close()
