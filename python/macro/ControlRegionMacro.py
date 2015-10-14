import sys
import argparse
import ROOT as rt

#local imports
import macro
from razorAnalysis import *

LUMI = 209 #in /pb
MCLUMI = 1 

SAMPLES_TTJ1L = ["SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["SingleTop", "TTJets", "WJets"]
SAMPLES_DYJ2L = ["VV", "SingleTop", "WJets", "TTJets", "DYJets"]

DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonFull_1p19/Tight30Skim"
PREFIX_1L = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
FILENAMES_1L = {
            "TTJets"   : DIR_1L+"/"+PREFIX_1L+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root_Tight30Skim.root_razorskim.root",
            "WJets"    : DIR_1L+"/"+PREFIX_1L+"_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root_Tight30Skim.root_razorskim.root",
            "SingleTop": DIR_1L+"/"+PREFIX_1L+"_SingleTop_1pb_weighted.root_Tight30Skim.root_razorskim.root",
            "Data"     : DIR_1L+"/"+PREFIX_1L+"_SingleLepton_Run2015D_GoodLumiGoldenUnblind.root_Tight30Skim_NoDuplicates.root_razorskim.root"
            }

#DIR_2L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17"
#FILENAMES_2L = {
#            "DYJets"   : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJetsToLL_M-50_HTBinned_1pb_weighted_razorskim.root",
#            "TTJets"   : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_razorskim.root",
#            "WJets"    : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_WJetsToLNu_HTBinned_1pb_weighted_razorskim.root",
#            "SingleTop": DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleTop_1pb_weighted_razorskim.root",
#            "VV"       : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_VV_1pb_weighted_razorskim.root",
#            "Data"     : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleLepton_Run2015C_GoodLumi_razorskim_noduplicates.root",
#            }

WEIGHTDIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors"
weightfilenames = {
        "muon": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D_Golden.root",
        "ele": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root",
        "muontrig": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015D_Golden.root",
        "eletrig": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015D_Golden.root",
        "pileup": WEIGHTDIR+"/PileupWeights/NVtxReweight_ZToMuMu_2015D_25ns_20150923.root",
        }
weighthistnames = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "muontrig": "ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight",
        "eletrig": "ScaleFactor_EleTriggerEleCombinedEffDenominatorTight",
        "pileup": "NVtxReweight",
        }

weightOpts = ["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"]

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
    weightHists = loadWeightHists(weightfilenames, weighthistnames, debugLevel)
    sfHists = {}

    #DYJets control sample
    #dyjetsDileptonHists = makeControlSampleHists("DYJetsDilepton", filenames=FILENAMES_2L, samples=SAMPLES_DYJ2L, 
    #            cutsMC=dyjetsDileptonCutsMC, cutsData=dyjetsDileptonCutsData, bins=dyjetsDileptonBins,
    #            lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, opts=weightOpts, debugLevel=debugLevel)
    #appendScaleFactors("DYJets", dyjetsDileptonHists, sfHists, debugLevel=debugLevel) 

    #TTJets control sample
    ttjetsSingleLeptonHists = makeControlSampleHists("TTJetsSingleLepton", filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                cutsMC=ttjetsSingleLeptonCutsMC, cutsData=ttjetsSingleLeptonCutsData, bins=ttjetsSingleLeptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, opts=weightOpts, debugLevel=debugLevel)
    appendScaleFactors("TTJets", ttjetsSingleLeptonHists, sfHists, debugLevel=debugLevel)

    #WJets control sample
    wjetsSingleLeptonHists = makeControlSampleHists("WJetsSingleLepton", filenames=FILENAMES_1L, samples=SAMPLES_WJ1L, 
                cutsMC=wjetsSingleLeptonCutsMC, cutsData=wjetsSingleLeptonCutsData, bins=wjetsSingleLeptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, opts=weightOpts, debugLevel=debugLevel)
    appendScaleFactors("WJets", wjetsSingleLeptonHists, sfHists, debugLevel=debugLevel)

    #write scale factors
    outfile = rt.TFile("RazorScaleFactors.root", "RECREATE")
    for name in sfHists:
        print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
        sfHists[name].Write()
    outfile.Close()
