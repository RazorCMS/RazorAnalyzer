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

DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2"
FILENAMES_1L = {
            "TTJets"   : DIR_1L+"/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_razorskim.root",
            "WJets"    : DIR_1L+"/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted_razorskim.root",
            "SingleTop": DIR_1L+"/SingleTop_1pb_weighted_razorskim.root",
            "Data"     :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/OneLeptonReduced_v2/SingleLepton_Run2015C_GOLDEN_razorskim_noduplicates.root"
            }

DIR_2L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17"
FILENAMES_2L = {
            "DYJets"   : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJetsToLL_M-50_HTBinned_1pb_weighted_razorskim.root",
            "TTJets"   : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_razorskim.root",
            "WJets"    : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_WJetsToLNu_HTBinned_1pb_weighted_razorskim.root",
            "SingleTop": DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleTop_1pb_weighted_razorskim.root",
            "VV"       : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_VV_1pb_weighted_razorskim.root",
            "Data"     : DIR_2L+"/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleLepton_Run2015C_GoodLumi_razorskim_noduplicates.root",
            }

weightfilenames = {
        "muon": "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20150924_PR_2015D/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D.root",
        "ele": "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20150924_PR_2015D/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D.root",
        "pileup": "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/NVtxReweight_ZToMuMu_2015D_25ns_20150923.root",
        }
weighthistnames = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "pileup": "NVtxReweight",
        }

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
    dyjetsDileptonHists = makeControlSampleHists("DYJetsDilepton", filenames=FILENAMES_2L, samples=SAMPLES_DYJ2L, 
                cutsMC=dyjetsDileptonCutsMC, cutsData=dyjetsDileptonCutsData, bins=dyjetsDileptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, debugLevel=debugLevel)
    appendScaleFactors("DYJets", dyjetsDileptonHists, sfHists, debugLevel=debugLevel) 

    #TTJets control sample
    ttjetsSingleLeptonHists = makeControlSampleHists("TTJetsSingleLepton", filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                cutsMC=ttjetsSingleLeptonCutsMC, cutsData=ttjetsSingleLeptonCutsData, bins=ttjetsSingleLeptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, debugLevel=debugLevel)
    appendScaleFactors("TTJets", ttjetsSingleLeptonHists, sfHists, debugLevel=debugLevel)

    #WJets control sample
    wjetsSingleLeptonHists = makeControlSampleHists("WJetsSingleLepton", filenames=FILENAMES_1L, samples=SAMPLES_WJ1L, 
                cutsMC=wjetsSingleLeptonCutsMC, cutsData=wjetsSingleLeptonCutsData, bins=wjetsSingleLeptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, debugLevel=debugLevel)
    appendScaleFactors("WJets", wjetsSingleLeptonHists, sfHists, debugLevel=debugLevel)

    #write scale factors
    outfile = rt.TFile("RazorScaleFactors.root", "RECREATE")
    for name in sfHists:
        print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
        sfHists[name].Write()
    outfile.Close()
