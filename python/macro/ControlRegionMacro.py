import sys
import argparse
import ROOT as rt

#local imports
import macro
from razorAnalysis import *

LUMI = 16 #in /pb
MCLUMI = 1 

SAMPLES_TTJ1L = ["SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["SingleTop", "TTJets", "WJets"]
SAMPLES_DYJ2L = ["VV", "SingleTop", "WJets", "TTJets", "DYJets"]

FILENAMES_1L = {
            "TTJets"   :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_razorskim.root",
            "WJets"    :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted_razorskim.root",
            "SingleTop":"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/SingleTop_1pb_weighted_razorskim.root",
            "Data"     :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/OneLeptonReduced_v2/SingleLepton_Run2015C_GOLDEN_razorskim_noduplicates.root"
            }
FILENAMES_2L = {
            "DYJets"   :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJetsToLL_M-50_HTBinned_1pb_weighted_razorskim.root",
            "TTJets"   :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_razorskim.root",
            "WJets"    :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_WJetsToLNu_HTBinned_1pb_weighted_razorskim.root",
            "SingleTop":"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleTop_1pb_weighted_razorskim.root",
            "VV"       :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_VV_1pb_weighted_razorskim.root",
            "Data"     :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p17/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleLepton_Run2015C_GoodLumi_razorskim_noduplicates.root",
            }

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="display detailed output messages",
                                action="store_true")
    args = parser.parse_args()
    debug = args.debug

    #initialize
    weightHists = loadWeightHists(debug)
    sfHists = {}

    #DYJets control sample
    dyjetsDileptonHists = makeControlSampleHists("DYJetsDilepton", filenames=FILENAMES_2L, samples=SAMPLES_DYJ2L, 
                cutsMC=dyjetsDileptonCutsMC, cutsData=dyjetsDileptonCutsData, bins=dyjetsDileptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, debug=debug)
    appendScaleFactors("DYJets", dyjetsDileptonHists, sfHists, debug=debug) 

    #TTJets control sample
    ttjetsSingleLeptonHists = makeControlSampleHists("TTJetsSingleLepton", filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                cutsMC=ttjetsSingleLeptonCutsMC, cutsData=ttjetsSingleLeptonCutsData, bins=ttjetsSingleLeptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, debug=debug)
    appendScaleFactors("TTJets", ttjetsSingleLeptonHists, sfHists, debug=debug)

    #WJets control sample
    wjetsSingleLeptonHists = makeControlSampleHists("WJetsSingleLepton", filenames=FILENAMES_1L, samples=SAMPLES_WJ1L, 
                cutsMC=wjetsSingleLeptonCutsMC, cutsData=wjetsSingleLeptonCutsData, bins=wjetsSingleLeptonBins,
                lumiMC=MCLUMI, lumiData=LUMI, weightHists=weightHists, sfHists=sfHists, debug=debug)
    appendScaleFactors("WJets", wjetsSingleLeptonHists, sfHists, debug=debug)

    #write scale factors
    outfile = rt.TFile("RazorScaleFactors.root", "RECREATE")
    for name in sfHists:
        print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
        sfHists[name].Write()
    outfile.Close()
