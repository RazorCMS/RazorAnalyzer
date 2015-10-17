import sys
import argparse
import ROOT as rt

#local imports
import macro
from razorAnalysis import *

LUMI_NONBLIND = 210 #in /pb
LUMI_FULL = 378 #in /pb
MCLUMI = 1 

SAMPLES_TTJ1L = ["SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["SingleTop", "TTJets", "WJets"]
SAMPLES_DYJ2L = ["VV", "SingleTop", "WJets", "TTJets", "DYJets"]
SAMPLES_SIGNAL = ["SingleTop", "WJets", "TTJets"]

DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonFull_1p19/Tight30RazorSkim"
PREFIX_1L = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
FILENAMES_1L = {
            "TTJets"   : DIR_1L+"/"+PREFIX_1L+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root_tight30razorskim.root",
            "WJets"    : DIR_1L+"/"+PREFIX_1L+"_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root_tight30razorskim.root",
            "SingleTop": DIR_1L+"/"+PREFIX_1L+"_SingleTop_1pb_weighted_tight30razorskim.root",
            "Data"     : DIR_1L+"/"+PREFIX_1L+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates.root_tight30razorskim.root"
            }

DIR_SIGNAL = "RazorInclusiveSpring15Backgrounds"
PREFIX_SIGNAL = "RazorInclusive"
FILENAMES_SIGNAL = {
        "TTJets"    : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_TTJets_1pb_weighted.root",
        "WJets"     : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop" : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_SingleTop_1pb_weighted.root",
        "Data"      : 'root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19/RazorInclusive_SingleLepton_Run2015D_GoodLumiUnblind_NoDuplicates.root_razorskim.root',
        }

WEIGHTDIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors"
weightfilenames = {
        "muon": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D_Golden.root",
        "ele": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root",
        "muontrig": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015D_Golden.root",
        "eletrig": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015D_Golden.root",
        "pileup": WEIGHTDIR+"/PileupWeights/NVtxReweight_ZToMuMu_2015Dv3_378ipb.root",
        }
weighthistnames = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "muontrig": "ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight",
        "eletrig": "ScaleFactor_EleTriggerEleCombinedEffDenominatorTight",
        "pileup": "NVtxReweight",
        }

weightOpts = ["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"]

ttjetsSingleLeptonBinsReduced = {
        "MR" : [300, 400, 550, 700, 4000],
        "Rsq": [0.15,0.25,0.41, 1.5]
        }
#ttjetsSingleLeptonBinsReduced = {
#        "MR" : [300, 350, 400, 450, 500, 550, 700, 900, 4000],
#        "Rsq": [0.15,0.175,0.20,0.225,0.25,0.30,0.41, 1.5]
#        }

wjetsSingleLeptonBinsReduced = {
        "MR" : [300, 400, 550, 700, 4000],
        "Rsq": [0.15,0.25,0.41, 1.5]
        }

#wjetsSingleLeptonBinsReduced = {
#        "MR" : [300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000],
#        "Rsq": [0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5]
#        }

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
    ttjetsSingleLeptonHists = makeControlSampleHists("TTJetsSingleLepton", 
                filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                cutsMC=ttjetsSingleLeptonCutsMC, cutsData=ttjetsSingleLeptonCutsData, 
                bins=ttjetsSingleLeptonBinsReduced, lumiMC=MCLUMI, lumiData=LUMI_FULL, 
                weightHists=weightHists, sfHists=sfHists, opts=weightOpts, debugLevel=debugLevel)
    appendScaleFactors("TTJets", ttjetsSingleLeptonHists, sfHists, debugLevel=debugLevel)

    #WJets control sample
    wjetsSingleLeptonHists = makeControlSampleHists("WJetsSingleLepton", 
                filenames=FILENAMES_1L, samples=SAMPLES_WJ1L, 
                cutsMC=wjetsSingleLeptonCutsMC, cutsData=wjetsSingleLeptonCutsData, 
                bins=wjetsSingleLeptonBinsReduced, lumiMC=MCLUMI, lumiData=LUMI_FULL, 
                weightHists=weightHists, sfHists=sfHists, opts=weightOpts, debugLevel=debugLevel)
    appendScaleFactors("WJets", wjetsSingleLeptonHists, sfHists, debugLevel=debugLevel)

    #write scale factors
    outfile = rt.TFile("RazorScaleFactors.root", "RECREATE")
    for name in sfHists:
        print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
        sfHists[name].Write()
    outfile.Close()

    #estimate yields in signal region
    for lepType in ["Mu", "Ele"]:
        for jets in ["MultiJet"]:
            boxName = lepType+jets
            for btags in [0,1,2,3]:
                #get correct cuts string
                thisBoxCutsData = razorCutsData[boxName]
                thisBoxCutsMC = razorCutsMC[boxName]
                if btags < 3:
                    thisBoxCutsData += " && nBTaggedJets == "+str(btags)
                    thisBoxCutsMC += " && nBTaggedJets == "+str(btags)
                else:
                    thisBoxCutsData += " && nBTaggedJets >= "+str(btags)
                    thisBoxCutsMC += " && nBTaggedJets >= "+str(btags)

                makeControlSampleHists((boxName+str(btags)+"BTag"), 
                        filenames=FILENAMES_SIGNAL, samples=SAMPLES_SIGNAL, 
                        cutsMC=thisBoxCutsMC, cutsData=thisBoxCutsData, 
                        bins=leptonicSignalRegionBins, lumiMC=MCLUMI, lumiData=LUMI_NONBLIND, 
                        weightHists=weightHists, sfHists=sfHists, opts=[], 
                        debugLevel=debugLevel)
