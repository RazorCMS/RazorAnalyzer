import sys
import argparse
import ROOT as rt

#local imports
import macro
from razorAnalysis import *

LUMI_NONBLIND = 149 #in /pb
MCLUMI = 1 

SAMPLES_SIGNAL = ["TTV", "VV", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]

DIR_SIGNAL = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC"
PREFIX_SIGNAL = "RazorInclusive"
FILENAMES_SIGNAL = {
        "TTJets"    : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted_razorskim.root",
        "WJets"     : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_WJetsToLNu_HTBinned_1pb_weighted_razorskim.root",
        "SingleTop" : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_SingleTop_1pb_weighted_razorskim.root",
        "VV" : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_VV_1pb_weighted_razorskim.root",
        "TTV" : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_TTV_1pb_weighted_razorskim.root",
        "DYJets"     : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_DYJetsToLL_M-50_HTBinned_1pb_weighted_razorskim.root",
        "ZInv"     : DIR_SIGNAL+"/"+PREFIX_SIGNAL+"_ZJetsToNuNu_HTBinned_1pb_weighted_razorskim.root",
        "Data"      : DIR_SIGNAL+'/'+PREFIX_SIGNAL+'_CombinedLeptonic_Run2015D_GoodLumiGolden_razorskim.root',
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

    #get scale factor histograms
    sfHists = loadScaleFactorHists(processNames=SAMPLES_SIGNAL, debugLevel=debugLevel)

    #estimate yields in signal region
    for lepType in ["Mu", "Ele"]:
        for jets in ["MultiJet"]:
            boxName = lepType+jets
            for btags in [0,1,2,3]:
                print "---",boxName,"Box,",btags,"B-tags ---"
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
                        weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", opts=[], 
                        debugLevel=debugLevel)
