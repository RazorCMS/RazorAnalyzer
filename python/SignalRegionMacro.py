import sys
import argparse
import ROOT as rt

#local imports
import macro.macro as macro
from macro.razorAnalysis import *
from macro.razorWeights import *
from macro.razorMacros import *

LUMI_NONBLIND = 133 #in /pb
MCLUMI = 1 

SAMPLES_SIGNAL = ["TTV", "VV", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]

DIR_MC= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC_FullRazorInclusive/combined"
DIR_DATA= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data"
PREFIX_SIGNAL = "RazorInclusive"
FILENAMES_SIGNAL = {
        "TTJets"    : DIR_MC+"/"+PREFIX_SIGNAL+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root",
        "WJets"     : DIR_MC+"/"+PREFIX_SIGNAL+"_WJets_1pb_weighted.root",
        "SingleTop" : DIR_MC+"/"+PREFIX_SIGNAL+"_SingleTop_1pb_weighted.root",
        "VV" : DIR_MC+"/"+PREFIX_SIGNAL+"_VV_1pb_weighted.root",
        "TTV" : DIR_MC+"/"+PREFIX_SIGNAL+"_TTV_1pb_weighted.root",
        "DYJets"     : DIR_MC+"/"+PREFIX_SIGNAL+"_DYJetsToLL_1pb_weighted.root",
        "ZInv"     : DIR_MC+"/"+PREFIX_SIGNAL+"_ZJetsToNuNu_1pb_weighted.root",
        "Data"      : DIR_DATA+'/'+PREFIX_SIGNAL+'_SingleLepton_Run2015D_GoodLumiUnblind_NoDuplicates_razorskim.root',
        }

WEIGHTDIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors"
weightfilenames = {
        "muon": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D_Golden.root",
        "ele": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root",
        "muontrig": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015D_Golden.root",
        "eletrig": WEIGHTDIR+"/LeptonEfficiencies/20151013_PR_2015D_GoldenUnblind/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015D_Golden.root",
        "pileup": WEIGHTDIR+"/PileupWeights/NVtxReweight_ZToMuMu_2015D_unblind.root",
        }
weighthistnames = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "muontrig": "ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight",
        "eletrig": "ScaleFactor_EleTriggerEleCombinedEffDenominatorTight",
        "pileup": "NVtxReweight",
        }
weightOpts = ["doNVtxWeights"]

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
    #sfHists = {}
    sfHists = loadScaleFactorHists(processNames=SAMPLES_SIGNAL, debugLevel=debugLevel)

    #estimate yields in signal region
    for lepType in ["Mu", "Ele"]:
        for jets in ["MultiJet"]:
            boxName = lepType+jets
            btaglist = [0]
            #btaglist = [0,1,2]
            for btags in btaglist:
                print "\n---",boxName,"Box,",btags,"B-tags ---"
                #get correct cuts string
                thisBoxCuts = razorCuts[boxName]
                if btags < len(btaglist)-1:
                    thisBoxCuts += " && nBTaggedJets == "+str(btags)
                else:
                    thisBoxCuts += " && nBTaggedJets >= "+str(btags)

                if len(btaglist) > 1:
                    extboxName = boxName+str(btags)+"BTag"
                else:
                    extboxName = boxName
                makeControlSampleHists(extboxName, 
                        filenames=FILENAMES_SIGNAL, samples=SAMPLES_SIGNAL, 
                        cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                        bins=leptonicSignalRegionBins, lumiMC=MCLUMI, lumiData=LUMI_NONBLIND, 
                        weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                        weightOpts=weightOpts, shapeErrors=["muoneff", "eleeff", "jes"], debugLevel=debugLevel)
