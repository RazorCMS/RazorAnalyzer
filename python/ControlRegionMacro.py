import sys
import argparse
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import *
from macro.razorWeights import *
from macro.razorMacros import *

LUMI_FULL = 1264 #in /pb
MCLUMI = 1 

SAMPLES_TTJ1L = ["DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["DYJets", "SingleTop", "TTJets", "WJets"]

DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonFull_1p20"
PREFIX_1L = "RazorControlRegions"
FILENAMES_1L = {
            "TTJets"   : DIR_1L+"/"+PREFIX_1L+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted_razorskim.root",
            "WJets"    : DIR_1L+"/"+PREFIX_1L+"_WJetsToLNu_HTBinned_1pb_weighted_razorskim.root",
            "SingleTop": DIR_1L+"/"+PREFIX_1L+"_SingleTop_1pb_weighted_razorskim.root",
            "DYJets"   : DIR_1L+"/"+PREFIX_1L+"_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted_razorskim.root",
            "Data"     : DIR_1L+"/Tight30Skim/"+PREFIX_1L+"_SingleLepton_Run2015D_Tight30Skim_GoodLumiGolden_NoDuplicates_razorskim.root"
            }

WEIGHTDIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors"
LEPTONWEIGHTDIR = "LeptonEfficiencies/20151013_PR_2015D_Golden_1264"
weightfilenames = {
        "muon": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D_Golden.root",
        "ele": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root",
        "muontrig": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015D_Golden.root",
        "eletrig": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015D_Golden.root",
        "pileup": WEIGHTDIR+"/PileupWeights/NVtxReweight_ZToMuMu_2015D_1264ipb.root",
        }
weighthistnames = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "muontrig": "ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight",
        "eletrig": "ScaleFactor_EleTriggerEleCombinedEffDenominatorTight",
        "pileup": "NVtxReweight",
        }
weightOpts = ["doNPVWeights", "doLep1Weights", "do1LepTrigWeights"]

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
    #NOTE: applying scale factors as a cross check.
    sfHists = loadScaleFactorHists(processNames=SAMPLES_WJ1L, debugLevel=debugLevel)

    #WJets control sample
    wjetsSingleLeptonHists = makeControlSampleHists("WJetsSingleLepton", 
                filenames=FILENAMES_1L, samples=SAMPLES_WJ1L, 
                cutsMC=wjetsSingleLeptonCutsMC, cutsData=wjetsSingleLeptonCutsData, 
                bins=wjetsSingleLeptonBinsReduced, lumiMC=MCLUMI, lumiData=LUMI_FULL, 
                weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, debugLevel=debugLevel)
    #appendScaleFactors("WJets", wjetsSingleLeptonHists, sfHists, lumiData=LUMI_FULL, debugLevel=debugLevel)

    #TTJets control sample
    ttjetsSingleLeptonHists = makeControlSampleHists("TTJetsSingleLepton", 
                filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                cutsMC=ttjetsSingleLeptonCutsMC, cutsData=ttjetsSingleLeptonCutsData, 
                bins=ttjetsSingleLeptonBinsReduced, lumiMC=MCLUMI, lumiData=LUMI_FULL, 
                weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, debugLevel=debugLevel)
    #appendScaleFactors("TTJets", ttjetsSingleLeptonHists, sfHists, lumiData=LUMI_FULL, debugLevel=debugLevel)

    #write scale factors
    #outfile = rt.TFile("RazorScaleFactors.root", "RECREATE")
    #for name in sfHists:
    #    print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
    #    sfHists[name].Write()
    #outfile.Close()
